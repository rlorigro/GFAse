from torch.nn import Sequential, Linear, SELU
from torch_geometric.loader import DataLoader, NeighborLoader
from torch_geometric.nn import MetaLayer
from torch_geometric.data import Data
from torch_scatter import scatter_mean
from matplotlib import pyplot
import torch


def parse_csv_as_torch_data(nodes_path, edges_path):
    node_data = list()
    node_truth = list()
    id_to_name = dict()
    name_to_id = dict()

    edge_data = list()
    edge_index = [list(), list()]

    with open(nodes_path, 'r') as file:
        for l,line in enumerate(file):

            # Skip header
            if l == 0:
                continue

            id,name,length,contact_coverage,hash_coverage = line.strip().split(',')

            is_bubble = False

            suffix = name[-2:]
            if (not name.startswith("UR")) and (suffix == ".0" or suffix == ".1"):
                is_bubble = True

            d = torch.zeros(4, dtype=torch.float32)
            d[0] = float(id)
            d[1] = float(length)
            d[2] = float(contact_coverage)
            d[3] = float(hash_coverage)

            y_i = torch.zeros(2)
            y_i[0] = int(id)
            y_i[1] = int(is_bubble)

            node_data.append(d)
            node_truth.append(y_i)
            id_to_name[int(id)] = name
            name_to_id[name] = int(id)

    node_data = sorted(node_data, key=lambda x: x[0])
    node_data = torch.stack(node_data)[:,1:]

    node_truth = sorted(node_truth, key=lambda x: x[0])
    node_truth = torch.stack(node_truth)[:,1:]

    with open(edges_path, 'r') as file:
        for l,line in enumerate(file):
            # skip header
            if l == 0:
                continue

            id_a,id_b,name_a,name_b,contact_weight,hash_weight = line.strip().split(',')

            edge_index[0].append(int(id_a))
            edge_index[1].append(int(id_b))
            edge_index[0].append(int(id_b))
            edge_index[1].append(int(id_a))

            d = torch.zeros(2, dtype=torch.float32)
            d[0] = float(contact_weight)
            d[1] = float(hash_weight)

            # Need to duplicate attribute for both directions of edge
            edge_data.append(d)
            edge_data.append(d)

    edge_index[0] = torch.LongTensor(edge_index[0])
    edge_index[1] = torch.LongTensor(edge_index[1])

    edge_index = torch.stack(edge_index)
    edge_data = torch.stack(edge_data)

    print("%s:\t[%d,%d]" % ("node_truth", node_truth.shape[0], node_truth.shape[1]))
    print("%s:\t[%d,%d]" % ("node_data", node_data.shape[0], node_data.shape[1]))
    print("%s:\t[%d,%d]" % ("edge_index", edge_index.shape[0], edge_index.shape[1]))
    print("%s:\t[%d,%d]" % ("edge_data", edge_data.shape[0], edge_data.shape[1]))

    return Data(x=node_data, y=node_truth, edge_index=edge_index, edge_attr=edge_data), name_to_id, id_to_name


class EdgeModel(torch.nn.Module):
    def __init__(self, node_size, edge_size, global_size, hidden_size=32):
        super().__init__()

        self.input_size = 2*node_size + edge_size + global_size

        self.edge_mlp = Sequential(
            Linear(self.input_size, hidden_size),
            SELU(),
            Linear(hidden_size, edge_size)
        )

    def forward(self, src, dest, edge_attr, u, batch):
        # src, dest: [E, F_x], where E is the number of edges.
        # edge_attr: [E, F_e]
        # u: [B, F_u], where B is the number of graphs.
        # batch: [E] with max entry B - 1.

        # print("--- forward ---")
        # print("src", src.shape, src.type())
        # print("dest", dest.shape, dest.type())
        # print("edge_attr", edge_attr.shape, edge_attr.type())
        # print("u", u.shape, u.type())
        # print("batch", batch.shape, batch.type())
        # print("u[batch]", u[batch].shape)

        out = torch.cat([src, dest, edge_attr, u[batch]], 1)
        return self.edge_mlp(out)


class NodeModel(torch.nn.Module):
    def __init__(self, node_size, edge_size, global_size, hidden_size=64):
        super().__init__()

        self.input_size_1 = node_size + edge_size

        self.node_mlp_1 = Sequential(
            Linear(self.input_size_1, hidden_size),
            SELU(),
            Linear(hidden_size, hidden_size)
        )

        self.input_size_2 = node_size + hidden_size + global_size

        self.node_mlp_2 = Sequential(
            Linear(self.input_size_2, hidden_size),
            SELU(),
            Linear(hidden_size, node_size)
        )

    def forward(self, x, edge_index, edge_attr, u, batch):
        # x: [N, F_x], where N is the number of nodes.
        # edge_index: [2, E] with max entry N - 1.
        # edge_attr: [E, F_e]
        # u: [B, F_u]
        # batch: [N] with max entry B - 1.

        row, col = edge_index
        out = torch.cat([x[row], edge_attr], dim=1)
        out = self.node_mlp_1(out)

        out = scatter_mean(out, col, dim=0, dim_size=x.size(0))
        out = torch.cat([x, out, u[batch]], dim=1)
        return self.node_mlp_2(out)


class GlobalModel(torch.nn.Module):
    def __init__(self, node_size, global_size, hidden_size=64):
        super().__init__()

        self.input_size = node_size + global_size

        self.global_mlp = Sequential(
            Linear(self.input_size, hidden_size),
            SELU(),
            Linear(hidden_size, global_size)
        )

    def forward(self, x, edge_index, edge_attr, u, batch):
        # x: [N, F_x], where N is the number of nodes.
        # edge_index: [2, E] with max entry N - 1.
        # edge_attr: [E, F_e]
        # u: [B, F_u]
        # batch: [N] with max entry B - 1.
        out = torch.cat([u, scatter_mean(x, batch, dim=0)], dim=1)
        return self.global_mlp(out)


class NodeClassifier(torch.nn.Module):
    '''
    A simple, general purpose, fully connected network
    '''
    def __init__(self, node_size):
        # Perform initialization of the pytorch superclass
        super(NodeClassifier, self).__init__()

        # Define output layer dimensions
        D_in, H1, H2, D_out = [node_size, 64, 64, 1]    # These numbers correspond to each layer: [input, hidden_1, output]

        # Define output layer types
        self.norm0 = torch.nn.BatchNorm1d(D_in)

        self.linear1 = Linear(D_in, H1)
        self.s1 = SELU()
        self.linear2 = Linear(H1, H2)
        self.s2 = SELU()
        self.linear3 = Linear(H2, D_out)

    def forward(self, x):
        '''
        This method defines the network layering and activation functions
        '''

        x = self.linear1(x) # hidden layer
        x = self.s1(x)

        x = self.linear2(x) # hidden layer
        x = self.s2(x)

        x = self.linear3(x) # hidden layer

        return x


class Graphnet(torch.nn.Module):
    '''
    A simple, general purpose graphnet which classifies nodes after creating embeddings
    '''
    def __init__(self, node_size, edge_size, global_size):
        # Perform initialization of the pytorch superclass
        super(Graphnet, self).__init__()

        self.g1 = MetaLayer(
            edge_model=EdgeModel(node_size, edge_size, global_size),
            node_model=NodeModel(node_size, edge_size, global_size),
            global_model=GlobalModel(node_size, global_size)
        )

        self.g2 = MetaLayer(
            edge_model=EdgeModel(node_size, edge_size, global_size),
            node_model=NodeModel(node_size, edge_size, global_size),
            global_model=GlobalModel(node_size, global_size)
        )

        self.g3 = MetaLayer(
            edge_model=EdgeModel(node_size, edge_size, global_size),
            node_model=NodeModel(node_size, edge_size, global_size),
            global_model=GlobalModel(node_size, global_size)
        )

        self.output_model = NodeClassifier(node_size)

    def forward(self, data):
        '''
        This method performs graph level aggregation of features and then uses node features to output a classification
        '''

        u = torch.mean(data.x, dim=0).type(torch.FloatTensor).unsqueeze(0)

        x, edge_attr, u = self.g1(x=data.x, edge_index=data.edge_index, edge_attr=data.edge_attr, u=u, batch=data.batch)
        x, edge_attr, u = self.g2(x=data.x, edge_index=data.edge_index, edge_attr=data.edge_attr, u=u, batch=data.batch)
        x, edge_attr, u = self.g3(x=data.x, edge_index=data.edge_index, edge_attr=data.edge_attr, u=u, batch=data.batch)
        x = self.output_model(x)

        return x


def train_batch(model, data, optimizer, loss_fn):
    # Run forward calculation
    y_predict = model.forward(data)

    # print()
    # print(x)
    # print(torch.sigmoid(y_predict).data.numpy().T[:10])
    # print(y[:10])

    # Compute loss.
    loss = loss_fn(y_predict, data.y.float())

    # Before the backward pass, use the optimizer object to zero all of the
    # gradients for the variables it will update (which are the learnable weights
    # of the model)
    optimizer.zero_grad()

    # Backward pass: compute gradient of the loss with respect to model
    # parameters
    loss.backward()

    # Calling the step function on an Optimizer makes an update to its
    # parameters
    optimizer.step()

    return loss.data.item()


def train(model, loader, optimizer, loss_fn, epochs=5):
    model.train()

    losses = list()

    batch_index = 0
    for e in range(epochs):
        for data in loader:

            if data.batch is None:
                data.batch = torch.LongTensor([0]*data.x.shape[0])

            # print("u", u.shape)
            print("x", data.x.shape)
            print("y", data.y.shape)
            print("edge_index", data.edge_index.shape)
            print("edge_attr", data.edge_attr.shape)
            print("batch", data.batch.shape)

            loss = train_batch(model=model, data=data, optimizer=optimizer, loss_fn=loss_fn)
            print(loss)

            losses.append(loss)

            batch_index += 1

        print("Epoch: ", e+1)
        print("Batches: ", batch_index)

    return losses


def plot_loss(losses, show=True):
    fig = pyplot.gcf()
    fig.set_size_inches(8,6)
    ax = pyplot.axes()
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Loss")
    x_loss = list(range(len(losses)))
    pyplot.plot(x_loss, losses)

    if show:
        pyplot.show()

    pyplot.close()


def main():
    nodes_path = "/home/ryan/data/test_gfase/paolo_ul_guppy6_run14/tri_phase_test_output/nodes.csv"
    edges_path = "/home/ryan/data/test_gfase/paolo_ul_guppy6_run14/tri_phase_test_output/edges.csv"

    data, name_to_id, id_to_name = parse_csv_as_torch_data(nodes_path=nodes_path, edges_path=edges_path)

    node_size = 3
    edge_size = 2
    global_size = node_size

    learning_rate = 1e-5

    model = Graphnet(node_size=node_size, edge_size=edge_size, global_size=global_size)

    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

    loss_fn = torch.nn.BCEWithLogitsLoss()

    batch_size = 1
    epochs = 10

    data_loader_train = NeighborLoader(data=data, batch_size=batch_size, num_neighbors=[-1] * 4)
    data_loader_test = DataLoader(dataset=[data], batch_size=batch_size)
    losses = train(model, data_loader_train, optimizer, loss_fn, epochs)

    print("Final loss:", sum(losses[-10:])/10)
    plot_loss(losses)

    model.eval()

    s = torch.nn.Sigmoid()
    with open("maybe_bubbles.csv", 'w') as file:
        file.write("name,y_true,y_predict,y_predict_logit,color\n")

        for data in data_loader_test:
            y_predict = model.forward(data)

            for i in range(y_predict.shape[0]):
                p = s(y_predict[i]).item()

                color = "gray"
                if p > 0.5:
                    color = "Tomato"

                file.write("%s,%.1f,%.3f,%f,%s\n" % (id_to_name[i], data.y[i].item(), p, y_predict[i].item(), color))


if __name__ == "__main__":
    main()

