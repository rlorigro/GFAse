from matplotlib import pyplot
import numpy
import argparse
import networkx
from collections import defaultdict


class IdMap:
    def __init__(self):
        self.name_to_id = dict()
        self.id_to_name = list()

    def try_add(self, name):
        if name not in self.name_to_id:
            id = len(self.name_to_id)

            self.name_to_id[name] = id
            self.id_to_name.append(name)

            return id
        else:
            return self.name_to_id.get(name)

    def get_id(self, name):
        return self.name_to_id.get(name)

    def get_name(self, id):
        return self.id_to_name[id]


def iterate_contacts(csv_path):
    with open(csv_path, 'r') as file:
        for l,line in enumerate(file):
            if l == 0:
                continue

            data = line.strip().split(',')
            data[-1] = int(data[-1])

            yield data


def main(csv_path):
    g = networkx.Graph()
    id_map = IdMap()

    edges = defaultdict(lambda: defaultdict(int))

    for name_0, name_1, n_matches in iterate_contacts(csv_path):
        id_0 = id_map.try_add(name_0)
        id_1 = id_map.try_add(name_1)

        g.add_node(id_0, color="C0")
        g.add_node(id_1, color="C0")

        edges[id_0][id_1] = n_matches

    for id_0 in edges.keys():
        for id_1,count in edges[id_0].items():
            name_0 = id_map.get_name(id_0)
            name_1 = id_map.get_name(id_1)

            g.add_edge(id_0, id_1, color="C1", weight=max(1,numpy.log10(count/10)))

    pos = networkx.spring_layout(g, iterations=100)
    node_colors = networkx.get_node_attributes(g,'color').values()
    edge_colors = networkx.get_edge_attributes(g,'color').values()
    weights = networkx.get_edge_attributes(g,'weight').values()

    networkx.draw(
        g,
        pos,
        edge_color=edge_colors,
        width=list(weights),
        node_size=10,
        # with_labels=True,
        node_color=node_colors)

    pyplot.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        required=True,
        type=str,
        help="Input file of alignments.csv"
    )

    args = parser.parse_args()

    main(csv_path=args.i)
