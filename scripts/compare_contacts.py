from modules.IterativeHistogram import *
from matplotlib import pyplot
from matplotlib import cm
import argparse
import math

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
            data = line.strip().split(',')
            data[-1] = int(data[-1])

            yield data


def iterate_bubbles_csv(csv_path):
    with open(csv_path, 'r') as file:
        for l,line in enumerate(file):
            if l == 0:
                continue

            data = line.strip().split(',')

            name = data[0]
            bubble_id = int(data[1])
            phase = bool(int(data[2]))

            yield name, bubble_id, phase


def main(path_a, path_b, bubbles_csv_path):
    allowed_names = {
        # "PR.1.1.1.0",
        # "PR.1.1.1.1",
        # "PR.1.3.54.0",
        # "PR.1.3.54.1",
        # "PR.1.5.5.0",
        # "PR.1.5.5.1",
        # "PR.60.3.159.0",
        # "PR.60.3.159.1",
        # "PR.60.5.131.0",
        # "PR.60.5.131.1",
        # "PR.60.7.1.0",
        # "PR.60.7.1.1"
    }

    g = networkx.Graph()

    id_map = IdMap()

    edges = defaultdict(lambda: defaultdict(lambda: [0,0]))
    intra_bubble_edges = defaultdict(lambda: set())

    phases = [set(), set()]

    total_contacts_a = 0
    total_contacts_b = 0
    max_contacts_a = 0
    max_contacts_b = 0

    colormap = cm.get_cmap('coolwarm', 256)

    prev_bubble_id = -1
    bubble_names = [None,None]
    for name,bubble_id,phase in iterate_bubbles_csv(bubbles_csv_path):
        print(name,bubble_id,phase)

        if bubble_id != prev_bubble_id:
            if len(allowed_names) == 0 or ((bubble_names[0] in allowed_names) and (bubble_names[1] in allowed_names)):
                intra_bubble_edges[bubble_names[0]].add(bubble_names[1])
                intra_bubble_edges[bubble_names[1]].add(bubble_names[0])

                phases[0].add(bubble_names[0])
                phases[1].add(bubble_names[1])

                print("adding bubble",bubble_names)

                bubble_names = [None,None]

        print(bubble_names)
        bubble_names[phase] = name
        print(bubble_names)
        prev_bubble_id = bubble_id

    intra_bubble_edges[bubble_names[0]].add(bubble_names[1])
    intra_bubble_edges[bubble_names[1]].add(bubble_names[0])

    for name_0,name_1,count in iterate_contacts(path_a):
        if name_0 == name_1:
            continue

        if len(allowed_names) > 0 and not ((name_0 in allowed_names) and (name_1 in allowed_names)):
            continue

        color_0 = 'grey' if name_0 in phases[0] else 'lightgrey'
        color_1 = 'grey' if name_1 in phases[0] else 'lightgrey'

        id_0 = id_map.try_add(name_0)
        id_1 = id_map.try_add(name_1)

        g.add_node(id_0, color=color_0)
        g.add_node(id_1, color=color_1)

        edges[id_0][id_1][0] = count
        edges[id_1][id_0][0] = count

        total_contacts_a += count

        if count > max_contacts_a:
            max_contacts_a = count

    for name_0,name_1,count in iterate_contacts(path_b):
        if name_0 == name_1:
            continue

        if len(allowed_names) > 0 and not ((name_0 in allowed_names) and (name_1 in allowed_names)):
            continue

        color_0 = 'darkgrey' if name_0 in phases[0] else 'lightgrey'
        color_1 = 'darkgrey' if name_1 in phases[0] else 'lightgrey'

        id_0 = id_map.try_add(name_0)
        id_1 = id_map.try_add(name_1)

        g.add_node(id_0, color=color_0)
        g.add_node(id_1, color=color_1)

        edges[id_0][id_1][1] = count
        edges[id_1][id_0][1] = count

        total_contacts_b += count

        if count > max_contacts_b:
            max_contacts_b = count

    max_normalized_contacts_a = max_contacts_a/total_contacts_a
    max_normalized_contacts_b = max_contacts_b/total_contacts_b

    max_normalized_contacts = max(max_normalized_contacts_a, max_normalized_contacts_b)

    n_bins = 100
    start = -5
    stop = 5
    interval = (stop-start)/n_bins

    signal_ratio_histogram_a = IterativeHistogram(start=start,stop=stop,n_bins=n_bins)
    signal_ratio_histogram_b = IterativeHistogram(start=start,stop=stop,n_bins=n_bins)

    signal_ratios_per_node_a = defaultdict(lambda: [0,0])
    signal_ratios_per_node_b = defaultdict(lambda: [0,0])

    e = 1e-12

    for id_0 in edges.keys():
        for id_1,[count_a,count_b] in edges[id_0].items():
            name_0 = id_map.get_name(id_0)
            name_1 = id_map.get_name(id_1)

            print(name_0, name_1)

            # Skip non-bubbles
            if (name_0 not in intra_bubble_edges) or (name_1 not in intra_bubble_edges):
                print("non-bubble")
                continue

            # Skip edges within a bubble
            if name_1 in intra_bubble_edges[name_0]:
                print("intra_bubble_edge", name_0, name_1)
                continue

            weight = float(math.log2(((float(count_a)+e)/float(total_contacts_a)) / ((float(count_b)+e)/float(total_contacts_b))))
            weight_scaled = weight/5
            weight_scaled = min(max(-1.0, weight_scaled), 1.0)

            color_weight = (weight_scaled + 1)/2.0
            color = colormap(color_weight)

            print(name_0, name_1, count_a, count_b, weight, weight_scaled, color_weight)

            is_cross_phase_edge = (name_0 in phases[0]) != (name_1 in phases[0])

            if is_cross_phase_edge:
                signal_ratios_per_node_a[name_0][0] += count_a
                signal_ratios_per_node_a[name_1][0] += count_a
                signal_ratios_per_node_b[name_0][0] += count_b
                signal_ratios_per_node_b[name_1][0] += count_b
            else:
                signal_ratios_per_node_a[name_0][1] += count_a
                signal_ratios_per_node_a[name_1][1] += count_a
                signal_ratios_per_node_b[name_0][1] += count_b
                signal_ratios_per_node_b[name_1][1] += count_b

            g.add_edge(id_0, id_1, color=color)
            # g.add_edge(id_0, id_1, weight=count_a+e, color=color)

    for v in signal_ratios_per_node_a.values():
        signal_ratio_histogram_a.update(math.log2((float(v[1])+e)/(float(v[0])+e)))

    for v in signal_ratios_per_node_b.values():
        signal_ratio_histogram_b.update(math.log2((float(v[1])+e)/(float(v[0])+e)))

    print(phases[0])
    print(phases[1])

    remapping = dict()

    for i in range(len(id_map.id_to_name)):
        remapping[i] = id_map.id_to_name[i]

    networkx.relabel_nodes(g, remapping, False)

    pos = networkx.bipartite_layout(g, nodes=phases[0])
    # pos = networkx.spring_layout(g)
    node_colors = networkx.get_node_attributes(g,'color').values()
    edge_colors = networkx.get_edge_attributes(g,'color').values()
    weights = networkx.get_edge_attributes(g,'weight').values()

    networkx.draw(
        g,
        pos,
        edge_color=edge_colors,
        with_labels=True,
        node_color=node_colors)

    fig = pyplot.figure()
    axes = pyplot.axes()

    x = signal_ratio_histogram_a.get_bin_centers()
    y = signal_ratio_histogram_a.get_histogram()
    print(y)
    pyplot.bar(x=x,height=y,color="C1",alpha=0.6,width=interval)

    x = signal_ratio_histogram_b.get_bin_centers()
    y = signal_ratio_histogram_b.get_histogram()
    print(y)
    pyplot.bar(x=x,height=y,color="C0",alpha=0.6,width=interval)

    axes.set_title("Signal ratio per node")
    axes.set_ylabel("Frequency (# of nodes)")
    axes.set_xlabel("Signal ratio: log2(consistent/inconsistent)")

    fig2 = pyplot.figure()
    axes2 = pyplot.axes()

    x = list()
    y = list()

    for node_name in signal_ratios_per_node_a.keys():
        n_inconsistent_a,n_consistent_a = signal_ratios_per_node_a[node_name]
        n_inconsistent_b,n_consistent_b = signal_ratios_per_node_b[node_name]

        ratio_a = math.log2((float(n_consistent_a)+e) / (float(n_inconsistent_a)+e))
        ratio_b = math.log2((float(n_consistent_b)+e) / (float(n_inconsistent_b)+e))

        print(n_inconsistent_a,n_consistent_a,ratio_a)
        print(n_inconsistent_b,n_consistent_b,ratio_b)

        x.append(ratio_a)
        y.append(ratio_b)

    axes2.scatter(x,y)

    pyplot.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-a",
        required=True,
        type=str,
        help="Input file of contacts.csv"
    )

    parser.add_argument(
        "-b",
        required=True,
        type=str,
        help="Input file of contacts.csv"
    )

    parser.add_argument(
        "--bubbles",
        required=True,
        type=str,
        help="Input file of bubbles.csv"
    )

    args = parser.parse_args()

    main(path_a=args.a, path_b=args.b, bubbles_csv_path=args.bubbles)
