import argparse


def get_path_lengths(gfa_path):
    paths = dict()
    node_lengths = dict()

    with open(gfa_path, 'r') as file:
        for l,line in enumerate(file):

            if line.startswith('S'):
                tokens = line.strip().split('\t')
                node_name = tokens[1]

                if len(tokens) < 3:
                    print("WARNING: empty sequence on line " + str(l) + " in file " + gfa_path)
                    node_lengths[node_name] = 0
                else:
                    node_lengths[node_name] = len(tokens[2])

            elif line.startswith('P'):

                tokens = line.strip().split()

                path_name = tokens[1]
                nodes = [x[:-1] for x in tokens[2].split(',')]

                paths[path_name] = nodes

        for item in sorted(paths.items()):
            l = 0
            for node_name in item[1]:
                l += node_lengths[node_name]

            print(item[0], l)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        required=True,
        type=str
    )

    args = parser.parse_args()

    get_path_lengths(args.i)
