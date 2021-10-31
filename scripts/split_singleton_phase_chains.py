import argparse


def split(fasta_path, csv_path):
    singletons = set()

    """
    Parse a file with the following format to find all the names that correspond to singleton paths (length=1)

    path_name,n_steps,nodes
    4.p,1,tip.0+
    1.m,2,a+ 0.1+
    4.m,1,tip.1+
    3.p,3,1.1+ m+ 4.0+
    3.m,3,1.0+ m+ 4.1+

    """
    with open(csv_path, 'r') as file:
        for l,line in enumerate(file):
            if l == 0:
                continue

            data = line.strip().split(',')

            length = int(data[1])

            if length == 1:
                singletons.add(data[0])

    singletons_output_path = '.'.join(fasta_path.split('.')[:-1]) + "_singletons.fasta"
    chains_output_path = '.'.join(fasta_path.split('.')[:-1]) + "_chains.fasta"

    # Iterate the FASTA and split depending on whether the sequence was a singleton chain or not
    with open(fasta_path, 'r') as file, open(singletons_output_path, 'w') as singletons_out_file, open(chains_output_path, 'w') as chains_out_file:
        is_singleton = False
        name = None

        for l,line in enumerate(file):
            if line[0] == '>':
                name = line.strip()[1:]
                is_singleton = (name in singletons)

            else:
                if is_singleton:
                    singletons_out_file.write(">")
                    singletons_out_file.write(name)
                    singletons_out_file.write('\n')
                    singletons_out_file.write(line)
                else:
                    chains_out_file.write(">")
                    chains_out_file.write(name)
                    chains_out_file.write('\n')
                    chains_out_file.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fasta","-f",
        required=True,
        type=str
    )

    parser.add_argument(
        "--csv","-c",
        required=True,
        type=str
    )

    args = parser.parse_args()

    split(args.fasta, args.csv)
