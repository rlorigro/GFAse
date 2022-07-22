import argparse
import sys


def main(path):
    output_path = '.'.join(path.split('.')[:-1]) + ".fasta"

    print(output_path)

    with open(path, 'r') as file, open(output_path,'w') as output_file:
        for l,line in enumerate(file):
            if line[0] == 'S' and line[1].isspace():
                data = line.strip().split()
                name = data[1]
                sequence = data[2]

                output_file.write('>')
                output_file.write(name)
                output_file.write('\n')
                output_file.write(sequence)
                output_file.write('\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        required=True,
        type=str,
        help="Input GFA to be phased"
    )

    args = parser.parse_args()

    main(path=args.i)
