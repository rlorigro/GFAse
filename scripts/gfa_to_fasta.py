import argparse
import sys
import os


def main(path, output_path):
    output_directory = os.path.dirname(output_path)

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    if not output_path.endswith(".fasta") or output_path.endswith(".fa"):
        exit("ERROR: output path does not have FASTA suffix: " + output_path)

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

    parser.add_argument(
        "-o",
        required=True,
        type=str,
        help="Output directory where [filename].fasta will be written"
    )

    args = parser.parse_args()

    main(path=args.i, output_path=args.o)
