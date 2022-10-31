#!/usr/bin/env python3
import argparse
import os


def iterate_fasta(path):
    name = ""
    sequence = ""

    with open(path, 'r') as file:
        for l,line in enumerate(file):
            if line.startswith('>'):
                # Output previous sequence
                if l > 0:
                    yield name, sequence

                name = line.strip()[1:].split(' ')[0]
                sequence = ""

            else:
                sequence += line.strip()

    yield name, sequence


def iterate_gfa(path):
    with open(path, 'r') as file:
        for l,line in enumerate(file):
            if line[0] == 'S' and line[1].isspace():
                data = line.strip().split()
                name = data[1]
                sequence = data[2]

                yield name, sequence


def write_sequence_to_fasta(name, sequence, file):
    file.write('>')
    file.write(name)
    file.write('\n')
    file.write(sequence)
    file.write('\n')


def main(path, output_directory):

    if not len(output_directory) == 0:
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

    output_path_0 = os.path.join(output_directory,"hap_0.fasta")
    output_path_1 = os.path.join(output_directory,"hap_1.fasta")

    iterator = iterate_fasta(path)
    if path.endswith(".gfa"):
        iterator = iterate_gfa(path)

    with open(output_path_0, 'w') as hap_0, open(output_path_1, 'w') as hap_1:
        for name,sequence in iterator:
            is_phased = (not name.startswith("UR") and name.count('.') > 0)

            print(name, is_phased)

            if is_phased:
                if name.endswith(".0"):
                    write_sequence_to_fasta(name,sequence,hap_0)
                elif name.endswith(".1"):
                    write_sequence_to_fasta(name,sequence,hap_1)
                else:
                    exit("Error: sequence name with phased format has invalid suffix: " + name)
            else:
                write_sequence_to_fasta(name,sequence,hap_0)
                write_sequence_to_fasta(name,sequence,hap_1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        required=True,
        type=str,
        help="Input fasta to be split"
    )

    parser.add_argument(
        "-o",
        required=True,
        type=str,
        help="Output directory"
    )

    args = parser.parse_args()

    main(path=args.i, output_directory=args.o)
