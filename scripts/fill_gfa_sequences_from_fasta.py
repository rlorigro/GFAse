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


def main(gfa_path, fasta_path, output_path):
    output_directory = os.path.dirname(output_path)

    if not len(output_directory) == 0:
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

    fasta_sequences = dict()
    for name,sequence in iterate_fasta(fasta_path):
        fasta_sequences[name] = sequence

    with open(output_path, 'w') as out_file, open(gfa_path, 'r') as gfa_file:
        for l,line in enumerate(gfa_file):

            # Substitute the Fasta sequence if this line contains sequence.
            # Otherwise don't modify the line before writing
            if line.startswith("S"):
                data = line.split() + ["\n"]
                name = data[1]
                data[2] = fasta_sequences[name]
                line = '\t'.join(data)

            # Write line
            out_file.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-g",
        required=True,
        type=str,
        help="Input GFA, a template GFA which will have its sequences replaced by any of the Fasta sequences with matching names"
    )

    parser.add_argument(
        "-f",
        required=True,
        type=str,
        help="Input Fasta, sequences from here will replace the output GFA sequences"
    )

    parser.add_argument(
        "-o",
        required=True,
        type=str,
        help="Output directory, where the output GFA is written"
    )

    args = parser.parse_args()

    main(gfa_path=args.g, fasta_path=args.f, output_path=args.o)
