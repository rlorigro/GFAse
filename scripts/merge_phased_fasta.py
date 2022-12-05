#!/usr/bin/env python3
import argparse
import sys
import os


def iterate_fasta_elements(path):
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


def iterate_fasta_names(path):
    with open(path, 'r') as file:
        for l,line in enumerate(file):
            if line.startswith('>'):
                name = line.strip()[1:].split(' ')[0]
                yield name


def write_sequence_to_fasta(name, sequence, file):
    file.write('>')
    file.write(name)
    file.write('\n')
    file.write(sequence)
    file.write('\n')


def main(path_a, path_b, output_directory):
    if not len(output_directory) == 0:
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

    output_csv_path = os.path.join(output_directory,"phases.csv")
    output_fasta_path = os.path.join(output_directory,"merged.fasta")

    with open(output_csv_path, 'w') as out_csv, open(output_fasta_path, 'w') as out_fasta:
        set_a = set()
        set_b = set()

        for name in iterate_fasta_names(path_a):
            set_a.add(name)

        for name in iterate_fasta_names(path_b):
            set_b.add(name)

        disjoint = set_a.isdisjoint(set_b)
        if not disjoint:
            sys.stderr.write("WARNING: non-disjoint contig names, adding suffixes '_hap_a' and '_hap_b' to fasta\n")

        for name,sequence in iterate_fasta_elements(path_a):
            if not disjoint:
                name += "_hap_a"

            write_sequence_to_fasta(name, sequence, out_fasta)
            out_csv.write("%s,%d\n" % (name, -1))

        for name,sequence in iterate_fasta_elements(path_b):
            if not disjoint:
                name += "_hap_b"

            write_sequence_to_fasta(name, sequence, out_fasta)
            out_csv.write("%s,%d\n" % (name, 1))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-a",
        required=True,
        type=str,
        help="Input fasta haplotype to be merged"
    )

    parser.add_argument(
        "-b",
        required=True,
        type=str,
        help="Input fasta haplotype to be merged"
    )

    parser.add_argument(
        "-o",
        required=True,
        type=str,
        help="Output directory"
    )

    args = parser.parse_args()

    main(path_a=args.a, path_b=args.b, output_directory=args.o)
