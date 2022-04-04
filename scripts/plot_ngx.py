#!/usr/bin/env python3
from matplotlib import pyplot
import argparse
import numpy
import sys
import os


def get_name_of_unit(n):
    if n == 1_000_000_000:
        return "Gbp"
    if n == 1_000_000:
        return "Mbp"
    if n == 1_000:
        return "Kbp"
    if n == 0:
        return "bp"


def plot_ngx(lengths, genome_size, output_dir):
    fig = pyplot.figure()
    axes = pyplot.axes()

    fig.set_size_inches(12,10)

    total = sum(lengths)

    if total > genome_size:
        sys.stderr.write("WARNING: observed total sequence length is greater than genome size")

    unit = 10**max(0,round((numpy.log10(genome_size) - 3) / 3)*3)

    ng50 = None

    x_prev = 0
    s_prev = None
    for s in sorted(lengths, reverse=True):
        s /= unit

        x0 = x_prev
        x1 = x_prev + s

        if x0*unit >= genome_size/2 and x1*unit > genome_size/2:
            ng50 = s*unit

        pyplot.plot([x0,x1],[s,s], color="C1")

        if s_prev is not None:
            pyplot.plot([x0,x0],[s_prev,s], color="C1")

        x_prev = x1
        s_prev = s

    axes.axvline(genome_size/2/unit, linestyle='--', linewidth=0.6, color='gray')
    axes.ticklabel_format(style='plain')

    pyplot.xticks(rotation=45)
    pyplot.yticks(rotation=45)

    axes.set_xlabel("Cumulative coverage (%s)" % get_name_of_unit(unit))
    axes.set_ylabel("Length (%s)" % get_name_of_unit(unit))
    axes.set_title("NGx")

    all_lengths_txt_path = os.path.join(output_dir, "all_lengths.txt")
    ng50_txt_path = os.path.join(output_dir, "ng50.txt")
    fig_path = os.path.join(output_dir, "ngx.png")

    pyplot.savefig(fig_path, dpi=200)

    pyplot.close()

    with open(all_lengths_txt_path,'w') as file:
        for l in lengths:
            file.write(str(l))
            file.write('\n')

    with open(ng50_txt_path,'w') as file:
        file.write(str(ng50))
        file.write('\n')


def main(input_path, output_dir, genome_size):
    if os.path.exists(output_dir):
        exit("ERROR: output directory already exists")
    else:
        os.makedirs(output_dir)

    lengths = list()

    with open(input_path, 'r') as file:
        for l,line in enumerate(file):
            if line[0] == '>':
                lengths.append(0)
            else:
                # Increment most recent length by size of line (without newline char)
                lengths[-1] += len(line) - 1

    lengths = sorted(lengths, reverse=True)

    plot_ngx(lengths=lengths, genome_size=genome_size, output_dir=output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        required=True,
        type=str,
        help="Input fasta file"
    )

    parser.add_argument(
        "-o",
        required=True,
        type=str,
        help="Output directory"
    )

    parser.add_argument(
        "-g",
        required=True,
        type=int,
        help="Size of genome"
    )

    args = parser.parse_args()
    main(input_path=args.i, output_dir=args.o, genome_size=args.g)
