#!/usr/bin/env python3
from matplotlib import pyplot
from matplotlib.patches import Rectangle
import vcf

from modules.IterativeHistogram import IterativeHistogram
from collections import defaultdict
import argparse
import numpy
import math
import sys
import os
import re


def get_chromosome_ordering(name):
    if name.startswith("chr"):
        name = name.split("chr")[-1]

    ordinal = 0

    has_alpha = False
    has_numeric = False
    for c in name:
        if c.isalpha():
            has_alpha = True
        elif c.isnumeric():
            has_numeric = True

    if has_alpha and not has_numeric:
        weight = 1
        for c in name:
            ordinal += ord(c) * weight
            weight *= 0.1

    elif has_numeric and not has_alpha:
        ordinal = int(name)

    else:
        print("WARNING: region name is both numeric and alphabetical: " + name)
        return 9999999

    return ordinal


def parse_contig_header_line(line):
    data = re.split("[#<>,= ]", line.strip())

    name = None
    length = None

    for i,item in enumerate(data):
        if item == "ID":
            name = data[i+1]
        if item == "length":
            length = int(data[i+1])

    return name,length


"""
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample
"""
class VCFElement:
    def __init__(self, line):
        self.chrom,self.pos,self.id,self.ref,self.alt,self.qual,self.filter,self.info,self.format,self.sample = line.strip().split()

        self.pos = int(self.pos)
        self.qual = int(self.qual) if self.qual.isalpha() else 0

    def __str__(self):
        return '\t'.join(list(map(str,[self.chrom,self.pos,self.id,self.ref,self.alt,self.qual,self.filter,self.info,self.format,self.sample])))


"""
    # chrom  chromStart      chromEnd        name    gieStain
    # gneg      688 ************************************************************
    # gpos50    121 ***********
    # gpos75    89  ********
    # gpos25    87  ********
    # gpos100   81  *******
    # acen      48  ****
    # gvar      17  *
    # stalk     5
"""
class IdeogramElement:
    def __init__(self, line):
        data = line.strip().split()

        if len(data) == 5:
            chrom,chromStart,chromEnd,name,gieStain = data
        elif len(data) == 4:
            chrom,chromStart,chromEnd,gieStain = data
            name = None

        self.chrom = chrom
        self.start = int(chromStart)
        self.stop = int(chromEnd)
        self.name = name
        self.stain = gieStain
        self.color = None

        self.assign_color()

    def assign_color(self):
        if self.stain == "gneg":
            self.color = [0,0,0,0]
        elif self.stain == "gpos50":
            self.color = [0,0,0,0.5]
        elif self.stain == "gpos75":
            self.color = [0,0,0,0.75]
        elif self.stain == "gpos25":
            self.color = [0,0,0,0.25]
        elif self.stain == "gpos100":
            self.color = [0,0,0,1.0]
        elif self.stain == "acen":
            self.color = "darkred"
        elif self.stain == "gvar":
            self.color = "lightcoral"
        elif self.stain == "stalk":
            self.color = "lightcoral"
        else:
            print("WARNING: no color assignment for stain value: " + self.stain)
            return None

    def __str__(self):
        return '\t'.join(list(map(str,[self.chrom, self.start, self.stop, self.name, self.stain, self.color])))


def parse_vcf(file, lengths, variants_per_chromosome, bin_size, vcf_index):
    for l,line in enumerate(file):
        if "##contig" in line:
            name,length = parse_contig_header_line(line)

            if name not in lengths:
                lengths[name] = length
                print(name,length)

        elif line[0] != "#":
            e = VCFElement(line)

            if e.chrom not in variants_per_chromosome:
                length = lengths[e.chrom]
                n_bins = int(math.ceil(float(length)/bin_size))

                variants_per_chromosome[e.chrom][0] = IterativeHistogram(start=0, stop=length, n_bins=n_bins)
                variants_per_chromosome[e.chrom][1] = IterativeHistogram(start=0, stop=length, n_bins=n_bins)

            variants_per_chromosome[e.chrom][vcf_index].update(e.pos)


def plot_ideograms(axes, ideograms_per_chromosome, n_rows):

    for i,item in enumerate(sorted(ideograms_per_chromosome.items(), key=lambda x: get_chromosome_ordering(x[0]))):
        name,_ = item

        row_index = int(math.floor(float(i)/n_rows))
        column_index = 3*(i % n_rows)

        print(n_rows, axes.shape)
        print(name, row_index, column_index)

        chromosome_ends = defaultdict(int)

        if name in ideograms_per_chromosome:
            for item in ideograms_per_chromosome[name]:
                rect = Rectangle((item.start,-1),(item.stop-item.start),2,linewidth=0,facecolor=item.color, zorder=0)

                # Add the patch to the Axes
                axes[column_index+1][row_index].add_patch(rect)

                if item.stop > chromosome_ends[name]:
                    chromosome_ends[name] = item.stop

        else:
            print("WARNING: name not found in ideograms: " + name)

        rect = Rectangle((0,-1),(chromosome_ends[name]), 2, linewidth=0.5, fill=None, zorder=-1, color="C0")
        axes[column_index+1][row_index].add_patch(rect)

        axes[column_index+1][row_index].set_ylabel(name, rotation=0, labelpad=30)

        axes[column_index][row_index].spines['top'].set_visible(False)
        axes[column_index][row_index].spines['right'].set_visible(False)
        axes[column_index][row_index].spines['bottom'].set_visible(False)
        axes[column_index][row_index].spines['left'].set_visible(False)

        axes[column_index+1][row_index].spines['top'].set_visible(False)
        axes[column_index+1][row_index].spines['right'].set_visible(False)
        axes[column_index+1][row_index].spines['bottom'].set_visible(False)
        axes[column_index+1][row_index].spines['left'].set_visible(False)
        axes[column_index+1][row_index].set_xticks([])
        axes[column_index+1][row_index].set_yticks([])

        axes[column_index+2][row_index].spines['top'].set_visible(False)
        axes[column_index+2][row_index].spines['right'].set_visible(False)
        axes[column_index+2][row_index].spines['bottom'].set_visible(False)
        axes[column_index+2][row_index].spines['left'].set_visible(False)

        if column_index != axes.shape[0]:
            axes[column_index][row_index].set_xticks([])
            axes[column_index+2][row_index].set_xticks([])
        else:
            axes[column_index][-1].set_xlabel("Coordinate (Mbp)")


def plot_variants(figure, axes, n_rows, bin_size, variants_per_chromosome, lengths):

    for i,item in enumerate(sorted(variants_per_chromosome.items(), key=lambda x: get_chromosome_ordering(x[0]))):
        name,[fn_iterative_histogram,fp_iterative_histogram] = item

        row_index = int(math.floor(float(i)/n_rows))
        column_index = 3*(i % n_rows)

        print(n_rows, axes.shape)

        print(name, row_index, column_index)

        fp_frequencies = fp_iterative_histogram.get_histogram()
        fn_frequencies = -fn_iterative_histogram.get_histogram()

        assert(len(fp_frequencies) == len(fn_frequencies))

        x = fp_iterative_histogram.get_bin_centers()

        # x = [float(i)/2 for i in range(len(fp_frequencies))]

        axes[column_index][row_index].bar(x=x, height=fp_frequencies, width=bin_size, color="C1")
        axes[column_index+2][row_index].bar(x=x, height=fn_frequencies, width=bin_size, color="C1")

        # axes[column_index][row_index].set_ylim([0,20])
        axes[column_index+1][row_index].set_ylim([-1.2,1.2])
        # axes[column_index+2][row_index].set_ylim([-20,0])

        axes[column_index][row_index].set_xlim([-0.1*(max(lengths.values())),max(lengths.values())*1.1])
        axes[column_index+1][row_index].set_xlim([-0.1*(max(lengths.values())),max(lengths.values())*1.1])
        axes[column_index+2][row_index].set_xlim([-0.1*(max(lengths.values())),max(lengths.values())*1.1])


def main(fn_vcf_path, fp_vcf_path, ideogram_path):
    variants_per_chromosome = defaultdict(lambda: [None,None])
    lengths = dict()

    bin_size = 200_000

    ideograms_per_chromosome = defaultdict(list)

    with open(fp_vcf_path, 'r') as fp_file, open(fn_vcf_path, 'r') as fn_file:
        parse_vcf(fn_file, lengths, variants_per_chromosome, bin_size, 0)
        parse_vcf(fp_file, lengths, variants_per_chromosome, bin_size, 1)

    if ideogram_path is not None:
        with open(ideogram_path, 'r') as file:
            for l,line in enumerate(file):
                if l == 0:
                    continue

                e = IdeogramElement(line)

                name = e.chrom.split("chr")[-1]
                if name in variants_per_chromosome:
                    ideograms_per_chromosome[e.chrom].append(e)

                    print(e)
                else:
                    print("WARNING: ideogram name not found in VCF: " + name)

    n_rows = int(math.ceil(float(len(variants_per_chromosome))/2))
    figure, axes = pyplot.subplots(nrows=n_rows*3, ncols=2, sharey=False, sharex=True, gridspec_kw={'height_ratios': [2, 1, 2]*n_rows})

    max_ylim = 0
    for i,item in enumerate(sorted(variants_per_chromosome.items(), key=lambda x: get_chromosome_ordering(x[0]))):
        name,[fn_iterative_histogram,fp_iterative_histogram] = item
        m = max([max(fp_iterative_histogram.get_histogram()), max(fn_iterative_histogram.get_histogram())])

        if m > max_ylim:
            max_ylim = m

    plot_variants(figure, axes, n_rows, bin_size, variants_per_chromosome, lengths)
    plot_ideograms(axes, ideograms_per_chromosome, n_rows)

    figure.set_size_inches(24,32)
    pyplot.savefig("vcf_distribution.png", dpi=200)

    pyplot.show()
    pyplot.close()

    print(max_ylim)


def parse_comma_separated_string(s):
    return s.strip().split(',')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-fp",
        required=True,
        type=str,
        help="Input vcf file containing false positive records, can be a comma separated list"
    )

    parser.add_argument(
        "-fn",
        required=True,
        type=str,
        help="Input vcf file containing false negative records, can be a comma separated list"
    )

    parser.add_argument(
        "--ideogram",
        required=False,
        type=str,
        default=None,
        help="Ideogram 'cytoBandIdeo' table format from UCSC genome table browser"
    )

    args = parser.parse_args()

    main(fn_vcf_path=args.fn, fp_vcf_path=args.fp, ideogram_path=args.ideogram)
