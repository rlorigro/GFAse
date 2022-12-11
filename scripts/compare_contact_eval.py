from collections import defaultdict
from matplotlib.lines import Line2D
from matplotlib import pyplot
import matplotlib.ticker
import argparse
import os.path
import pandas
import seaborn
import numpy
import json


def iterate_length_csv(path):
    with open(path,'r') as file:
        for l,line in enumerate(file):
            l,n = list(map(int, line.strip().split(',')))
            yield l,n


def iterate_signal_csv(path):
    defaultdict(list)
    with open(path,'r') as file:
        for l,line in enumerate(file):
            tokens = line.strip().split(',')

            if len(tokens[-1]) == 0:
                tokens = tokens[:-1]

            name = tokens[0]
            data = list(map(float,tokens[1:]))

            yield name, data


def plot_subread_length_distribution(directory_path, ax, color, coverage_weighted):
    path = os.path.join(directory_path,"subread_lengths.csv")

    lengths = defaultdict(int)
    x = list()
    y = list()

    for l,n in iterate_length_csv(path):
        lengths[l] = n

    for a,b in sorted(lengths.items(), key=lambda x: x[0]):
        x.append(a)

        if coverage_weighted:
            y.append(a*b)
        else:
            y.append(b)

    x.append(x[-1])
    y.append(0)

    ax.plot(x,y,color=color)

    ax.set_xlim([0,4000])

    if coverage_weighted:
        ax.set_ylabel("Coverage")
    else:
        ax.set_ylabel("Frequency (#)")

    ax.set_xlabel("Alignment length")


def plot_subread_count_distribution(directory_path, ax, color):
    path = os.path.join(directory_path,"subread_counts.csv")

    lengths = defaultdict(int)
    x = list()
    y = list()

    for l,n in iterate_length_csv(path):
        lengths[l] = n

    for a,b in sorted(lengths.items(), key=lambda x: x[0]):
        x.append(a)
        y.append(a*b)

    ax.plot(x,y,color=color)

    ax.set_xlim([0,60])
    ax.set_xlabel("# Alignments per read (mapq>0)")
    ax.set_ylabel("Frequency")


def plot_signal_distribution(directory_path, ax, color):
    path = os.path.join(directory_path,"phasing_summary.csv")

    x = list()
    y = list()

    for name,data in iterate_signal_csv(path):
        if name == "":
            x = data

        if name == "signal_ratio":
            y = list(map(numpy.log10, data))

    ax.plot(x[1:], y[1:], color)
    ax.set_xlabel("mapQ")
    ax.set_ylabel("Signal ratio (log10)")


def plot_contact_distribution(directory_path, ax, color):
    path = os.path.join(directory_path,"phasing_summary.csv")

    x = list()
    y = list()

    print(path)

    for name,data in iterate_signal_csv(path):
        if name == "":
            x = data

        if name == "n_consistent_contacts":
            y = numpy.array(data)
            print(y)

        if name == "n_inconsistent_contacts":
            y += numpy.array(data)
            print(data)
            print(y)

    ax.plot(x[3:], y[3:], color)
    ax.set_xlabel("mapQ")
    ax.set_ylabel("# contacts")


def main():
    paths = [
        "/home/ryan/data/test_gfase/terra/porec_qc/HG002_PAM_10290_contact_eval_results/",
        "/home/ryan/data/test_gfase/terra/porec_qc/HG002_PAM_10354_contact_eval_results/",
        "/home/ryan/data/test_gfase/terra/porec_qc/HG002_PAM_10393_contact_eval_results/",
        "/home/ryan/data/test_gfase/terra/hic_qc/HG002.HiC_1_S1_contact_eval_results",
        "/home/ryan/data/test_gfase/terra/hic_qc/HG002.HiC_1_S2_contact_eval_results",
        "/home/ryan/data/test_gfase/terra/hic_qc/HG002.HiC_1_S3_contact_eval_results",
    ]

    labels = [
        "PoreC_ONT_PAM_10290",
        "PoreC_ONT_PAM_10354",
        "PoreC_ONT_PAM_10393",
        "HiC_HPRC_1_S1",
        "HiC_HPRC_1_S2",
        "HiC_HPRC_1_S3",
    ]

    colors = [
        "#3496A5",
        "#3496A5",
        "#3496A5",
        "#D6418D",
        "#D6418D",
        "#D6418D",
        "#D6418D"
    ]

    f,axes = pyplot.subplots(2,2)

    for i,p in enumerate(paths):
        plot_subread_length_distribution(p, axes[0,0], colors[i], False)
        plot_subread_count_distribution(p, axes[0,1], colors[i])
        plot_signal_distribution(p, axes[1,0], colors[i])
        plot_contact_distribution(p, axes[1,1], colors[i])

    f.set_size_inches(12,8)

    custom_lines = list()
    for i in range(len(colors)):
        custom_lines.append(Line2D([0], [0], color=colors[i], lw=4))

    axes[0,1].legend(custom_lines, labels)

    f.tight_layout()
    pyplot.savefig("contact_comparison.png", dpi=200)
    pyplot.savefig("contact_comparison.pdf", dpi=200)

    pyplot.show()


if __name__ == "__main__":
    main()
