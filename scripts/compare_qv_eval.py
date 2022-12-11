from matplotlib.lines import Line2D
from matplotlib import pyplot
import matplotlib.ticker
import argparse
import os.path
import pandas
import seaborn
import numpy
import json


def main():

    paths = [
        "/home/ryan/data/test_gfase/terra/qv_eval/shasta_HG002_standard_m1_p30.txt",
        "/home/ryan/data/test_gfase/terra/qv_eval/shasta_HG002_UL_m1_p50.txt",
        "/home/ryan/data/test_gfase/terra/qv_eval/shasta_HG002_UUL_p50_m1.txt",
        "/home/ryan/data/test_gfase/terra/qv_eval/shasta_HG002_r10_fast_mode_m1_p20.txt",
        "/home/ryan/data/test_gfase/terra/qv_eval/shasta_hg002_r10_slow_mode_m1_p10.txt",
        "/home/ryan/data/test_gfase/terra/qv_eval/HG002_verkko_trio_production.txt",
        # "/home/ryan/data/test_gfase/terra/qv_eval/HG002_hifiasm_hic_v016.txt",
        "/home/ryan/data/test_gfase/terra/qv_eval/HG002_hifiasm_trio_v016.txt",
    ]

    labels = [
        "Shasta R9 standard",
        "Shasta R9 UL",
        "Shasta R9 UUL",
        "Shasta R10 400bps 1fc",
        "Shasta R10 260bps 2fc",
        "Verkko Trio 'production'",
        # "Hifiasm HiC",
        "Hifiasm Trio",
        # "Hifiasm HiC",
        # "Hifiasm Trio",
    ]

    colors = [
        "#3496A5",
        "#3496A5",
        "#3496A5",
        "#6CDB9F",
        "#6CDB9F",
        "#f6b26bff",
        # "#F27E7E",
        "#F27E7E",
    ]

    qvs = list()

    y_max = 0

    for path in paths:
        with open(path,'r') as file:
            for l,line in enumerate(file):
                if line.startswith("QV"):
                    data = line.strip().split('\t')
                    qv = float(data[2])
                    qvs.append(qv)

                    if qv > y_max:
                        y_max = qv

                    break

    width = 1  # the width of the bars

    fig, ax = pyplot.subplots()

    for i,value in enumerate(qvs):
        x = float(i)*width
        ax.text(x+width*0.1, 2, labels[i], rotation="vertical", color="white", ha="center")
        ax.bar(x,
               value,
               width,
               label=labels[i],
               color=colors[i],
               edgecolor="gray")

    custom_lines = list()
    for i in range(len(colors)):
        custom_lines.append(Line2D([0], [0], color=colors[i], lw=4))

    ax.axhline(30, linestyle="--", linewidth=0.5, zorder=-1)
    ax.axhline(40, linestyle="--", linewidth=0.5, zorder=-1)
    ax.axhline(50, linestyle="--", linewidth=0.5, zorder=-1)
    ax.axhline(60, linestyle="--", linewidth=0.5, zorder=-1)

    ax2 = ax.twinx()
    ax2.set_ylabel("Identity")
    y_ticks = [10,20,30,40,50]
    ax2.set_yticks(y_ticks)
    ax2.set_yticklabels([1-(10**(-p/10)) for p in y_ticks])

    ax.set_ylim([0,y_max*1.1])
    ax2.set_ylim([0,y_max*1.1])

    pyplot.xticks([])

    ax.set_title("GIAB HG002 QV (Yak Illumina)")
    ax.set_xlabel("Assembly")
    ax.set_ylabel("QV")

    fig.set_size_inches(1 + len(paths)*0.5,4)

    fig.tight_layout()
    pyplot.savefig("qv_comparison.png", dpi=200)
    pyplot.savefig("qv_comparison.pdf", dpi=200)

    pyplot.show()


if __name__ == "__main__":
    main()
