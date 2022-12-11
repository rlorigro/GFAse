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
        "/home/ryan/data/test_gfase/terra/sv_eval/shasta_HG002_standard_m1_p30.json",
        "/home/ryan/data/test_gfase/terra/sv_eval/shasta_HG002_UL_m1_p50.json",
        "/home/ryan/data/test_gfase/terra/sv_eval/shasta_HG002_UUL_p50_m1.json",
        "/home/ryan/data/test_gfase/terra/sv_eval/shasta_HG002_r10_fast_mode_m1_p20.json",
        "/home/ryan/data/test_gfase/terra/sv_eval/shasta_hg002_r10_slow_mode_m1_p10.json",
        "/home/ryan/data/test_gfase/terra/sv_eval/HG002_verkko_trio_full_coverage.json",
        "/home/ryan/data/test_gfase/terra/sv_eval/HG002_verkko_trio_production.json",
        "/home/ryan/data/test_gfase/terra/sv_eval/HG002_hifiasm_trio_v016.json",
    ]

    labels = [
        "Shasta R9 standard",
        "Shasta R9 UL",
        "Shasta R9 UUL",
        "Shasta R10 400bps 1FC ",
        "Shasta R10 260bps 2FC",
        "Verkko Trio 'full coverage'",
        "Verkko Trio 'production'",
        "Hifiasm Trio",
    ]

    colors = [
        "#3496A5",
        "#3496A5",
        "#3496A5",
        "#6CDB9F",
        "#6CDB9F",
        "#f6b26bff",
        "#f6b26bff",
        "#F27E7E",
        "#F27E7E",
    ]

    precisions = list()
    recalls = list()
    f1s = list()

    for path in paths:
        with open(path,'r') as file:
            data = json.load(file)

            print(os.path.basename(path))
            print("precision: ", data["precision"])
            print("recall: ", data["recall"])
            print("f1: ", data["f1"])

            precisions.append(data["precision"])
            recalls.append(data["recall"])
            f1s.append(data["f1"])

            # TODO: why don't these sum to the same value for all results?
            # print(data["TP-base"] + data["FP"] + data["FN"])

    width = 0.8/float(len(paths))  # the width of the bars

    fig, ax = pyplot.subplots()

    offset = -(float(len(paths) - 1)/2.0)*width
    for i,value in enumerate(precisions):
        x = offset + float(i)*width
        ax.text(x+width*0.1, 0.05, labels[i], rotation="vertical", color="white", ha="center")
        ax.bar(x,
               value,
               width,
               label=labels[i],
               color=colors[i],
               edgecolor="gray")

    offset += 1
    for i,value in enumerate(recalls):
        x = offset + float(i)*width
        ax.text(x+width*0.1, 0.05, labels[i], rotation="vertical", color="white", ha="center")
        ax.bar(x,
               value,
               width,
               label=labels[i],
               color=colors[i],
               edgecolor="gray")

    offset += 1
    for i,value in enumerate(f1s):
        x = offset + float(i)*width
        ax.text(x+width*0.1, 0.05, labels[i], rotation="vertical", color="white", ha="center")
        ax.bar(x,
               value,
               width,
               label=labels[i],
               color=colors[i],
               edgecolor="gray")

    custom_lines = list()
    for i in range(len(colors)):
        custom_lines.append(Line2D([0], [0], color=colors[i], lw=4))

    x_tick_labels = ["Precision", "Recall", "F1"]
    x_ticks = [0,1,2]

    minor_ticks = matplotlib.ticker.MultipleLocator(base=0.05)
    major_ticks = matplotlib.ticker.MultipleLocator(base=0.1)
    ax.yaxis.set_minor_locator(minor_ticks)
    ax.yaxis.set_major_locator(major_ticks)

    ax.set_xticklabels(x_tick_labels)
    ax.set_xticks(x_ticks)

    ax.axhline(1, linestyle="--", linewidth=0.5, zorder=-1)
    ax.axhline(0.95, linestyle="--", linewidth=0.5, zorder=-1)
    ax.axhline(0.9, linestyle="--", linewidth=0.5, zorder=-1)

    ax.set_title("GIAB HG002 Structural Variant Accuracy (Truvari)")

    fig.set_size_inches(3 + len(paths)*0.5,4)

    fig.tight_layout()
    pyplot.savefig("sv_comparison.png", dpi=200)
    pyplot.savefig("sv_comparison.pdf", dpi=200)

    pyplot.show()


if __name__ == "__main__":
    main()
