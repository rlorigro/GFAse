from modules.IterativeHistogram import *
from collections import defaultdict
from matplotlib.lines import Line2D
from matplotlib import pyplot
import os.path


def iterate_length_csv(path):
    with open(path,'r') as file:
        for l,line in enumerate(file):
            if l == 0:
                continue

            l,n = list(map(int, line.strip().split('\t')))
            yield l,n


def plot_subread_length_distribution(directory_path, ax, color, coverage_weighted, normalized):
    path = os.path.join(directory_path,"unabridged_distribution.tsv")
    x_max = 400_000

    lengths = defaultdict(int)
    histogram = IterativeHistogram(0,1_200_000,600)

    for l,n in iterate_length_csv(path):
        lengths[l] = n

    for x,y in sorted(lengths.items(), key=lambda x: x[0]):
        histogram.update(x,y)

    x = histogram.get_bin_centers()
    y = histogram.get_histogram()

    y_label = ""

    if normalized:
        y_label += "Normalized "

    if coverage_weighted:
        y *= x
        y_label += "Coverage"
    else:
        y_label += "Frequency (#)"

    if normalized:
        y /= numpy.sum(y)

    ax.plot(x,y,color=color,linewidth=0.7)

    ax.set_xlabel("Read length (bp)")
    ax.set_ylabel(y_label)

    ax.set_xlim([0,x_max])


def plot_cumulative_subread_length_distribution(directory_path, ax, color, coverage_weighted):
    path = os.path.join(directory_path,"unabridged_distribution.tsv")
    x_max = 400_000
    coverage_constant = 3_150_000_000

    lengths = defaultdict(int)
    histogram = IterativeHistogram(0,1_200_000,500)

    for l,n in iterate_length_csv(path):
        lengths[l] = n

    for x,y in sorted(lengths.items(), key=lambda x: x[0]):
        histogram.update(x,y)

    x = histogram.get_bin_centers()
    y = histogram.get_histogram()

    y_label = ""
    if coverage_weighted:
        y *= x
        y_label += "Cumulative coverage (1x = %dbp)" % coverage_constant
    else:
        y_label += "Cumulative frequency (#)"

    total_coverage = numpy.sum(y,dtype=numpy.float64)/coverage_constant
    y = total_coverage - numpy.cumsum(y, dtype=numpy.float64)/coverage_constant
    ax.plot(x,y,color=color,linewidth=0.7)

    ax.set_xlabel("Read length (bp)")
    ax.set_ylabel(y_label)

    ax.set_xlim([0,x_max])

    ax.grid(visible=True, which='major', axis='both')


def main():
    paths = [
        "/home/ryan/data/test_gfase/HG002_UL_R10/HG002_UL-20221211T015417Z-001/HG002_UL/ULCIR_HG002_1",
        "/home/ryan/data/test_gfase/HG002_UL_R10/HG002_UL-20221211T015417Z-001/HG002_UL/ULCIR_HG002_2",
        "/home/ryan/data/test_gfase/HG002_UL_R10/HG002_UL-20221211T015417Z-001/HG002_UL/UL_HG002_1",
        "/home/ryan/data/test_gfase/HG002_UL_R10/HG002_UL-20221211T015417Z-001/HG002_UL/UL_HG002_2",
        "/home/ryan/data/test_gfase/HG002_UL_R10/HG002_UL-20221211T015417Z-001/HG002_UL/UL_HG002_3",
        "/home/ryan/data/test_gfase/HG002_UL_R10/HG002_UL-20221211T015417Z-001/HG002_UL/UL_HG002_4",
        "/home/ryan/data/test_gfase/HG002_UL_R10/HG002_UL-20221211T015417Z-001/HG002_UL/ULNEB_HG002_1",
        "/home/ryan/data/test_gfase/HG002_UL_R10/HG002_UL-20221211T015417Z-001/HG002_UL/ULNEB_HG002_2",
    ]

    labels = [
        "ULCIR_HG002_1",
        "ULCIR_HG002_2",
        "UL_HG002_1",
        "UL_HG002_2",
        "UL_HG002_3",
        "UL_HG002_4",
        "ULNEB_HG002_1",
        "ULNEB_HG002_2",
    ]

    colors = [
        "#3496A5",
        "#3496A5",
        "#f6b26bff",
        "#f6b26bff",
        "#f6b26bff",
        "#f6b26bff",
        "#F27E7E",
        "#F27E7E",
    ]

    f,axes = pyplot.subplots(ncols=2)

    for i,p in enumerate(paths):
        plot_subread_length_distribution(p, axes[0], colors[i], True, True)
        plot_cumulative_subread_length_distribution(p, axes[1], colors[i], True)

    f.set_size_inches(12,5)

    custom_lines = list()
    for i in range(len(colors)):
        custom_lines.append(Line2D([0], [0], color=colors[i], lw=4))

    axes[0].legend(custom_lines, labels)

    f.tight_layout()
    pyplot.savefig("length_comparison.png", dpi=200)
    pyplot.savefig("length_comparison.pdf", dpi=200)

    pyplot.show()


if __name__ == "__main__":
    main()
