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


def iterate_asmgene_results(paths):
    for path in paths:
        with open(path,'r') as file:
            ref_counts = defaultdict(int)
            sample_counts = defaultdict(int)

            for l,line in enumerate(file):
                if l == 0:
                    continue

                data = line.strip().split()

                key = data[1]
                ref_count = int(data[2])
                sample_count = int(data[3])

                print(key, ref_count, sample_count)

                ref_counts[key] = ref_count
                sample_counts[key] = sample_count

            yield ref_counts, sample_counts


def main():

    paths_diploid_only = [
        # "/home/ryan/data/test_gfase/terra/asmgene_eval_diploid_only/shasta_HG01993_p50.txt",
        # "/home/ryan/data/test_gfase/terra/asmgene_eval_diploid_only/shasta_HG02132_p50.txt",
        # "/home/ryan/data/test_gfase/terra/asmgene_eval_diploid_only/shasta_HG02647_p50.txt",
        # "/home/ryan/data/test_gfase/terra/asmgene_eval_diploid_only/shasta_HG03669_p50.txt",
        "/home/ryan/data/test_gfase/terra/asmgene_eval_diploid_only/shasta_HG002_standard_m1_p30.txt",
        "/home/ryan/data/test_gfase/terra/asmgene_eval_diploid_only/shasta_HG002_UL_m1_p50.txt",
        "/home/ryan/data/test_gfase/terra/asmgene_eval_diploid_only/shasta_HG002_UUL_p50_m1.txt",
        "/home/ryan/data/test_gfase/terra/asmgene_eval_diploid_only/shasta_HG002_r10_fast_mode_m1_p20.txt",
        "/home/ryan/data/test_gfase/terra/asmgene_eval_diploid_only/shasta_hg002_r10_slow_mode_m1_p10.txt",
        "/home/ryan/data/test_gfase/terra/asmgene_eval_diploid_only/verkko_trio_HG002_production.txt",
        "/home/ryan/data/test_gfase/terra/asmgene_eval_diploid_only/hifiasm_trio_HG002_v016.txt",
    ]

    paths_with_unphased = [
        # "/home/ryan/data/test_gfase/terra/asmgene_eval_w_unphased/shasta_HG01993_p50.txt",
        # "/home/ryan/data/test_gfase/terra/asmgene_eval_w_unphased/shasta_HG02132_p50.txt",
        # "/home/ryan/data/test_gfase/terra/asmgene_eval_w_unphased/shasta_HG02647_p50.txt",
        # "/home/ryan/data/test_gfase/terra/asmgene_eval_w_unphased/shasta_HG03669_p50.txt",
        "/home/ryan/data/test_gfase/terra/asmgene_eval_w_unphased/shasta_HG002_standard_m1_p30.txt",
        "/home/ryan/data/test_gfase/terra/asmgene_eval_w_unphased/shasta_HG002_UL_m1_p50.txt",
        "/home/ryan/data/test_gfase/terra/asmgene_eval_w_unphased/shasta_HG002_UUL_p50_m1.txt",
        "/home/ryan/data/test_gfase/terra/asmgene_eval_w_unphased/shasta_HG002_r10_fast_mode_m1_p20.txt",
        "/home/ryan/data/test_gfase/terra/asmgene_eval_w_unphased/shasta_hg002_r10_slow_mode_m1_p10.txt",
        "/home/ryan/data/test_gfase/terra/asmgene_eval_diploid_only/verkko_trio_HG002_production.txt",
        "/home/ryan/data/test_gfase/terra/asmgene_eval_diploid_only/hifiasm_trio_HG002_v016.txt",
    ]

    labels = [
        # "Shasta HG01993 R9 UL",
        # "Shasta HG02132 R9 UL",
        # "Shasta HG02647 R9 UL",
        # "Shasta HG03669 R9 UL",
        "Shasta R9 standard",
        "Shasta R9 UL",
        "Shasta R9 UUL",
        "Shasta R10 fast 1FC",
        "Shasta R10 slow 2FC",
        "Verkko Trio 'production'",
        "Hifiasm Trio",
    ]

    colors = [
        # "#3496A5",
        # "#3496A5",
        # "#3496A5",
        # "#3496A5",
        "#3496A5",
        "#3496A5",
        "#3496A5",
        "#6CDB9F",
        "#6CDB9F",
        "#f6b26bff",
        "#F27E7E",
    ]

    """
    H	Metric	genesToRef	phase_0
    X	full_sgl	34173	33884
    X	full_dup	0	60
    X	frag	0	16
    X	part50+	0	85
    X	part10+	0	7
    X	part10-	0	121
    X	dup_cnt	2034	1453
    X	dup_sum	6763	4591
    """

    full_single_proportions = [list(),list()]
    frag_proportions = [list(),list()]
    dupe_proportions = [list(),list()]
    multicopy_proportions = [list(),list()]

    for ref_counts,sample_counts in iterate_asmgene_results(paths_diploid_only):
        total_single_ref = ref_counts["full_sgl"] + ref_counts["full_dup"] + ref_counts["frag"]

        full_single_proportion = float(sample_counts["full_sgl"]) / float(total_single_ref)
        full_single_proportions[0].append(full_single_proportion)

        dupe_proportion = float(sample_counts["full_dup"]) / float(total_single_ref)
        dupe_proportions[0].append(dupe_proportion)

        frag_proportion = float(sample_counts["frag"]) / float(total_single_ref)
        frag_proportions[0].append(frag_proportion)

        multicopy_proportion = float(sample_counts["dup_cnt"]) / float(ref_counts["dup_cnt"])
        multicopy_proportions[0].append(multicopy_proportion)

    for ref_counts,sample_counts in iterate_asmgene_results(paths_with_unphased):
        total_single_ref = ref_counts["full_sgl"] + ref_counts["full_dup"] + ref_counts["frag"]

        full_single_proportion = float(sample_counts["full_sgl"]) / float(total_single_ref)
        full_single_proportions[1].append(full_single_proportion)

        dupe_proportion = float(sample_counts["full_dup"]) / float(total_single_ref)
        dupe_proportions[1].append(dupe_proportion)

        frag_proportion = float(sample_counts["frag"]) / float(total_single_ref)
        frag_proportions[1].append(frag_proportion)

        multicopy_proportion = float(sample_counts["dup_cnt"]) / float(ref_counts["dup_cnt"])
        multicopy_proportions[1].append(multicopy_proportion)

    print(full_single_proportions)
    print(dupe_proportions)
    print(frag_proportions)
    print(multicopy_proportions)

    width = 0.8/float(len(paths_diploid_only))  # the width of the bars

    fig, ax = pyplot.subplots()

    DIPLOID_ONLY = 0
    WITH_UNPHASED = 1

    offset = -(float(len(paths_diploid_only) - 1)/2.0)*width
    for i in range(len(full_single_proportions[0])):
        x = offset + float(i)*width
        ax.text(x+width*0.1, 0.05, labels[i], rotation="vertical", color="white", ha="center")
        ax.bar(x,
               full_single_proportions[WITH_UNPHASED][i],
               width,
               label=labels[i],
               color=colors[i],
               alpha=0.7)
        ax.bar(x,
               full_single_proportions[DIPLOID_ONLY][i],
               width,
               label=labels[i],
               color=colors[i])
        ax.bar(x,
               full_single_proportions[WITH_UNPHASED][i],
               width,
               label=labels[i],
               color=(0,0,0,0),
               edgecolor="gray",
               )

    offset += 1
    for i in range(len(full_single_proportions[0])):
        x = offset + float(i)*width
        ax.text(x+width*0.1, 0.05, labels[i], rotation="vertical", color="white", ha="center")
        ax.bar(x,
               multicopy_proportions[WITH_UNPHASED][i],
               width,
               label=labels[i],
               color=colors[i],
               alpha=0.7)
        ax.bar(x,
               multicopy_proportions[DIPLOID_ONLY][i],
               width,
               label=labels[i],
               color=colors[i])
        ax.bar(x,
               multicopy_proportions[WITH_UNPHASED][i],
               width,
               label=labels[i],
               color=(0,0,0,0),
               edgecolor="gray",
               )

    custom_lines = list()
    custom_lines.append(Line2D([0], [0], color=[0,0,0,0.5], lw=4))
    custom_lines.append(Line2D([0], [0], color=[0,0,0,0.3], lw=4))

    ax.legend(custom_lines, ["Phased","Unphased"], bbox_to_anchor=(1.05, 1.025))

    x_tick_labels = ["Full single copy genes", "Multi copy genes"]
    x_ticks = [0,1]

    minor_ticks = matplotlib.ticker.MultipleLocator(base=0.05)
    major_ticks = matplotlib.ticker.MultipleLocator(base=0.1)
    ax.yaxis.set_minor_locator(minor_ticks)
    ax.yaxis.set_major_locator(major_ticks)

    ax.set_xticklabels(x_tick_labels)
    ax.set_xticks(x_ticks)

    # ax.set_ylim([0,1.15])

    ax.axhline(1, linestyle="-", linewidth=0.5, zorder=-1)
    ax.axhline(0.9, linestyle="--", linewidth=0.5, zorder=-1)
    ax.axhline(0.8, linestyle="--", linewidth=0.5, zorder=-1)
    ax.axhline(0.7, linestyle="--", linewidth=0.5, zorder=-1)
    ax.axhline(0.6, linestyle="--", linewidth=0.5, zorder=-1)
    ax.axhline(0.5, linestyle="--", linewidth=0.5, zorder=-1)

    ax.set_title("HG002 gene completeness (asmgene)")
    ax.set_ylabel("Proportion matching (>99%)")

    fig.set_size_inches(4 + len(paths_diploid_only)*0.5,4)

    fig.tight_layout()
    pyplot.savefig("asmgene_comparison.png", dpi=200)
    pyplot.savefig("asmgene_comparison.pdf", dpi=200)

    pyplot.show()


if __name__ == "__main__":
    main()
