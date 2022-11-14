from matplotlib.lines import Line2D
from matplotlib import pyplot
import argparse
import os.path
import pandas
import seaborn


# def main(input_paths, output_dir, color_indexes):
def main():

    # -----------------------------------------------------------------------------------------------------------------

    vertical_line_pos = 0

    paths = [
        "/home/ryan/data/test_gfase/terra/verkko_hg002_prod_v1_1/HG002_verkko_production_v1_1_hic_mc_manual.tsv",
        "/home/ryan/data/test_gfase/terra/verkko_hg002_prod_v1_1/HG002_verkko_production_v1_1_hic_mc_m3_manual.tsv",
        "/home/ryan/data/test_gfase/terra/verkko_hg002_prod_v1_1/HG002_verkko_production_v1_1_hic_mc_m5_manual.tsv",
        "/home/ryan/data/test_gfase/terra/verkko_hg002_prod_v1_1/HG002_verkko_production_v1_1_hic_mc_m7_manual.tsv",
    ]

    palette = seaborn.color_palette([
        "#3496A5",
        "#3496A5",
        "#3496A5",
        "#3496A5",
    ],
        len(paths))

    # seaborn.set_palette(palette)

    legend_colors = [
        "#3496A5",
        "#3496A5",
        "#3496A5",
        "#3496A5",
    ]

    legend = [
        "Verkko + GFAse 'production' m1 (HiC)",
        "Verkko + GFAse 'production' m3 (HiC)",
        "Verkko + GFAse 'production' m5 (HiC)",
        "Verkko + GFAse 'production' m7 (HiC)",
    ]

    # # -----------------------------------------------------------------------------------------------------------------
    # vertical_line_pos = 4.5
    #
    # paths = [
    #     "/home/ryan/data/test_gfase/terra/shasta_standard_2fc/shasta_hg002_standard_2fc_hic.tsv",
    #     "/home/ryan/data/test_gfase/terra/linked_read_experiments/HG002_shasta_v0.10.0_UUL_p50_porec_1fc.tsv",
    #     "/home/ryan/data/test_gfase/terra/parameter_search/p50.tsv",
    #     "/home/ryan/data/test_gfase/terra/shasta_r10_slow_mode/shasta_hg002_r10_slow_mode_porec_mc_phase.tsv",
    #     "/home/ryan/data/test_gfase/terra/shasta_r10_slow_mode/shasta_hg002_r10_slow_mode_mc_phase.tsv",
    #     "/home/ryan/data/test_gfase/terra/linked_read_experiments/HG002_trio_hifiasm_v0.14.1_HPRC.tsv",
    #     "/home/ryan/data/test_gfase/terra/linked_read_experiments/HG002_verkko_v1.1_full_coverage.tsv",
    #     "/home/ryan/data/test_gfase/terra/verkko_hg002_prod_v1_1/HG002_verkko_production_v1_1_porec_2fc_mc_manual.tsv",
    #     "/home/ryan/data/test_gfase/terra/verkko_hg002_prod_v1_1/HG002_verkko_production_v1_1_hic_mc_manual.tsv"
    # ]
    #
    # palette = seaborn.color_palette([
    #     "#3496A5",
    #     "#3496A5",
    #     "#3496A5",
    #     "#6CDB9F",
    #     "#6CDB9F",
    #     "#F27E7E",
    #     "#f6b26bff",
    #     "#f6b26bff",
    #     "#f6b26bff",
    # ],
    #     len(paths))
    #
    # # seaborn.set_palette(palette)
    #
    # legend_colors = [
    #     "#3496A5",
    #     "#3496A5",
    #     "#3496A5",
    #     "#6CDB9F",
    #     "#6CDB9F",
    #     "#F27E7E",
    #     "#f6b26bff",
    #     "#f6b26bff",
    #     "#f6b26bff",
    # ]
    #
    # legend = [
    #     "Shasta R9 Standard + GFAse (2fc HiC)",
    #     "Shasta R9 UUL + GFAse (PoreC)",
    #     "Shasta R9 UUL + GFAse (HiC)",
    #     "Shasta R10 + GFAse (PoreC)",
    #     "Shasta R10 + GFAse (HiC)",
    #     "Trio Hifiasm",
    #     "Verkko + GFAse 'full coverage' (HiC)",
    #     "Verkko + GFAse 'production' (PoreC)",
    #     "Verkko + GFAse 'production' (HiC)",
    # ]
    #
    # # -----------------------------------------------------------------------------------------------------------------
    #
    # vertical_line_pos = 6.5
    #
    # paths = [
    #     "/home/ryan/data/test_gfase/terra/linked_read_experiments/HG002_shasta_v0.10.0_UUL.tsv",
    #     "/home/ryan/data/test_gfase/terra/parameter_search/p20.tsv",
    #     "/home/ryan/data/test_gfase/terra/parameter_search/p30.tsv",
    #     "/home/ryan/data/test_gfase/terra/parameter_search/p40.tsv",
    #     "/home/ryan/data/test_gfase/terra/parameter_search/p50.tsv",
    #     "/home/ryan/data/test_gfase/terra/parameter_search/p60.tsv",
    #     "/home/ryan/data/test_gfase/terra/parameter_search/p70.tsv",
    #     "/home/ryan/data/test_gfase/terra/linked_read_experiments/HG002_trio_hifiasm_v0.14.1_HPRC.tsv",
    #     "/home/ryan/data/test_gfase/terra/linked_read_experiments/HG002_verkko_v1.1_full_coverage.tsv"
    # ]
    #
    # palette = seaborn.color_palette([
    #     "#7EF29D",
    #     "#6CDB9F",
    #     "#59C4A1",
    #     "#47ADA3",
    #     "#3496A5",
    #     "#227FA7",
    #     "#0F68A9",
    #     "#F27E7E",
    #     "#f6b26bff"],
    #     len(paths))
    #
    # # seaborn.set_palette(palette)
    #
    # legend_colors = [
    #     "#7EF29D",
    #     "#6CDB9F",
    #     "#59C4A1",
    #     "#47ADA3",
    #     "#3496A5",
    #     "#227FA7",
    #     "#0F68A9",
    #     "#F27E7E",
    #     "#f6b26bff"]
    #
    # legend = [
    #     "minLogP = 10",
    #     "minLogP = 20",
    #     "minLogP = 30",
    #     "minLogP = 40",
    #     "minLogP = 50",
    #     "minLogP = 60",
    #     "minLogP = 70",
    #     "Trio Hifiasm",
    #     "Verkko + GFAse"]

    # -----------------------------------------------------------------------------------------------------------------
    # paths = [
    #     "/home/ryan/data/test_gfase/terra/linked_read_experiments/HG002_shasta_v0.10.0_Standard.tsv",
    #     "/home/ryan/data/test_gfase/terra/linked_read_experiments/HG002_shasta_v0.10.0_UL.tsv",
    #     "/home/ryan/data/test_gfase/terra/linked_read_experiments/HG002_shasta_v0.10.0_UUL.tsv",
    #     "/home/ryan/data/test_gfase/terra/linked_read_experiments/HG002_trio_hifiasm_v0.14.1_HPRC.tsv",
    #     "/home/ryan/data/test_gfase/terra/linked_read_experiments/HG002_verkko_v1.1_full_coverage.tsv",
    #     "/home/ryan/data/test_gfase/terra/linked_read_experiments/HG01993_shasta_v0.10.0_UL.tsv",
    #     "/home/ryan/data/test_gfase/terra/linked_read_experiments/HG02132_shasta_v0.10.0_UL.tsv",
    #     "/home/ryan/data/test_gfase/terra/linked_read_experiments/HG02647_shasta_v0.10.0_UL.tsv",
    #     "/home/ryan/data/test_gfase/terra/linked_read_experiments/HG03669_shasta_v0.10.0_UL.tsv"]
    #
    # palette = seaborn.color_palette([
    #     "#6fa8dcff",
    #     "#6fa8dcff",
    #     "#6fa8dcff",
    #     "#f6b26bff",
    #     "#93c47dff",
    #     "#6fa8dcff",
    #     "#6fa8dcff",
    #     "#6fa8dcff",
    #     "#6fa8dcff"],
    #     9)
    #
    # # seaborn.set_palette(palette)
    #
    # legend_colors = [
    #     "#6fa8dcff",
    #     "#93c47dff",
    #     "#f6b26bff"]
    #
    # legend = [
    #     "Shasta + GFAse",
    #     "Verkko + GFAse",
    #     "Trio Hifiasm"]

    hamming_per_sample = pandas.DataFrame()
    switch_per_sample = pandas.DataFrame()
    covered_variants_per_sample = pandas.DataFrame()
    names = list()

    for p,path in enumerate(paths):
        name = os.path.basename(path).replace(".tsv","")
        names.append(name)

        df = pandas.read_csv(path, sep='\t', header=0).set_index("chromosome")

        switch = df["all_switch_rate"].rename(name)
        hamming = df["blockwise_hamming_rate"].rename(name)
        covered_variants = df["covered_variants"].rename(name)

        if p == 0:
            switch_per_sample = switch
            hamming_per_sample = hamming
            covered_variants_per_sample = covered_variants
        else:
            switch_per_sample = pandas.concat([switch_per_sample,switch], axis=1)
            hamming_per_sample = pandas.concat([hamming_per_sample,hamming], axis=1)
            covered_variants_per_sample = pandas.concat([covered_variants_per_sample,covered_variants], axis=1)

    # print(covered_variants_per_sample)
    # print(covered_variants_per_sample.keys())
    # print(covered_variants_per_sample.index)

    if "ALL" in hamming_per_sample:
        hamming_per_sample.drop("ALL")

    if "chrY" in hamming_per_sample:
        hamming_per_sample.drop("chrY")

    if "chrX" in hamming_per_sample:
        hamming_per_sample.drop("chrX")

    covered_variants_per_sample = covered_variants_per_sample.loc[["ALL"]]
    print(covered_variants_per_sample)

    if "ALL" in switch_per_sample:
        switch_per_sample.drop("ALL")

    if "chrY" in switch_per_sample:
        switch_per_sample.drop("chrY")

    if "chrX" in switch_per_sample:
        switch_per_sample.drop("chrX")

    print(covered_variants_per_sample)

    f,a = pyplot.subplots(nrows=3, gridspec_kw={'height_ratios': [2, 2, 1]})
    f.set_size_inches(9,12)

    custom_lines = list()
    for i in range(len(legend_colors)):
        custom_lines.append(Line2D([0], [0], color=legend_colors[i], lw=4))

    sub_index = 0
    seaborn.swarmplot(switch_per_sample, s=3, ax=a[sub_index], palette=palette)
    seaborn.boxplot(switch_per_sample, boxprops={'facecolor':'None'}, showfliers=False, ax=a[sub_index])
    a[sub_index].get_xaxis().set_visible(False)
    a[sub_index].set_ylabel("Switch rate")
    a[sub_index].axvline(vertical_line_pos, linestyle="--", linewidth=0.5)
    a[sub_index].legend(custom_lines, legend, bbox_to_anchor=(1.05, 1.025))

    sub_index = 1
    a[sub_index].set_yscale("log")
    a[sub_index].get_yaxis().set_major_formatter(pyplot.ScalarFormatter())
    seaborn.swarmplot(hamming_per_sample, s=3, ax=a[sub_index], palette=palette)
    seaborn.boxplot(hamming_per_sample, boxprops={'facecolor':'None'}, showfliers=False, ax=a[sub_index])
    a[sub_index].get_xaxis().set_visible(False)
    a[sub_index].set_ylabel("Hamming error")
    a[sub_index].axvline(vertical_line_pos, linestyle="--", linewidth=0.5)

    sub_index = 2
    seaborn.barplot(covered_variants_per_sample, ax=a[sub_index], errorbar=None, palette=palette)
    pyplot.xticks(fontsize=9, rotation=70, ha='right')
    a[sub_index].set_ylabel("Covered variants", fontsize=11)
    a[sub_index].axvline(vertical_line_pos, linestyle="--", linewidth=0.5)
    a[sub_index].set_xticklabels(legend)
    pyplot.tight_layout()

    pyplot.savefig("whatshap_comparison.png", dpi=200)
    pyplot.show()
    pyplot.close()


def parse_comma_separated_string(s):
    return s.strip().split(',')


def parse_comma_separated_int_string(s):
    if s is None:
        return []
    else:
        return list(map(int,s.strip().split(',')))


if __name__ == "__main__":
    main()

    # parser = argparse.ArgumentParser()
    #
    # parser.add_argument(
    #     "-i",
    #     required=True,
    #     type=str,
    #     help="Comma separated list of: input TSV file paths"
    # )
    #
    # parser.add_argument(
    #     "-o",
    #     required=True,
    #     type=str,
    #     help="Output directory"
    # )
    #
    # parser.add_argument(
    #     "-c","--colors",
    #     required=False,
    #     default=None,
    #     type=str,
    #     help="Comma separated list of color indexes to use for each item in the plot: "
    #          "0=red, "
    #          "1=orange, "
    #          "2=yellow, "
    #          "3=light green, "
    #          "4=green, "
    #          "5=green-blue, "
    #          "6=turquoise, "
    #          "7=blue, "
    #          "8=indigo, "
    #          "9=purple, "
    #          "10=pink, "
    # )
    #
    # args = parser.parse_args()
    #
    # args.i = parse_comma_separated_string(args.i)
    #
    # if args.colors is not None:
    #     args.colors = parse_comma_separated_int_string(args.colors)
    #
    #     if len(args.colors) != len(args.i):
    #         exit("ERROR: color list length doesn't match list of input files")
    #
    #     for c in args.colors:
    #         if 0 > c > 10:
    #             exit("ERROR: color argument must be between 0 and 10 (inclusive), see help for details")
    #
    # main(input_paths=args.i, output_dir=args.o, color_indexes=args.colors)
