from modules.Paf import *
from collections import defaultdict
from matplotlib import pyplot
import numpy
import sys
import os



def get_alignments_by_contig(maternal_ref_align_path, paternal_ref_align_path):
    alignments_by_contig = defaultdict(lambda: [list(), list()])

    with open(maternal_ref_align_path, 'r') as file:
        for l,line in enumerate(file):
            paf_element = PafElement(line, store_tags=True)

            name = paf_element.get_query_name()

            alignments_by_contig[name][1].append(paf_element)

    with open(paternal_ref_align_path, 'r') as file:
        for l,line in enumerate(file):
            paf_element = PafElement(line, store_tags=True)

            name = paf_element.get_query_name()

            alignments_by_contig[name][0].append(paf_element)

    return alignments_by_contig


def plot_dual_alignment(name, paternal_elements, maternal_elements):
    query_length = 0

    if len(paternal_elements) > 0:
        query_length = paternal_elements[0].get_query_length()

    if len(maternal_elements) > 0:
        query_length = maternal_elements[0].get_query_length()

    mat_mismatches = numpy.zeros(query_length)
    pat_mismatches = numpy.zeros(query_length)
    mat_inserts = numpy.zeros(query_length)
    pat_inserts = numpy.zeros(query_length)
    mat_deletes = numpy.zeros(query_length)
    pat_deletes = numpy.zeros(query_length)

    fig,axes = pyplot.subplots(nrows=3, sharex=True)

    print(name)

    print("--- PATERNAL ELEMENTS ---")
    for element in sorted(paternal_elements, key=lambda x: x.get_query_start()):
        for coord,operation,length in iterate_cigar(element):
            if operation == 'X':
                pat_mismatches[coord - 1] += length
            if operation == 'I':
                pat_inserts[coord - 1] += length
            if operation == 'D':
                pat_deletes[coord - 1] += length

        print(element)
        print()

    print("--- MATERNAL ELEMENTS ---")
    for element in sorted(maternal_elements, key=lambda x: x.get_query_start()):
        for coord,operation,length in iterate_cigar(element):
            if operation == 'X':
                mat_mismatches[coord - 1] += length
            if operation == 'I':
                mat_inserts[coord - 1] += length
            if operation == 'D':
                mat_deletes[coord - 1] += length

        print(element)
        print()

    print()

    window_size = 1000

    x = numpy.arange(window_size/2, query_length - window_size/2 + 1)

    pat_mismatches_binned = numpy.convolve(pat_mismatches, numpy.ones(window_size), 'valid')
    mat_mismatches_binned = numpy.convolve(mat_mismatches, numpy.ones(window_size), 'valid')
    pat_inserts_binned = numpy.convolve(pat_inserts, numpy.ones(window_size), 'valid')
    mat_inserts_binned = numpy.convolve(mat_inserts, numpy.ones(window_size), 'valid')
    pat_deletes_binned = numpy.convolve(pat_deletes, numpy.ones(window_size), 'valid')
    mat_deletes_binned = numpy.convolve(mat_deletes, numpy.ones(window_size), 'valid')

    max_mismatches = max(numpy.max(pat_mismatches_binned), numpy.max(mat_mismatches_binned))
    max_inserts = max(numpy.max(pat_inserts_binned), numpy.max(mat_inserts_binned))
    max_deletes = max(numpy.max(pat_deletes_binned), numpy.max(mat_deletes_binned))

    axes[0].set_ylim([-max_mismatches,max_mismatches])
    axes[1].set_ylim([-max_inserts,max_inserts])
    axes[2].set_ylim([-max_deletes,max_deletes])

    mat_mismatches_binned = -mat_mismatches_binned
    mat_inserts_binned = -mat_inserts_binned
    mat_deletes_binned = -mat_deletes_binned

    axes[0].plot(x, pat_mismatches_binned, linewidth=0.5, color="C0")
    axes[0].plot(x, mat_mismatches_binned, linewidth=0.5, color="C1")
    axes[1].plot(x, pat_inserts_binned, linewidth=0.5, color="C0")
    axes[1].plot(x, mat_inserts_binned, linewidth=0.5, color="C1")
    axes[2].plot(x, pat_deletes_binned, linewidth=0.5, color="C0")
    axes[2].plot(x, mat_deletes_binned, linewidth=0.5, color="C1")

    # Too costly to plot this way
    # axes[0].fill_between(x, pat_mismatches_binned, color="C0", alpha=0.3)
    # axes[0].fill_between(x, mat_mismatches_binned, color="C1", alpha=0.3)
    # axes[1].fill_between(x, pat_inserts_binned, color="C0", alpha=0.3)
    # axes[1].fill_between(x, mat_inserts_binned, color="C1", alpha=0.3)
    # axes[2].fill_between(x, pat_deletes_binned, color="C0", alpha=0.3)
    # axes[2].fill_between(x, mat_deletes_binned, color="C1", alpha=0.3)

    fig.set_size_inches(16,9)

    axes[0].set_title(name)

    axes[0].set_ylabel("Mismatches")
    axes[1].set_ylabel("Inserts")
    axes[2].set_ylabel("Deletes")

    pyplot.show()
    pyplot.close()


def assign_phase_by_alignment(name, paternal_elements, maternal_elements):
    query_length = 0

    if len(paternal_elements) > 0:
        query_length = paternal_elements[0].get_query_length()
    elif len(maternal_elements) > 0:
        query_length = maternal_elements[0].get_query_length()

    pat_contigs = set()
    mat_contigs = set()

    mat_matches = 0
    pat_matches = 0
    mat_mismatches = 0
    pat_mismatches = 0
    mat_inserts = 0
    pat_inserts = 0
    mat_deletes = 0
    pat_deletes = 0

    # print(name)

    max_indel_length = 50

    # print("--- PATERNAL ELEMENTS ---")
    mat_longest_alignment = 0
    mat_longest_contig = None
    pat_longest_alignment = 0
    pat_longest_contig = None

    for element in sorted(paternal_elements, key=lambda x: x.get_query_start()):
        l = element.get_ref_length()
        if l > pat_longest_alignment:
            pat_longest_alignment = l
            pat_longest_contig = element.get_ref_name()

        pat_contigs.add(element.get_ref_name())
        for coord,operation,length in iterate_cigar(element):
            if operation == '=':
                pat_matches += length
            if operation == 'X':
                pat_mismatches += length
            if operation == 'I' and length < max_indel_length:
                pat_inserts += length
            if operation == 'D' and length < max_indel_length:
                pat_deletes += length

        # print(element)
        # print()

    # print("--- MATERNAL ELEMENTS ---")
    for element in sorted(maternal_elements, key=lambda x: x.get_query_start()):
        l = element.get_ref_length()
        if l > mat_longest_alignment:
            mat_longest_alignment = l
            mat_longest_contig = element.get_ref_name()

        mat_contigs.add(element.get_ref_name())
        for coord,operation,length in iterate_cigar(element):
            if operation == '=':
                mat_matches += length
            if operation == 'X':
                mat_mismatches += length
            if operation == 'I' and length < max_indel_length:
                mat_inserts += length
            if operation == 'D' and length < max_indel_length:
                mat_deletes += length

        # print(element)
        # print()

    # print()

    mat_non_matches = float(mat_mismatches + mat_inserts + mat_deletes) + 1e-9
    pat_non_matches = float(pat_mismatches + pat_inserts + pat_deletes) + 1e-9

    score = numpy.log2(mat_non_matches/pat_non_matches)

    ref_name = None

    if score > 0:
        ref_name = pat_longest_contig
    elif score < 0:
        ref_name = mat_longest_contig
    elif pat_longest_contig == mat_longest_contig:
        ref_name = pat_longest_contig
    else:
        # If the scores are equal, see which alignment is longer
        pat_length = sum([c.get_ref_length() for c in paternal_elements])
        mat_length = sum([c.get_ref_length() for c in maternal_elements])

        if pat_length > mat_length and pat_length > 0:
            ref_name = pat_longest_contig
        elif mat_length > pat_length and mat_length > 0:
            ref_name = mat_longest_contig

    # if any("X" in c for c in mat_contigs):
    #     print("%s,%d,%.4f,%.4f,%.4f" % (name, query_length, pat_non_matches, mat_non_matches, score))
    #     print("pat_contigs", pat_contigs)
    #     print("mat_contigs", mat_contigs)

    # if "Y" in pat_contigs:
    #     print("%s,%d,%.4f,%.4f,%.4f" % (name, query_length, pat_non_matches, mat_non_matches, score))
    #     print("pat_contigs", pat_contigs)
    #     print("mat_contigs", mat_contigs)

    # if score == 0 and len(mat_contigs) == 0:
    #     print("%s,%d,%.4f,%.4f,%.4f" % (name, query_length, pat_non_matches, mat_non_matches, score))
    #     print("pat_contigs", pat_contigs)
    #     print("mat_contigs", mat_contigs)
    #     print("ref_name", ref_name)

    # if name.endswith(".p") and pat_identity < mat_identity:
    #     plot_dual_alignment(name=name, paternal_elements=paternal_elements, maternal_elements=maternal_elements)
    #
    # if name.endswith(".m") and mat_identity < pat_identity:
    #     plot_dual_alignment(name=name, paternal_elements=paternal_elements, maternal_elements=maternal_elements)

    return score, ref_name


def main():
    maternal_ref_align_path = "/home/ryan/data/test_gfase/paolo_ul_run11/gfase_het_hom_unphased_binning/align/unphased_initial_VS_HG002_mat_cur_20211005.paf"
    paternal_ref_align_path = "/home/ryan/data/test_gfase/paolo_ul_run11/gfase_het_hom_unphased_binning/align/unphased_initial_VS_HG002_pat_cur_20211005.paf"
    bandage_labels_path = "/home/ryan/data/test_gfase/paolo_ul_run11/gfase_het_hom_unphased_binning/unphased_parental_counts.csv"

    output_path = bandage_labels_path.replace(".csv", "_aligned.csv")

    alignments_by_contig = get_alignments_by_contig(maternal_ref_align_path, paternal_ref_align_path)

    mat_total_length = 0
    pat_total_length = 0
    unphased_total_length = 0

    names = set()

    with open(bandage_labels_path, 'r') as file, open(output_path, 'w') as output_file:
        for l,line in enumerate(file):
            if l == 0:
                output_file.write(line.strip() + ",ref_name,ref_score,query_length\n")

            name = line.strip().split(',')[0]

            if name in names:
                exit("ERROR: duplicate name in CSV")

            names.add(name)

            if name in alignments_by_contig:
                paternal_elements, maternal_elements = alignments_by_contig[name]

    # for name,[paternal_elements, maternal_elements] in alignments_by_contig.items():
        # if not "37.67035665" in name:
        #     continue
        #
        # plot_dual_alignment(name=name, paternal_elements=paternal_elements, maternal_elements=maternal_elements)

                score, ref_name = assign_phase_by_alignment(
                    name=name,
                    paternal_elements=paternal_elements,
                    maternal_elements=maternal_elements)

                query_length = 0
                if len(paternal_elements) > 0:
                    query_length = paternal_elements[0].get_query_length()
                elif len(maternal_elements) > 0:
                    query_length = maternal_elements[0].get_query_length()
                else:
                    exit("ERROR: no entries for contig: " + name)

                threshold = numpy.log2(0.01)

                if score < -threshold:
                    mat_total_length += query_length
                elif score > threshold:
                    pat_total_length += query_length
                else:
                    unphased_total_length += query_length

                output_file.write(line.strip())
                output_file.write(',')
                output_file.write(str(ref_name))
                output_file.write(',')
                output_file.write("%.4f" % score)
                output_file.write(',')
                output_file.write(str(query_length))
                output_file.write('\n')

    print("mat_total_length", mat_total_length)
    print("pat_total_length", pat_total_length)
    print("unphased_total_length", unphased_total_length)


if __name__ == "__main__":
    main()
