from collections import defaultdict
from matplotlib import pyplot
import numpy
import sys
import os


def parse_cigar_as_tuples(cigar_string):
    operations = list()

    length_string = ""
    prev_is_numeric = True

    for i,c in enumerate(cigar_string):
        if c.isalpha() or c == '=':
            if i > 0 and not prev_is_numeric:
                exit("ERROR: cigar string contains impossible sequence of numeric and alphabetic characters")

            operations.append((c, int(length_string)))
            length_string = ""
            prev_is_numeric = False

        if c.isnumeric():
            length_string += c
            prev_is_numeric = True

    return operations


def is_reference_move(cigar_type):
    if cigar_type == 'M':
        return True
    elif cigar_type == 'I':
        return False
    elif cigar_type == 'D':
        return True
    elif cigar_type == '=':
        return True
    elif cigar_type == 'X':
        return True
    elif cigar_type == 'S':
        return False
    elif cigar_type == 'H':
        return False
    else:
        exit("ERROR: unrecognized cigar type: " + cigar_type)


def is_query_move(cigar_type):
    if cigar_type == 'M':
        return True
    elif cigar_type == 'I':
        return True
    elif cigar_type == 'D':
        return False
    elif cigar_type == '=':
        return True
    elif cigar_type == 'X':
        return True
    elif cigar_type == 'S':
        return True
    elif cigar_type == 'H':
        return False
    else:
        exit("ERROR: unrecognized cigar type: " + cigar_type)


def get_ref_alignment_length(cigar_operations):
    l = 0
    for item in cigar_operations:
        l += item[1]*is_reference_move(item[0])

    return l


class PafElement:
    def __init__(self, paf_line, store_tags=False):
        if store_tags:
            self.tokens = paf_line.strip().split()
        else:
            self.tokens = paf_line.strip().split()[:12]

        self.tokens[1] = int(self.tokens[1])
        self.tokens[6] = int(self.tokens[6])
        self.tokens[2] = int(self.tokens[2])
        self.tokens[3] = int(self.tokens[3])
        self.tokens[7] = int(self.tokens[7])
        self.tokens[8] = int(self.tokens[8])
        self.tokens[11] = int(self.tokens[11])
        self.tokens[4] = (self.tokens[4] == '-')

        self.cigar_index = None
        self.cigar = None

    def __str__(self):
        s = "query_name: " + str(self.get_query_name()) + '\n' + \
            "ref_name: " + str(self.get_ref_name()) + '\n' + \
            "query_length: " + str(self.get_query_length()) + '\n' + \
            "ref_length: " + str(self.get_ref_length()) + '\n' + \
            "query_start: " + str(self.get_query_start()) + '\n' + \
            "query_stop: " + str(self.get_query_stop()) + '\n' + \
            "ref_start: " + str(self.get_ref_start()) + '\n' + \
            "ref_stop: " + str(self.get_ref_stop()) + '\n' + \
            "map_quality: " + str(self.get_map_quality()) + '\n' + \
            "reversal:" + str(self.get_reversal())

        return s

    def get_data_by_column(self, c):
        data = self.tokens[c]

        # Parse differently if index indicates a tag field
        if c > 11:
            data = data.split(":")

            if len(data) < 2:
                exit("ERROR: specified a tag that has no colon delimiter")

            # Only keep what follows the last ":" delimiter
            data = data[-1]

        return data

    def get_query_name(self):
        return self.tokens[0]

    def get_ref_name(self):
        return self.tokens[5]

    def get_query_length(self):
        return self.tokens[1]

    def get_ref_length(self):
        return self.tokens[6]

    def get_query_start(self):
        return self.tokens[2]

    def get_query_stop(self):
        return self.tokens[3]

    def get_ref_start(self):
        return self.tokens[7]

    def get_ref_stop(self):
        return self.tokens[8]

    def get_map_quality(self):
        return self.tokens[11]

    def get_reversal(self):
        return self.tokens[4]

    def get_cigar(self):
        if self.cigar_index is None:
            self.find_cigar_index()

        if self.cigar is None:
            self.cigar = parse_cigar_as_tuples(self.tokens[self.cigar_index].split("cg:Z:")[-1])
            self.tokens[self.cigar_index] = ""

        return self.cigar

    def find_cigar_index(self):
        self.cigar_index = self.find_tag_index("cg:Z:")

    def find_tag_index(self, tag_prefix):
        for t,token in enumerate(self.tokens[12:]):
            if token.startswith(tag_prefix):
                return 12 + t

    def find_tag(self, tag_substring):
        for token in self.tokens[12:]:
            if token.startswith(tag_substring):
                return token


def iterate_cigar(paf_element):
    coord = paf_element.get_query_start()

    for operation,length in paf_element.get_cigar():
        yield coord,operation,length

        if is_query_move(operation):
            coord += length


def main():
    maternal_ref_align_path = "/home/ryan/data/test_gfase/paolo_ul_run9/phased_k31_double_coverage/align/hprc_hg002/maternal_VS_HG002_mat_cur_20211005.paf"
    paternal_ref_align_path = "/home/ryan/data/test_gfase/paolo_ul_run9/phased_k31_double_coverage/align/hprc_hg002/maternal_VS_HG002_pat_cur_20211005.paf"

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

    for name,[paternal_elements, maternal_elements] in alignments_by_contig.items():
        # if not "0.1.m" in name:
        #     continue

        query_length = 0

        if len(paternal_elements) > 0:
            query_length = paternal_elements[0].get_query_length()

        if len(maternal_elements) > 0:
            query_length = maternal_elements[0].get_query_length()

        if query_length < 1_000_000 or query_length > 30_000_000:
            continue

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


if __name__ == "__main__":
    main()
