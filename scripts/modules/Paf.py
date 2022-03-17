from collections import defaultdict
from matplotlib import pyplot
import numpy


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
    cigar_operations = paf_element.get_cigar()

    ref_coord = paf_element.get_ref_start()
    query_coord = paf_element.get_query_start()

    if paf_element.get_reversal():
        cigar_operations = reversed(cigar_operations)
        ref_coord = paf_element.get_ref_stop()
        query_coord = paf_element.get_query_start()

    for operation,length in cigar_operations:
        yield ref_coord,query_coord,operation,length

        if is_reference_move(operation):
            ref_coord += length*(1-2*int(paf_element.get_reversal()))

        if is_query_move(operation):
            query_coord += length
