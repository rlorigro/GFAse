from modules.Paf import *
import argparse


def locate_cigar_events(paf_path, indel_threshold):

    n_inserts = 0
    n_deletes = 0

    reversal_char = ['+','-']

    print(','.join(["line", "ref_name", "ref_coord", "query_name", "query_coord", "reversal", "operation", "length", "map_quality"]))

    with open(paf_path, 'r') as file:
        for l,line in enumerate(file):
            element = PafElement(line, store_tags=True)

            for ref_coord,query_coord,operation,length in iterate_cigar(element):
                if length > indel_threshold:
                    if operation == 'I':
                        n_inserts += 1
                        print(','.join(list(map(str,[l,element.get_ref_name(), ref_coord, element.get_query_name(), query_coord, reversal_char[element.get_reversal()], operation, length, element.get_map_quality()]))))
                    if operation == 'D':
                        n_deletes += 1
                        print(','.join(list(map(str,[l,element.get_ref_name(), ref_coord, element.get_query_name(), query_coord, reversal_char[element.get_reversal()], operation, length, element.get_map_quality()]))))

    # print("n_inserts", n_inserts)
    # print("n_deletes", n_deletes)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        required=True,
        type=str,
        help="Input file of PAF alignment, containing cigar strings"
    )

    parser.add_argument(
        "-l",
        required=True,
        type=int,
        help="Minimum length of indel cigar operation (I/D) to report"
    )

    args = parser.parse_args()

    locate_cigar_events(paf_path=args.i, indel_threshold=args.l)

