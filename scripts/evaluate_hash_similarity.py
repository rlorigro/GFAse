from collections import defaultdict
from modules.IterativeHistogram import *
import numpy


def get_lengths_from_gfa(gfa_path):
    lengths = dict()

    with open(gfa_path,'r') as file:
        for l,line in enumerate(file):
            if line.startswith("S"):
                print(line[:100])

                n = 0
                name = ""
                length = 0
                for i,c in enumerate(line):
                    if c.isspace():
                        n += 1
                    elif n == 1:
                        name += c
                    elif n == 2:
                        length += 1

                print(name, length)
                lengths[name] = length

    return lengths


class Match:
    def __init__(self, name, score, total_hashes, similarity):
        self.name = name
        self.score = score
        self.total_hashes = total_hashes
        self.similarity = similarity

    def __str__(self):
        return '\t'.join([str(self.name), str(self.score), str(self.total_hashes), str(self.similarity)])


def read_hash_similarities(hash_path):
    matches_per_name = defaultdict(list)

    with open(hash_path, 'r') as file:
        for l,line in enumerate(file):
            if l == 0:
                continue

            name, other_name, score, total_hashes, similarity = line.strip().split(',')
            score = int(score)
            total_hashes = int(total_hashes)
            similarity = float(similarity)

            match = Match(other_name, score, total_hashes, similarity)

            matches_per_name[name].append(match)

    return matches_per_name


def is_shasta_bubble(name, other_name):
    success = False

    if name.startswith("PR") and other_name.startswith("PR"):
        a = name.split('.')
        b = other_name.split('.')

        for i in range(len(a) - 1):
            if a[i] != b[i]:
                break

            if int(a[-1]) == 1 - int(b[-1]):
                success = True

    return success


def get_confusion_matrix(matches_per_name, certainty_threshold, size_threshold, lengths):
    confusion_matrix = numpy.zeros([2,2])
    length_matrix = numpy.zeros([2,2])

    for name,matches in matches_per_name.items():
        # print(name, "---")

        max_score = 0
        max_name = None
        total_hit_count = 0
        for m,match in enumerate(matches):
            # print(match)
            if match.score > max_score:
                max_score = match.score
                max_name = match.name

            # if m > 4:
            #     break

            total_hit_count += match.score

        certainty = float(max_score) / float(total_hit_count)

        success = certainty > certainty_threshold and max_score > size_threshold
        is_true = is_shasta_bubble(name, max_name)

        # print("RESULT: ", name, max_name, max_score, certainty, success, is_true)

        l = lengths.get(name)

        confusion_matrix[int(is_true)][int(success)] += 1
        length_matrix[int(is_true)][int(success)] += l

    return confusion_matrix, length_matrix


def get_precision(confusion_matrix):
    return confusion_matrix[1][1] / (confusion_matrix[1][1] + confusion_matrix[0][1] + 1e-12)


def get_recall(confusion_matrix):
    return confusion_matrix[1][1] / (confusion_matrix[1][1] + confusion_matrix[1][0] + 1e-12)


def main():
    gfa_path = "/home/ryan/data/test_gfase/paolo_ul_guppy6_run14/run14_uul_test_subset.gfa"
    hash_path = "/home/ryan/code/GFAse/build/test_minhash/overlaps.csv"

    lengths = get_lengths_from_gfa(gfa_path)

    true_positives = 0
    false_positives = 0
    true_negatives = 0
    false_negatives = 0

    matches_per_name = read_hash_similarities(hash_path)

    for i in range(10):
        threshold = float(i)/10.0
        confusion_matrix, length_matrix = get_confusion_matrix(
            matches_per_name=matches_per_name,
            certainty_threshold=threshold,
            size_threshold=0,
            lengths=lengths)

        print(confusion_matrix)

        precision = get_precision(confusion_matrix)
        recall = get_recall(confusion_matrix)

        print(length_matrix)

        length_precision = get_precision(length_matrix)
        length_recall = get_recall(length_matrix)

        print(threshold, precision, recall, length_precision, length_recall)


if __name__ == "__main__":
    main()
