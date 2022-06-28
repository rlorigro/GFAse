import os

from matplotlib import pyplot
from collections import defaultdict
import argparse
import numpy
import math


def iterate_bubbles_csv(csv_path):
    with open(csv_path, 'r') as file:
        for l,line in enumerate(file):
            if l == 0:
                continue

            data = line.strip().split(',')

            name = data[0]
            bubble_id = int(data[1])
            phase = bool(int(data[2]))

            yield name, bubble_id, phase


def iterate_contacts_csv(csv_path):
    with open(csv_path, 'r') as file:
        for l,line in enumerate(file):
            data = line.strip().split(',')

            name_a = data[0]
            name_b = data[1]

            mapq_distribution = dict()

            for item in data[2].split(' '):
                data = item.split(':')
                mapq_distribution[int(data[0])] = int(data[1])

            yield name_a, name_b, mapq_distribution


def main(bubbles_path, contacts_path, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    intra_bubble_mapq = defaultdict(int)
    consistent_phase_mapq = defaultdict(int)
    inconsistent_phase_mapq = defaultdict(int)

    allowed_names = {}

    intra_bubble_edges = defaultdict(lambda: set())
    phases = [set(), set()]

    prev_bubble_id = -1
    bubble_names = [None,None]
    for name,bubble_id,phase in iterate_bubbles_csv(bubbles_path):
        print(name,bubble_id,phase)

        if bubble_id != prev_bubble_id:
            if len(allowed_names) == 0 or ((bubble_names[0] in allowed_names) and (bubble_names[1] in allowed_names)):
                intra_bubble_edges[bubble_names[0]].add(bubble_names[1])
                intra_bubble_edges[bubble_names[1]].add(bubble_names[0])

                phases[0].add(bubble_names[0])
                phases[1].add(bubble_names[1])

                bubble_names = [None,None]

        bubble_names[phase] = name
        prev_bubble_id = bubble_id

    output_csv_path = os.path.join(output_dir, "labeled_contacts.csv")

    max_frequency = 0
    with open(output_csv_path, 'w') as file:
        for name_a, name_b, mapq_distribution in iterate_contacts_csv(contacts_path):
            file.write("%s,%s," % (name_a, name_b))

            if name_b in intra_bubble_edges[name_a]:
                file.write("2,")
                for a,b in mapq_distribution.items():
                    file.write("%d:%d " % (a,b))
                    intra_bubble_mapq[a] += b

                    if b > max_frequency:
                        max_frequency = b

            else:
                a_in_0 = name_a in phases[0]
                a_in_1 = name_a in phases[1]
                b_in_0 = name_b in phases[0]
                b_in_1 = name_b in phases[1]

                if (a_in_0 != a_in_1) and (b_in_0 != b_in_1):
                    if a_in_0 == b_in_0:
                        file.write("0,")
                        for a,b in mapq_distribution.items():
                            file.write("%d:%d " % (a,b))
                            consistent_phase_mapq[a] += b

                            if b > max_frequency:
                                max_frequency = b

                    else:
                        file.write("1,")
                        for a,b in mapq_distribution.items():
                            file.write("%d:%d " % (a,b))
                            inconsistent_phase_mapq[a] += b

                            if b > max_frequency:
                                max_frequency = b

            file.write("\n")

    intra_bubble_histogram = numpy.zeros(61)
    inconsistent_phase_histogram = numpy.zeros(61)
    consistent_phase_histogram = numpy.zeros(61)

    # Mapq cant be greater than 60
    for a,b in intra_bubble_mapq.items():
        intra_bubble_histogram[a] += b

    for a,b in consistent_phase_mapq.items():
        consistent_phase_histogram[a] += b

    for a,b in inconsistent_phase_mapq.items():
        inconsistent_phase_histogram[a] += b

    # intra_bubble_histogram /= sum(intra_bubble_histogram)
    # inconsistent_phase_histogram /= sum(inconsistent_phase_histogram)
    # consistent_phase_histogram /= sum(consistent_phase_histogram)

    fig = pyplot.figure()
    axes = pyplot.axes()

    x = numpy.arange(0,61,1)
    print(x)
    print(intra_bubble_histogram)

    pyplot.plot(intra_bubble_histogram, color="C0", alpha=0.6, label="intra_bubble")
    pyplot.plot(consistent_phase_histogram, color="C1", alpha=0.6, label="consistent_phase")
    pyplot.plot(inconsistent_phase_histogram, color="C2", alpha=0.6, label="inconsistent_phase")
    axes.legend()

    axes.set_ylabel("Frequency (#)")
    axes.set_xlabel("Mapq")

    pyplot.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--bubbles",
        required=True,
        type=str,
        help="Input file of true_phases.csv"
    )

    parser.add_argument(
        "--contacts",
        required=True,
        type=str,
        help="Input file of contacts.csv containing mapq distributions for each contact"
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        type=str,
        help="Directory path where output will be written"
    )

    args = parser.parse_args()

    main(bubbles_path=args.bubbles, contacts_path=args.contacts, output_dir=args.output_dir)
