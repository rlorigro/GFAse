import argparse
import os


def strip_gfa_sequence(gfa_path):
    output_path = '.'.join(gfa_path.split('.')[:-1]) + "_no-sequence.gfa"

    print(output_path)

    with open(gfa_path, 'r') as file, open(output_path, 'w') as output_file:
        for line in file:
            if line.startswith("S"):
                line = line.split()
                line[2] = ''
                line = '\t'.join(line) + '\n'

                output_file.write(line)

            else:
                output_file.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        required=True,
        type=str
    )

    args = parser.parse_args()

    strip_gfa_sequence(args.i)
