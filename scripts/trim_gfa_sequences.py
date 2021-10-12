#!python3
# -*- coding: utf-8 -*-

import sys


def main(file_path):
    length = 30

    with open(file_path) as file:
        for l,line in enumerate(file):
            if line[0] == "S":
                line_type, name, sequence = line.split('\t')[0:3]

                if len(sequence) > length*2:
                    sequence = sequence[0:length] + "GGGGGGGG" + sequence[-length:]

                sys.stdout.write('\t'.join([line_type, name, sequence]))
                sys.stdout.write('\n')
            else:
                sys.stdout.write(line)


if __name__ == "__main__":

    if len(sys.argv) != 2:
        exit("ERROR: need to provide 1 argument: path of input GFA")

    main(sys.argv[1])

