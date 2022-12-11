from collections import defaultdict


def count_runlengths(sequence, runlengths):
    prev_char = "_"
    l = 0

    # print(sequence[:100])

    for char in sequence:
        if char != prev_char:
            # print(prev_char,l)

            if l > 40:
                runlengths[prev_char][l] += 1

            prev_char = char
            l = 1

        else:
            l += 1


def main():
    path = "/home/ryan/data/human/reference/chm13.draft_v1.1.fasta"

    runlengths = defaultdict(lambda: defaultdict(int))

    chunks = list()

    with open(path, 'r') as file:
        for l,line in enumerate(file):
            if line[0] == '>':
                print(line)

                if (l > 0):
                    print(len(chunks))
                    sequence = ''.join(chunks)
                    count_runlengths(sequence, runlengths)
                    chunks = list()

            else:
                # Increment most recent length by size of line (without newline char)
                chunks.append(line.strip())

    for c in runlengths.keys():
        print(c)
        for l,count in sorted(runlengths[c].items(), key=lambda x: x[0], reverse=True):
            print(l,count)


if __name__ == "__main__":
    main()
