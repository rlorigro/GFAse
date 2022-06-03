import sys


def main():
    path = "/home/labmate/data/test/gfase/verkko/hg002_production_beta2_deepcns/assembly.homopolymer-compressed.gfa"
    output_path = '.'.join(path.split('.')[:-1]) + ".fasta"

    print(output_path)

    with open(path, 'r') as file, open(output_path,'w') as output_file:
        for l,line in enumerate(file):
            if line[0] == 'S' and line[1].isspace():
                data = line.strip().split()
                name = data[1]
                sequence = data[2]

                output_file.write('>')
                output_file.write(name)
                output_file.write('\n')
                output_file.write(sequence)
                output_file.write('\n')


if __name__ == "__main__":
    main()
