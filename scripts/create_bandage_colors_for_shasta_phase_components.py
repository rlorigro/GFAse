#!/usr/bin/env python3
import argparse
import os
import random


def iterate_fasta(path):
    name = ""
    sequence = ""

    with open(path, 'r') as file:
        for l,line in enumerate(file):
            if line.startswith('>'):
                # Output previous sequence
                if l > 0:
                    yield name, sequence

                name = line.strip()[1:].split(' ')[0]
                sequence = ""

            else:
                sequence += line.strip()

    yield name, sequence


def iterate_gfa(path):
    with open(path, 'r') as file:
        for l,line in enumerate(file):
            if line[0] == 'S' and line[1].isspace():
                data = line.strip().split()
                name = data[1]
                sequence = data[2]

                yield name, sequence


def write_sequence_to_fasta(name, sequence, file):
    file.write('>')
    file.write(name)
    file.write('\n')
    file.write(sequence)
    file.write('\n')


def main(path, output_directory):

    if not len(output_directory) == 0:
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

    output_path = os.path.join(output_directory,"colored_components.csv")

    iterator = iterate_fasta(path)
    if path.endswith(".gfa"):
        iterator = iterate_gfa(path)

    # colors = [
    #     "#BC8F8F","#CD5C5C","#A52A2A","#B22222","#F08080","#800000","#8B0000","#FF0000","#FA8072","#FF6347",
    #     "#FF4500","#FF7F50","#FFA07A","#FF8C00","#FFA500","#DAA520","#B8860B","#FFD700","#F0E68C","#FFFACD",
    #     "#808000","#FFFF00","#6B8E23","#9ACD32","#556B2F","#ADFF2F","#7CFC00","#7FFF00","#8FBC8F","#228B22",
    #     "#32CD32","#90EE90","#98FB98","#006400","#008000","#00FF00","#3CB371","#2E8B57","#00FF7F","#00FA9A",
    #     "#66CDAA","#7FFFD4","#40E0D0","#20B2AA","#48D1CC","#2F4F4F","#AFEEEE","#008080","#008B8B","#00CED1",
    #     "#00FFFF","#00FFFF","#5F9EA0","#B0E0E6","#ADD8E6","#00BFFF","#87CEEB","#87CEFA","#4682B4","#1E90FF",
    #     "#708090","#778899","#B0C4DE","#6495ED","#4169E1","#191970","#000080","#00008B","#0000CD","#0000FF",
    #     "#483D8B","#6A5ACD","#7B68EE","#9370DB","#663399","#8A2BE2","#4B0082","#9932CC","#9400D3","#BA55D3",
    #     "#DDA0DD","#EE82EE","#800080","#8B008B","#FF00FF","#FF00FF","#DA70D6","#C71585","#FF1493","#FF69B4",
    #     "#DB7093","#DC143C","#FFC0CB","#FFB6C1","#BC8F8F","#CD5C5C","#A52A2A","#B22222","#F08080","#800000",
    #     "#8B0000","#FF0000"
    # ]

    colors = [
        "#800000",
        "#808000",
        "#3cb371",
        "#000080",
        "#ff0000",
        "#ff8c00",
        "#ffd700",
        "#7fff00",
        "#ba55d3",
        "#00ff7f",
        "#e9967a",
        "#00ffff",
        "#0000ff",
        "#ff00ff",
        "#1e90ff",
        "#eee8aa",
        "#dda0dd",
        "#ff1493",
    ]

    component_colors = dict()

    with open(output_path, 'w') as file:
        file.write("Name,Component,Color\n")

        for name,sequence in iterator:
            is_phased = name.startswith("PR")

            if is_phased:
                component = name.split('.')[-2]
                print(name, component)

                color = ""
                if component not in component_colors:
                    color = random.choice(colors)
                    component_colors[component] = color
                else:
                    color = component_colors[component]

                file.write(name)
                file.write(',')
                file.write(component)
                file.write(',')
                file.write(color)
                file.write('\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        required=True,
        type=str,
        help="Input fasta to be split"
    )

    parser.add_argument(
        "-o",
        required=True,
        type=str,
        help="Output directory"
    )

    args = parser.parse_args()

    main(path=args.i, output_directory=args.o)
