#!/usr/bin/env python3
import subprocess
import argparse
import sys
import os


def align(ref_path, query_path, n_threads, as_bam=False):
    if not ref_path.endswith(".fasta") and not ref_path.endswith(".fa"):
        exit("ERROR: input file is not in FASTA format: " + ref_path)

    if not query_path.endswith(".fasta") and not query_path.endswith(".fa"):
        exit("ERROR: input file is not in FASTA format: " + ref_path)

    ref_prefix = '_'.join(os.path.basename(ref_path).split('.')[:-1])
    query_prefix = '_'.join(os.path.basename(query_path).split('.')[:-1])

    output_arg = "-c"
    output_extension = ".paf"

    if as_bam:
        output_arg = "-a"
        output_extension = ".sam"

    output_name = query_prefix + '_' + "VS" + '_' + ref_prefix + output_extension

    command = ["minimap2",
               "-x", "asm20",
               "-t", str(n_threads),
               "-K", "10g",
               "--secondary=no",
               "-n", "10",
               output_arg,
               "--eqx",
               ref_path,
               query_path]

    print("Running: " + ' '.join(command))
    print("Redirecting to: " + output_name)

    with open(output_name, 'w') as file:
        result = subprocess.run(command, check=True, stdout=file, stderr=sys.stderr, universal_newlines=True)
        print(result.stderr)

    if as_bam:
        input_name = output_name

        command = ["samtools", "view",
                   "-b",
                   "-h",
                   "-@", str(n_threads),
                   input_name,
                   ]

        output_name = input_name.replace(output_extension,".bam")

        print("Running: " + ' '.join(command))
        print("Redirecting to: " + output_name)

        with open(output_name, 'w') as file:
            result = subprocess.run(command, check=True, stdout=file, stderr=sys.stderr, universal_newlines=True)
            print(result.stderr)

        os.remove(input_name)

    return output_name


def cross_align(paternal_ref_path, maternal_ref_path, paternal_query_path, maternal_query_path, query_path, n_threads, as_bam):
    if n_threads is None:
        n_threads = os.cpu_count()

    if paternal_query_path is not None and paternal_ref_path is not None:
        output_name = align(ref_path=paternal_ref_path, query_path=paternal_query_path, n_threads=n_threads, as_bam=as_bam)
        print("Output file name: " + output_name)

    if paternal_query_path is not None and maternal_ref_path is not None:
        output_name = align(ref_path=maternal_ref_path, query_path=paternal_query_path, n_threads=n_threads, as_bam=as_bam)
        print("Output file name: " + output_name)

    if maternal_query_path is not None and paternal_ref_path is not None:
        output_name = align(ref_path=paternal_ref_path, query_path=maternal_query_path, n_threads=n_threads, as_bam=as_bam)
        print("Output file name: " + output_name)

    if maternal_query_path is not None and maternal_ref_path is not None:
        output_name = align(ref_path=maternal_ref_path, query_path=maternal_query_path, n_threads=n_threads, as_bam=as_bam)
        print("Output file name: " + output_name)

    if query_path is not None and paternal_ref_path is not None:
        output_name = align(ref_path=paternal_ref_path, query_path=query_path, n_threads=n_threads, as_bam=as_bam)
        print("Output file name: " + output_name)

    if query_path is not None and maternal_ref_path is not None:
        output_name = align(ref_path=maternal_ref_path, query_path=query_path, n_threads=n_threads, as_bam=as_bam)
        print("Output file name: " + output_name)

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--ref_pat",
        required=True,
        type=str
    )

    parser.add_argument(
        "--ref_mat",
        required=True,
        type=str
    )

    parser.add_argument(
        "--query_pat",
        required=False,
        default=None,
        type=str
    )

    parser.add_argument(
        "--query_mat",
        required=False,
        default=None,
        type=str
    )

    parser.add_argument(
        "--query",
        required=False,
        default=None,
        type=str
    )

    parser.add_argument(
        "--threads",
        required=False,
        default=None,
        type=int
    )

    parser.add_argument(
        "--as_bam",
        required=False,
        action="store_true",
    )

    args = parser.parse_args()

    cross_align(paternal_ref_path=args.ref_pat,
                maternal_ref_path=args.ref_mat,
                paternal_query_path=args.query_pat,
                maternal_query_path=args.query_mat,
                query_path=args.query_unphased,
                n_threads=args.threads,
                as_bam=args.as_bam)


