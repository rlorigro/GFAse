import subprocess
import argparse
import sys
import os


def align(ref_path, query_path, n_threads):
    if not ref_path.endswith(".fasta") and not ref_path.endswith(".fa"):
        exit("ERROR: input file is not in FASTA format: " + ref_path)

    if not query_path.endswith(".fasta") and not query_path.endswith(".fa"):
        exit("ERROR: input file is not in FASTA format: " + ref_path)

    ref_prefix = '_'.join(os.path.basename(ref_path).split('.')[:-1])
    query_prefix = '_'.join(os.path.basename(query_path).split('.')[:-1])
    output_name = query_prefix + '_' + "VS" + '_' + ref_prefix + ".paf"

    command = ["minimap2",
               "-x", "asm20",
               "-t", str(n_threads),
               "--secondary=no",
               "-n", "10",
               "-c",
               "--eqx",
               ref_path,
               query_path]

    print("Running: " + ' '.join(command))

    with open(output_name, 'w') as file:
        result = subprocess.run(command, check=True, stdout=file, stderr=sys.stderr, universal_newlines=True)
        print(result.stderr)

    return output_name


def cross_align(paternal_ref_path, maternal_ref_path, paternal_query_path, maternal_query_path, unphased_query_path, n_threads):
    if n_threads is None:
        n_threads = os.cpu_count()

    if paternal_query_path is not None and paternal_ref_path is not None:
        output_name = align(ref_path=paternal_ref_path, query_path=paternal_query_path, n_threads=n_threads)
        print("Output file name: " + output_name)

    if paternal_query_path is not None and maternal_ref_path is not None:
        output_name = align(ref_path=maternal_ref_path, query_path=paternal_query_path, n_threads=n_threads)
        print("Output file name: " + output_name)

    if maternal_query_path is not None and paternal_ref_path is not None:
        output_name = align(ref_path=paternal_ref_path, query_path=maternal_query_path, n_threads=n_threads)
        print("Output file name: " + output_name)

    if maternal_query_path is not None and maternal_ref_path is not None:
        output_name = align(ref_path=maternal_ref_path, query_path=maternal_query_path, n_threads=n_threads)
        print("Output file name: " + output_name)

    if unphased_query_path is not None and paternal_ref_path is not None:
        output_name = align(ref_path=paternal_ref_path, query_path=unphased_query_path, n_threads=n_threads)
        print("Output file name: " + output_name)

    if unphased_query_path is not None and maternal_ref_path is not None:
        output_name = align(ref_path=maternal_ref_path, query_path=unphased_query_path, n_threads=n_threads)
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
        "--query_unphased",
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

    args = parser.parse_args()

    cross_align(paternal_ref_path=args.ref_pat,
                maternal_ref_path=args.ref_mat,
                paternal_query_path=args.query_pat,
                maternal_query_path=args.query_mat,
                unphased_query_path=args.query_mat,
                n_threads=args.threads)


'''
python3 ../scripts/cross_align_phases.py \
--ref_pat [ref_pat].fasta \
--ref_mat [ref_mat].fasta \
--query_pat [query_pat].fasta \
--query_mat [query_mat].fasta

'''
