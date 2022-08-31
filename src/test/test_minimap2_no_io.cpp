#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include <vector>
#include <string>

#include "Filesystem.hpp"
#include "minimap.h"
#include "kseq.h"

using ghc::filesystem::path;

#include <exception>
#include <iostream>
#include <fstream>

using std::cerr;
using std::string;
using std::vector;
using std::ifstream;
using std::runtime_error;

void map(
        mm_idxopt_t index_options,
        mm_mapopt_t map_options,
        const vector<string>& targets,
        const vector<string>& target_names,
        const vector<string>& queries,
        const vector<string>& query_names){

    std::vector<const char*> c_targets;
    std::vector<const char*> c_names;

    c_targets.reserve(targets.size());
    c_names.reserve(targets.size());

    for(size_t i=0; i<targets.size(); i++) {
        c_targets.push_back(targets[i].c_str());
        c_names.push_back(target_names[i].c_str());
    }

    mm_verbose = 3;

    mm_idx_t *mi = mm_idx_str(
            index_options.w,
            index_options.k,
            int(0),
            index_options.bucket_bits,
            int(targets.size()),
            c_targets.data(),
            c_names.data()
    );

    mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
    for (size_t i=0; i<queries.size(); i++) {
        auto& query = queries[i];
        auto& name = query_names[i];

        int n_reg;
        mm_reg1_t *reg;

        mm_mapopt_update(&map_options, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
        reg = mm_map(mi, query.size(), query.c_str(), &n_reg, tbuf, &map_options, name.c_str()); // get all hits for the query

        cerr << '\n';
        cerr << "id" << ' ' << reg->id << '\n';
        cerr << "n_reg" << ' ' << n_reg << '\n';
        cerr << "score" << ' ' << reg->score << '\n';
        cerr << "cnt" << ' ' << reg->cnt << '\n';
        cerr << "mlen" << ' ' << reg->mlen << '\n';
        cerr << "blen" << ' ' << reg->blen << '\n';

        for (int j = 0; j < n_reg; ++j) { // traverse hits and print them out
            mm_reg1_t *r2 = &reg[j];

            assert(r2->p); // with MM_F_CIGAR, this should not be NULL
            printf("%s\t%d\t%d\t%d\t%c\t", name.c_str(), int(query.size()), r2->qs, r2->qe, "+-"[r2->rev]);
            printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\tcg:Z:", mi->seq[r2->rid].name, mi->seq[r2->rid].len, r2->rs, r2->re, r2->mlen, r2->blen, r2->mapq);
            for (int k = 0; k < r2->p->n_cigar; ++k) // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
                printf("%d%c", r2->p->cigar[k] >> 4, MM_CIGAR_STR[r2->p->cigar[k] & 0xf]);
            putchar('\n');
            free(r2->p);
        }
        free(reg);
    }
    mm_tbuf_destroy(tbuf);
    mm_idx_destroy(mi);
}


void read_mono_fasta(path fasta_path, string& seq){
    if (not (fasta_path.extension() == ".fa" or fasta_path.extension() == ".fasta")){
        throw runtime_error("ERROR: cannot read non fasta file: " + fasta_path.string());
    }

    ifstream file(fasta_path);

    string line;
    size_t n_header = 0;
    while (getline(file, line)){
        if (line[0] == '>'){
            n_header++;
            if (n_header > 1){
                throw runtime_error("ERROR: can't test minimap with multi-sequence fasta");
            }

            continue;
        }
        else {
            seq += line;
        }
    }
}


int main(int argc, char *argv[])
{
    if (argc < 3) {
        fprintf(stderr, "Usage: minimap2-lite <target.fa> <query.fa>\n");
        return 1;
    }

    string a;
    string b;

    read_mono_fasta(argv[1], a);
    read_mono_fasta(argv[2], b);

    const vector<string>& targets = {a};
    const vector<string>& target_names = {"ref"};
    const vector<string>& queries = {b};
    const vector<string>& query_names = {"query"};

    mm_idxopt_t iopt;
    mm_mapopt_t mopt;

    mm_verbose = 3; // disable message output to stderr
    mm_set_opt(0, &iopt, &mopt);
    mm_set_opt("asm20", &iopt, &mopt);

    iopt.k = 21;
    mopt.flag |= MM_F_CIGAR; // perform alignment
    mopt.flag |= MM_F_EQX;

    cerr << "k=" << iopt.k << '\n';
    cerr << "min_mid_occ=" << mopt.min_mid_occ << '\n';
    cerr << "max_mid_occ=" << mopt.max_mid_occ << '\n';
    cerr << "mid_occ=" << mopt.mid_occ << '\n';
    cerr << "mid_occ_frac=" << mopt.mid_occ_frac << '\n';

    map(iopt, mopt, targets, target_names, queries, query_names);
}