#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "minimap.h"
#include "kseq.h"

#include <iostream>
#include <string>

using std::string;
using std::cout;
using std::cerr;


KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;
    int n_threads = 1;

//    iopt.k = 19;
//    iopt.w = 10;
//    mopt.min_mid_occ = 50;
//    mopt.max_mid_occ = 500;
//    mopt.flag |= MM_F_RMQ;     // long indel alignment (asm20)
//    mopt.bw = 100000;
//    mopt.max_gap = 10000;
//    mopt.a = 1;
//    mopt.b = 4;
//    mopt.q = 6;
//    mopt.q2 = 26;
//    mopt.e = 2;
//    mopt.e2 = 1;
//    mopt.min_dp_max = 200;
//    mopt.zdrop = mopt.zdrop_inv = 200;
//    mopt.best_n = 50;

    mm_verbose = 3; // disable message output to stderr
    mm_set_opt(0, &iopt, &mopt);
    mm_set_opt("asm20", &iopt, &mopt);

    mopt.min_cnt = 10;

    iopt.k = 21;
    mopt.flag |= MM_F_CIGAR; // perform alignment
//    mopt.flag |= MM_F_EQX;

    cerr << "k=" << iopt.k << '\n';
    cerr << "min_mid_occ=" << mopt.min_mid_occ << '\n';
    cerr << "max_mid_occ=" << mopt.max_mid_occ << '\n';

    if (argc < 3) {
        fprintf(stderr, "Usage: minimap2-lite <target.fa> <query.fa>\n");
        return 1;
    }

    // open query file for reading; you may use your favorite FASTA/Q parser
    gzFile f = gzopen(argv[2], "r");
    assert(f);
    kseq_t *ks = kseq_init(f);

    // open index reader
    mm_idx_reader_t *r = mm_idx_reader_open(argv[1], &iopt, 0);
    mm_idx_t *mi;
    while ((mi = mm_idx_reader_read(r, n_threads)) != 0) { // traverse each part of the index
        mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
        mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
        gzrewind(f);
        kseq_rewind(ks);
        while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
            mm_reg1_t *reg;
            int j, i, n_reg;
            reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &mopt, 0); // get all hits for the query

            cerr << '\n';
            cerr << "id" << ' ' << reg->id << '\n';
            cerr << "n_reg" << ' ' << n_reg << '\n';
            cerr << "score" << ' ' << reg->score << '\n';
            cerr << "cnt" << ' ' << reg->cnt << '\n';
            cerr << "mlen" << ' ' << reg->mlen << '\n';
            cerr << "blen" << ' ' << reg->blen << '\n';

            for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
                mm_reg1_t *r2 = &reg[j];

                string type;
                if (r2->id == r2->parent) type = r2->inv? 'I' : 'P';
                else type = r2->inv? 'i' : 'S';

                assert(r2->p); // with MM_F_CIGAR, this should not be NULL
                printf("%s\t%d\t%d\t%d\t%c\t", ks->name.s, ks->seq.l, r2->qs, r2->qe, "+-"[r2->rev]);
                printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\tTP:%s\tcg:Z:", mi->seq[r2->rid].name, mi->seq[r2->rid].len, r2->rs, r2->re, r2->mlen, r2->blen, r2->mapq, type.c_str());

                for (i = 0; i < r2->p->n_cigar; ++i) // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
                    printf("%d%c", r2->p->cigar[i] >> 4, MM_CIGAR_STR[r2->p->cigar[i] & 0xf]);
                putchar('\n');
                free(r2->p);
            }
            free(reg);
        }
        mm_tbuf_destroy(tbuf);
        mm_idx_destroy(mi);
    }
    mm_idx_reader_close(r); // close the index reader
    kseq_destroy(ks); // close the query file
    gzclose(f);
    return 0;
}
