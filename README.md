# GFAse

Tool for phasing genomic graph data using parental or proximity ligation data. Currently specialized for Shasta assembly graphs.

## Installation

Tested on Ubuntu20.04 

```
git clone https://github.com/rlorigro/GFAse.git
cd GFAse
mkdir build
cd build
cmake ..
make -j [n_threads]
```

### Dependencies

- cmake 3.0+
- c++17
- gcc
- OpenMP
- zlib


## Example usage

### Phase with Hi-C
A BAM file containing alignments of Hi-C data to the assembly contigs is required as input.
```
/home/ubuntu/software/GFAse/build/phase_contacts \
-i /extra/data/human/hg002/align/hg002_hic_s1-2_to_shasta_guppy6_run14.q1.bam \
-g /home/ubuntu/data/human/hg002/assembly/Assembly-Phased.gfa \
-o /home/ubuntu/data/human/hg002/assembly/hic_phase_run2/ \
-m 2 \
-p PR \
-t 46
```

![image](https://user-images.githubusercontent.com/28764332/169709430-4bac6bbd-441a-4cf6-91cb-a54041679a81.png)


### Phase with parental k-mers
A list of unique parental (separated into maternal/paternal fastas) kmers is required as input.
```
/home/ubuntu/software/GFAse/build/phase_haplotype_paths \
-i /home/ubuntu/data/human/hg002/assembly/Assembly-Detailed.gfa \
-p /home/ubuntu/data/human/kmer/hg03.48_55.unique.k31.fa \
-m /home/ubuntu/data/human/kmer/hg04.54_61.unique.k31.fa \
-k 31
```

## Example output

### Hi-C phasing
The plot below shows phase assignment of bubbles w.r.t. a diploid reference alignment. 

- Blue/Green = 01/10 (trans)
- Red/Orange = 00/11 (cis)

![image](https://user-images.githubusercontent.com/28764332/169707905-7d6e688f-e07a-4cdc-a9e7-e19893132d80.png)

### K-mer based phasing

A recent ultra-long nanopore Shasta assembly with high phase-block contiguity after k-mer phasing (NG50 = 33.1Mbp)
![HG002_UL_Shasta_Phased_maternal_ideogram](https://user-images.githubusercontent.com/28764332/169709071-0d3696c2-8ffb-4cbd-b7af-4dd73ad83734.png)

Three examples of shifts in NGx in Shasta assemblies after applying GFAse with k-mers. Note that initial contiguity is not necessarily a predictor of phased contiguity.
![image](https://user-images.githubusercontent.com/28764332/169709283-db012bc4-5fc7-4eee-9901-59fe83293fd6.png)

