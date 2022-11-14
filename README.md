# GFAse

Tool for phasing genomic graph data using parental or proximity ligation data. 

**Recent Updates:**
  - Proximity linkage phasing working consistently on Shasta and Verkko graphs with low and high bubble N50s
  - PoreC tested and working

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
- zlib1g-dev
- libbz2-dev
- libcurl4-openssl-dev


## Usage

### Phase with proximity linkage reads (HiC or PoreC)

```
App description
Usage: ./phase_contacts_with_monte_carlo [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -i,--input TEXT REQUIRED    Path to BAM containing proximity linked reads. MUST be grouped by read name. Does not need index or regional sorting.
  -g,--gfa TEXT               Path to GFA containing assembly graph to be phased
  -o,--output_dir TEXT REQUIRED
                              Path to (nonexistent) directory where output will be stored
  -m,--min_mapq INT           (Default = 1)     Minimum required mapq value for mapping to be counted.
  -c,--core_iterations UINT   (Default = 200)   Number of iterations to use for each shallow convergence in the sampling process. The final phasing round uses 3*core_iterations.
  -s,--sample_size UINT       (Default = 30)    How many shallowly converged phase states to sample from. This is also the maximum usable concurrency (n_threads) for this stage of the pipeline.
  -r,--n_rounds UINT          (Default = 2)     How many rounds to sample and merge.
  -t,--threads UINT           (Default = 1)     Maximum number of threads to use.
  --use_homology              (Default = 0)     Use sequence homology to find alts. For whenever the GFA does not have Shasta node labels.
  --skip_unzip                (Default = 0)     After phasing nodes in the graph, DON'T unzip/concatenate haplotypes before writing to fasta. Unzipping should be skipped when using overlapped GFAs because no stitching is performed.
```

## Examples

### Prepare alignments (HiC)

```
bwa index Assembly-Phased.fasta \
&& \
bwa mem -t 46 -5 -S -P \
Assembly-Phased.fasta \
HG002.HiC_1_S1_R1_001.fastq \
HG002.HiC_1_S1_R2_001.fastq \
| samtools sort -n -@ 24 - -o hic_to_assembly.sorted_by_read.bam \
```
Once you have the BAM, pass it as an argument to GFAse, along with the GFA. Note that samtools and BWA are competing for threads.

![image](https://user-images.githubusercontent.com/28764332/201423346-b2077b90-7f96-42fc-be8a-5b655315ff3c.png)

x = Time (min)

### Prepare alignments (PoreC)

```
minimap2 \
-a \
-x map-ont \
-k 17 \
-t 56 \
-K 10g \
-I 8g \
Assembly-Phased.fasta \
porec_reads.fastq.gz \
| samtools view -bh -@ 8 -q 1 - \
> /extra2/data/human/hg002/porec/PAM_10290_and_10354_r10_porec_vs_shasta_r10_slow.q1.bam
```
PoreC does not used paired files. It is strongly recommended that you filter the mapq with a pipe to avoid massive output BAMs. No sorting is needed under the assumption that minimap2 groups output by read name. Once you have the BAM, pass it as an argument to GFAse, along with the GFA.

![image](https://user-images.githubusercontent.com/28764332/201423543-274923c9-6a2d-4f3c-93cd-dcd1e350e6e7.png)

x = Time (min)


### Phase Shasta with proximity data

```
/home/ubuntu/software/GFAse/build/phase_contacts_with_monte_carlo \
-i hic_to_assembly.sorted_by_read.bam \
-g Assembly-Phased.gfa \
-o /path/to/output/directory/ \
-m 1 \
-t 62
```

Run time depends almost entirely on the size of the BAM to be loaded. Filtering by map quality first (q>0) will save loading time. 

[Resource usage plots to be added]

### Phase Verkko (or any unnanotated, overlapped GFA) with proximity data

```
/home/ubuntu/software/GFAse/build/phase_contacts_with_monte_carlo \
-i hic_to_assembly.sorted_by_read.bam \
-g Assembly-Phased.gfa \
-o /path/to/output/directory/ \
--use_homology \
--skip_unzip \
-m 1 \
-t 62
```

Run time depends mostly on the size of the BAM to be loaded, however, additional run time and memory usage is incurred as a result of needing to rediscover alts/homologs in the GFA (with flag `--use_homology`). ~128GB and 48 threads should be sufficient. 

Note: Verkko should be homopolymer decompressed before running alignments and GFAse. Verkko produces overlapped GFAs which are not stitchable by GFAse, so `--skip_unzip` should be used.

![image](https://user-images.githubusercontent.com/28764332/201499328-75395981-d8fd-42d3-883e-b7c35be81e56.png)


### Phase Shasta with parental k-mers
A list of unique parental (separated into maternal/paternal fastas) kmers is required as input.
```
/home/ubuntu/software/GFAse/build/phase_haplotype_paths \
-i /home/ubuntu/data/human/hg002/assembly/Assembly-Detailed.gfa \
-p /home/ubuntu/data/human/kmer/hg03.48_55.unique.k31.fa \
-m /home/ubuntu/data/human/kmer/hg04.54_61.unique.k31.fa \
-k 31
```
![image](https://user-images.githubusercontent.com/28764332/169711827-7f84d3c6-51e8-465d-9620-f2da047a15a1.png)


## Example output

### Hi-C phasing
The plot below shows phase assignment of bubbles w.r.t. a diploid reference alignment. 
![image](https://user-images.githubusercontent.com/28764332/201426111-2941f038-9015-4abe-b649-b7cd59580051.png)

### K-mer based phasing

A recent ultra-long nanopore Shasta assembly with high phase-block contiguity after k-mer phasing (NG50 = 33.1Mbp)
![HG002_UL_Shasta_Phased_maternal_ideogram](https://user-images.githubusercontent.com/28764332/169709071-0d3696c2-8ffb-4cbd-b7af-4dd73ad83734.png)

Three examples of shifts in NGx in Shasta assemblies after applying GFAse with k-mers. Note that initial phased contiguity in Shasta is not necessarily a predictor of phased contiguity after GFAse.
![image](https://user-images.githubusercontent.com/28764332/169709283-db012bc4-5fc7-4eee-9901-59fe83293fd6.png)

