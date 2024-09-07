# DeConStruct pipeline

## Description
The DeConStruct pipeline constitutes a generalized framework for the identification of candidate RNA structural switches by ensemble deconvolution analysis of chemical probing data read out by mutational profiling (MaP), and the prioritization of functionally-relevant switches by evolutionary conservation assessment via covariation analysis.


## Requirements
The pipeline relies upon the following components:

1. [RNA Framework](https://github.com/dincarnato/RNAFramework) (v2.8 or greater)
2. [DRACO](https://github.com/dincarnato/draco) (v1.2 or greater)
3. [IncaRNAto lab tools](https://github.com/dincarnato/labtools)
4. [Infernal](https://github.com/EddyRivasLab/infernal) (v1.1.3 or greater)
5. [R-scape](https://github.com/EddyRivasLab/R-scape) (v2.0.0q or greater)

Make sure that the `lib/` folder of RNA Framework is added to the `$PERL5LIB` environment variable:

```bash
$ export PERL5LIB=$PERL5LIB:/path/to/RNA/Framework
```

## Running the pipeline
The DeConStruct pipeline involves a number of discrete steps to go from raw MaP data to prioritized RNA structure elements:

1. Read QC and preprocessing, alignment, mutation counting (using the RNA Framework)
2. Ensemble deconvolution (using DRACO)
3. Automated construction of structure-informed alignments and covariation analysis (using cm-builder, Infernal and R-scape)

The various steps and possible customizations are described in detail below. The sample command lines exemplify how to replicate the analyses performed in Borovska *et al.*, 2024.<br/>


### Read pre-processing &amp; alignment to reference
Aligment to the reference transcriptome can be performed using any tool. We typically use [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) through the `rf-map` module of the RNA Framework. First of all, build the reference transcriptome index:

```bash
# Obtain the E. coli transcriptome reference
$ wget https://www.incarnatolab.com/downloads/datasets/EcoliEnsembles_Borovska_2024/reference.tar.gz
$ tar -xzvf reference.tar.gz
$ cd reference

# Build the Bowtie2 index
# Replace <n> with the number of cores available on your machine
$ bowtie2-build --threads <n> Ecoli_TUs.fasta Ecoli_TUs
```

If working with paired-end reads, first clip sequencing adapters and merge pairs into long reads. For this purpose we typically use [Cutadapt](https://github.com/marcelm/cutadapt/) and [PEAR](https://cme.h-its.org/exelixis/web/software/pear/), but any other program (e.g., [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) and [COPEread](https://sourceforge.net/projects/coperead/)) can be used as well. <br/>
Note that, if working with single reads, the following steps can be omitted as `rf-map` will take care of clipping the sequencing adapters prior to read mapping:

```bash
# Clip the adapters
# Replace <n> with the number of cores available on your machine and specify the paths to the R1 and R2 fastq files
$ cutadapt -j <n> -A AGATCGGAAG -a AGATCGGAAG -m 100:100 -O 1 -o R1.fq -p R2.fq <path to fastq for R1> <path to fastq for R2>

# Merge pair
# Replace <n> with the number of cores available on your machine and <ram> with the amount of RAM to be used (e.g., 10G to use 10 gigabytes)
$ pear -f R1.fq -r R2.fq -o merged_reads -n 100 -q 20 -u 0 -e -y <ram> -j <n> -z
$ cat merged_reads.discarded.fastq >> merged_reads.assembled.fastq
```

In the above example, the two output files from PEAR are merged. The `merged_reads.discarded.fastq` file contains reads from pairs that could not be assembled, possibly because the sequenced fragments are longer than the sum of the lengths of the two mates. Nonetheless, these reads are long enough to be used for the downstream analysis, therefore they can be retained. Please note that PEAR will automatically reverse complement the R2 reads before writing them to this file.<br/>

Reads can then be aligned to the reference transcriptome:

```bash
# Replace <f> with the number of fastq files to be processed in parallel (in this example 1) and <n> with the number of cores available on your machine
# Note: at least <p> * <n> cores must be available on the machine
$ rf-map -b2 -cq5 20 -ctn -cmn 0 -cl 90 -mp "--very-sensitive-local --nofw" -b5 6 -p <f> -wt <n> merged_reads.assembled.fastq
```

This will generate the `rf_map/` output folder, containing the aligned reads in BAM format.<br/>

The above command line assumes that libraries have been generated using a protocol that preserves the strandedness of the sequenced reads, particularly a *first-strand* protocol, meaning that reads are expected to align to the reverse complement of the transcript sequences. Therefore, the parameter `--nofw` tells to Bowtie2 to ignore the forward strand. In case of libraries generated using a *second-strand* protocol, meaning that reads are expected to align to the actual transcript sequences, the parameter `--norc` must be specified instead, while both parameters have to be omitted for *unstranded* protocols in which reads have a nearly equal probability of aligning to the transcript sequences, or their reverse complement.<br/>
The parameter `-b5 6` further removes the 6 5'-most bases of each read. Assuming that reverse transcription has been performed using random hexamers, these bases are discarded to remove potential mispriming artifacts. If using a  *second-strand* protocol instead, this has to be replaced with `-b3 6`. If reverse transcription is non-random primed (i.e., adapters are directly ligated to the RNA, and a primer complementary to the 3' adapter is used to drive reverse transcription), both parameters can be omitted.<br/>


### Mutation counting &amp; clean-up
This step involves processing the alignments into vectors of mutations that will be used by DRACO to perform ensemble deconvolution, and it is performed through the `rf-count` module of the RNA Framework:

```bash
# Replace <f> with the number of BAM files to be processed in parallel (in this example 1) and <n> with the number of cores available on your machine
# Note: at least <p> * <n> cores must be available on the machine
$ rf-count -p <f> -wt <n> -m -mm -wl 2000 -ds 100 -es -na -ni -md 1 -dc 3 -me 0.1 rf_map/*.bam
```

A number of filtering steps are applied to only retain high-quality reads/mutations. Please note that these steps have been optimized for [DMS-MaPseq experiments](https://pubmed.ncbi.nlm.nih.gov/27819661/) read out using group II intron reverse transcriptases (e.g., TGIRT-III), and might need to be carefully optimized for SHAPE-MaP experiments: 

- The `-ds 100` parameter causes reads covering &lt; 100 nt of the reference transcript to be discarded. This allows discarding low-quality alignments in which several bases have been clipped by Bowtie2.
- The `-es` parameter causes mutations residing within low-quality segments of the read (so, mutations residing between bases whose quality score is &lt; 20) to be discarded.
- The `-na` parameter causes ambiguously aligned deletions (i.e., deletions residing within repetitive regions, for which it's impossible to determine which of the repeats were actually deleted) to be discarded.
- The `-ni` parameter causes insertions to be ignored (as these account for &lt; 1% of mutations in DMS-MaPseq experiments).
- The `-md 1` parameter only retains deletions of up to 1 nt.
- The `-dc 3` causes mutations residing within 3 nt from each other to be discarded.

The above comand will generate the `rf_count/` output folder containing the Mutation Map (MM) files to be used for ensemble deconvolution. MM files will contain a number of non-informative reads, as well as regions with overly high sequencing coverage, which need to be removed prior to analysis with DRACO:

```bash
$ filterMM -ac -mpr 2 -xc 500000 -o dataset.mm rf_count/merged_reads.assembled.mm
```

In the above command:

- The `-ac` parameter causes only A/C mutations to be retained
- The `-mpr 2` parameter causes reads harboring &lt; 2 mutations to be discarded
- The `-xc 500000` causes regions exceeding a coverage of 500,000X to be randomly downsampled to achieve a maximum coverage of 500,000X<br/>


### Ensemble deconvolution
Ensemble deconvolution involves identifying regions of the transcriptome populating multiple conformations and reconstructing their individual reactivity profiles, and it is perfomed using `draco`:

```bash
# Replace <n> with the number of cores available on your machine
$ draco --absWinLen 100 --absWinOffset 5 --minBaseCoverage 2000 --minFilteredReads 2000 --minPermutations 10 --maxPermutations 50 --firstEigengapShift 0.95 --lookaheadEigengaps 1 --softClusteringIters 30 --softClusteringInits 500 --softClusteringWeightModule 0.005 --mm dataset.mm --output dataset.json --processors <n> --whitelist rf_count/whitelists/merged_reads.assembled.txt
```

In the above command line:

- The `--absWinLen 100` and `--absWinOffset 5` parameters specify the size (in nt) of the sliding window for the ensemble deconvolution analysis, and the sliding offset. Decreasing the offset will increase the sensitivity of the analysis, but it will significantly increse runtimes.
- The `--minBaseCoverage 2000` and `--minFilteredReads 2000` parameters cause bases or windows with coverage &lt; 2,000X not to be used for ensemble deconvolution. It is important to keep in mind that, as in the previous step we filtered out non-informative reads from the MM file, this coverage is only calculated on informative reads, meaning that the effective coverage of these regions is significantly higher.
- The `--minPermutations 10` and `--maxPermutations 50` parameters respectively define the mimum and maximum number of permutations to be performed to generate the null model that will be used to determine the number of coexisting conformations in the ensemble.
- The `--softClusteringIters 30` and `--softClusteringInits 500` parameters control the number of soft clustering iterations to be performed to cluster the reads across the identified conformations. Higher values will increase the accuracy of the analysis, further ensuring that the optimal clustering is achieved, but will also significantly impact runtimes. The same principle applies to decreasing values of the `--softClusteringWeightModule` parameter.

DRACO will generate the output JSON file `dataset.json`. This file needs to be processed using the `rf-json2rc` module of the RNA Framework to produce RNA Count (RC) files:

```bash
$ rf-json2rc -j dataset.json -r rf_count/merged_reads.assembled.rc -ec 1000 -mom 0.35 -ow -e 20 -cf 0.1 -i 0.1 -mcm 0.65 -mcr 0.65 -ki
```

where `dataset.json` and `rf_count/merged_reads.assembled.rc` are respectively the JSON output generated by DRACO and the RC output generated by `rf-count`.

If working with multiple datasets or replicates (e.g., DH5&#593; and TOP10 *E. coli* as in Borovska *et al.*, 2024), multiple files can be specified as a comma-separated list:

```bash
$ rf-json2rc -j DH5a.json,TOP10.json -r DH5a.rc,TOP10.rc -ec 1000 -mom 0.35 -e 20 -cf 0.1 -i 0.1 -mcm 0.65 -mcr 0.65 -ki
```

In the above command line:

- The `-ec 1000` parameter causes windows with a median cumulative coverage (the sum of the coverage across all the conformations for that window) &lt; 1,000X to be discarded
- The `-mom 0.35` parameter causes windows overlapping by at least 35% of their length to be merged (intra-experiment). Merging requires overlapping windows to populate the same number of conformations and to exceed a certain correlation threshold (specified by the `-mcm 0.65` parameter, corresponding to a Pearson correlation coefficient &ge; 0.65)
- The `-mcr 0.65` parameter controls the mimimum Pearson correlation coefficient between matching windows across different experiments/replicates
- The `-cf 0.1` parameter sets value to cap raw reactivities prior to correlation calculation 
- The `-i 0.1` parameter specifies that the first and last 10% of the bases in each window must be excluded from correlation calculation
- The `-ki` paramter causes reactivity values for bases ignored through the `-i` parameter to be reported in the output file (by default they would be set to NaN)
- The `-e 20` parameter causes the identified windows to be enlarged by &plusmn; 20 nt

The output `rf_json2rc/` folder will contain:

- The `stoichiometries.txt` file, containing the list of regions populating &ge; 2 conformations that were matched between the various experiments/replicates, along with the individual stoichiometries of the different conformations across each experiment/replicate
- An RC file for each analyzed experiment/replicate, containing per-base mutation counts and coverage of the deconvolved conformations for the windows that were matched<br/>


### Reactivity normalization &amp; structure modeling
RC files have to be normalized to obtain reactivities that can be used to constrain structure predictions. This step is performed using the `rf-norm` module of the RNA Framework:

```bash
# Replace <n> with the number of cores available on your machine
$ cd rf_json2rc/
$ for f in *.rc; do rf-norm -sm 4 -nm 3 -rb AC -mm 1 -n 100 -t $f -p <n>; done
```

In the above command line:

- The `-sm 4` and `-nm 3` parameters respectively specify the reactivity calculation and normalization methods to be used. For additional details please refer to the [RNA Framework docs](https://rnaframework-docs.readthedocs.io/en/latest/rf-norm/)
- The `-rb AC` parameter defines the reactive bases (A and C, as in this case, for DMS-MaPseq experiments)
- The `-n 100` parameter defines the minimum coverage of a base to be retained (bases with coverage &lt; 100X are reported as NaNs)

The above command will generate a folder of normalized reactivity XML files for each experiment/replicate, characterized by the suffix `_norm`. Each XML file will correspond to a different window/conformation. Files are named after the respective transcript and start-end coordinates. Files corresponding to alternative conformations for the same window share the same name, but are characterized by a different suffix (e.g., `_c0`, `_c1`, etc.).

Structure modeling can then be performed using the `consensusFold` tool. This tool can simultaneously use reactivity data from multiple experiments/replicates to yield a consensus secondary structure prediction:

```bash
# Replace <n> with the number of cores available on your machine
$ consensusFold -sl 4.8 -in -0.8 -md 600 -p <n> -g *_norm/
```

The `-md 600` parameter defines the maximum allowed base-pairing distance, while the `-sl 4.8` and `-in -0.8` respectively define the *slope* and *intercept* folding parameters. These parameters are highly reagent- and protocol-specific, and ideally they should be optimized for every experiment. We found that 4.8 and -0.8 tend to be optimal parameters across a variety of DMS-MaPseq experiments.

The output `consensusFold_out/` folder will contain two sub-directories:

- `structures/`, containing the modeled secondary structures for each conformation, in dot-bracket format
- `images/`, containg SVG plots of DMS reactivities and base-pairing probabilities<br/>


### Evolutionary conservation assessment by covariation analysis
This final step allows prioritizing regions of the transcriptome that are under strong purifying selection. The presence of significant covariation is usually a strong evidence for functionality of RNA structures, therefore this step can help mining functional RNA structural switches out of a plethora of structurally-heterogeneous regions identified by DRACO.

This analysis step relies on a modified version of `cm-builder`, which is capable of analyzing entire bacterial genomes. Before proceeding with the covariation analysis, all structures modeled in the previous step need to be consolidated into a single dot-bracket file. It is important to point out that, in some cases, two DRACO-deconvolved reactivity profiles for a same transcriptome region, corresponding to two alternative conformations making up the ensemble for that transcript, might converge on the same (or a very similar) structure during modeling. This can happen for different reasons, such as limits in the thermodynamic model, or conformations sharing a very similar secondary structure but differing at the tertiary structure level.

Therefore, prior to performing the covariation analysis, alternative conformations showing very similar secondary structures are collapsed into a single representative structure (to reduce runtimes) using the `removeRedundantStructs` script included in this repository:

```bash
# Consolidate a non-redundant set of structures
$ cat consensusFold_out/structures/*.db > all_structs.db
$ removeRedundantStructs -f 0.9 -o non_redundant.db all_structs.db
```

The `-f 0.9` parameter defines the Fowlkes-Mallows index (FMI) threshold to define two structures as being *similar*, hence to be collapsed. The FMI is essentially the geometric mean of positive predictive value (PPV) and sensitivity, and it is a measure of the fraction of base-pairs shared between two structures. It ranges between 0 (totally dissimilar structures) and 1 (identical structures).<br/>
By default, a random one between the two structures is selected. If the `-c` parameter is specified, the consensus structure (so the structure having only the base-pairs shared by both structures) is reported instead.

After having generated a non-redundant set of strucures, it is possible to proceed with to the covariation analysis:

```bash
# Replace:
# <path to FASTA file of bacterial genomes> with the path to a multi-FASTA file containing the genomic sequences of bacteria to be searched for structural homologs
# <path to reference organism genome> with the path to a FASTA file containg the sequence of the genome of the organism the structures have been extracted from (in this case, E. coli U00096.3, available from: https://www.ncbi.nlm.nih.gov/nuccore/545778205)
# <n> with the number of cores available on your machine
$ cm-builder -m non_redundant.db -d <path to FASTA file of bacterial genomes> -s <path to reference organism genome> -q 0.75 -t 1 -k -c <n>
-o cmbuilder/ -M cmbuilder/tmp/ -g -B -T 20 -I 10 -N -S
```

In the above command line:

- `-g` enables the genome mode, hence allowing the search for putative homologs on both genomic strands. It is possible to omit this parameter if a set of transcript sequences is passed via the `-d` parameter, instead of genome sequences
- `-B` causes the program to use bit-scores instead of E-values to select candidate homologs
- `-T 20` sets a bit-score inclusion threshold of 20
- `-I 10` causes the bit-score inclusion threshold to be increased by 10 at every iteration of the algorithm to refine the selected homologs
- `-N` allows estimating a noise-threshold for the bit-score. This is done by taking 10% of the database of bacterial genomes, reversing it, and using it to look for putative homologs. The match with the highest bit-score is used to define the noise threshold. If the noise threshold exceeds the value defined by `-T` (and increased by `-I` at every iteration), the bit-score threshold is then set to the noise threshold
- `-q 0.75` sets the minimum query coverage, so that truncated matches covering &lt; 75% of the query structure are discarded
- `-t 1` defines the positional tolerance. This value must __ALWAYS__ be set to 1 when working with full genomes as the relative position of the homologous structure might vary from genome to genome
- `-k` retains the constructed Stockholm alignments

The database of bacterial genomes to be scanned for homologs can be obtained from NCBI (see [here](https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/) for example).
Depending on the size of the provided database of bacterial genomes, and the number of structures to be evaluated, this process might require up to several days. At the end of the analysis, the `cmbuilder/` output directory will contain a number of constructed alignments in Stockholm format.<br/>
Alignments can optionally be polished using the `stockholmPolish` tool:

```bash
# Replace <n> with the number of cores available on your machine
$ stockholmPolish -p <n> cmbuilder/
```

The polished alignments can then be prioritized using the `selectAln` script included in this repository:

```bash
# Replace <n> with the number of cores available on your machine
$ selectAln -p <n> -he 0.05 -be 0.1 -hn 1 -bn 5 -hf 0.25 -bf 0.125 stockholmPolish_out/
```

In the above command line:

- `-he 0.05` defines the E-value threshold to define a helix as significantly covarying
- `-be 0.1` defines the E-value threshold to define a base-pair as significantly covarying
- `-hn 1` defines the minimum number of helices that need to covary to keep an alignment
- `-bn 5` defines the minimum number of base-pairs that need to covary to keep an alignment
- `-hf 0.25` defines the minimum fraction of helices that need to covary (according to the `-he` threshold) to keep an alignment
- `-bf 0.125` defines the minimum fraction of base-pairs that need to covary (according to the `-be` threshold) to keep an alignment

The `selected_aln/` output directory will contain:

- `stats.out`, with a list of all alignments processed, the calculated values and whether the alignment was selected or not
- `alignments/`, a folder containing all the selected alignments

The selected structures/alignments can then be pursued for further functional characterization.
