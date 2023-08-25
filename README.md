# LocoGSE : Low coverage based Genome Size Estimator

Many pipelines or algorithms can estimate the size of a genome. However, they are either based on kmers, a method that requires high sequencing depth (>20X) or mapping on an assembly, which is not always available.   
Moreover, assembling and annotating genomes just to get some statistics about their theoric size is not very time-efficient.  
Here we propose a novel method (calibrated for most Angiosperm lineages) based on mapping of reads on single copy genes, that is linearily related to the sequencing depth.  

This link can be expressed like that :


$$ DEPTH = \beta *  average~Coverage  * 1 000 000 $$

When the sequencing depth is known, the genome size can be estimated with:


$$ SIZE (~in~Mb) = \frac{total~Nucleotides}{\beta * average~Coverage * 1 000 000} $$


## Installation

### Conda environment 

It is recommended to use a conda environment to use LocoGSE.

```bash
git clone https://github.com/institut-de-genomique/LocoGSE.git
cd LocoGSE/
conda env create -f environment.yml
conda activate LocoGSE
```


## Arguments
  - `--reads FASTQ_PATH`: Input fastq file. **Required if there is no `--list_fastq` argument**
  - `--list_fastq TXT_PATH`: Text file with, on each line, the list of fastq files to be treated together (same sample, see [wiki](https://github.com/institut-de-genomique/LocoGSE/wiki/4.LocoGSE-tutorial)). The first column can be the name of the sample. **Required if there is no `--reads` argument**
  - `--ref_prot DB_PREFIX`: Path to a monocopy protein database to be used with [DIAMOND](https://github.com/bbuchfink/diamond). The prefix of these 2 files must be given. **Required**. By Default: OneKP consensus obtained from https://github.com/smirarab/1kp/tree/master/alignments
  - `--slope NUMBER`: **Optional**. Slope (regression factor) used to estimate sequencing depth from depth on monocopy proteins. It is specific to each plant lineage. Pre computed slopes are available for families listed in `--list_families` and lineages in `--list_lineages`. There is no need to provide a slope if the species of interest is in the list, you can just provide either the family or the lineage.
  - `--threads NUMBER`: **Optional**. Number of CPUs to use durig the mapping step.
  - `--recovery`: **Optional**. If present, continue an interrupted run that was started in the output directory provided with `--output`
  - `--family NAME`: **Optional if `--slope` is present**.  Specify the plant family in order to use a pre-computed slope
  - `--list-families`: **Optional**. Print all families with available pre-computed slope
  - `--lineage NAME`: **Optional if `--slope` or `--family` is present**.  Specify the plant lineage in order to use a pre-computed slope
  - `--list-lineages`: **Optional**. Print all plant lineages with available pre-computed slopes
  - `--length_trim NUMBER`: **Optional**. Trim sequences to this size. Default: 100, highly recommended since the training step was performed with 100nt reads
  - `--no_trim`: **Optional**. Deactivates the trimming step.
  - `--pegasus yes|no`: **Optional**. Option to write a Pegasus script that can be manually launched to submit multiple mapping commands at the same time.
  - `--picog yes|no`: **Optional**. Converts default units (Mb) to picograms.
  - `--lgprot TSV_PATH`: **Optional**. A TSV file with protein name in the first column and their lengths in the second column. If not provided, LocoGSE will compute it.
  - `--cleaning_output yes|no`: **Optional**. Remove temporary files and only keep results (list of deviant genes, depth on monocopy gene set, and estimated genome size)

## An example

```
LocoGSE --list_fastq merge_fastq.txt \
  --ref_prot database --output output_dir/ \
  -l name_prot_length.tsv -f Asteraceae --threads 4 \
  --picog yes --cleaning_output yes
```

With this command, LocoGSE maps each sample in `merge_fasts.txt` on the protein database with 4 cpus. After a filtering step, it searches in the database the slope which is associated at Asteracea family. Moreover, it estimates the size of the genome sample in picogram and it cleans the output directory (`output_dir`) to keep only the list of deviant genes and the size.

## For information on how to calibrate the method, please read the [wiki](https://github.com/institut-de-genomique/LocoGSE/wiki/Home)


