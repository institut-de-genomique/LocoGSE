# LocoGSE : Low coverage based Genome Size Estimator

Many pipelines or algorithms can estimate the size of a genome. However, they are either based on kmers, a method that requires high sequencing depth (>20X) or mapping on an assembly, which is not always available.   
Moreover, assembling and annotating genomes just to get some statistics about their theoric size is not very time-efficient.  
Here we propose a novel method (calibrated for most Angiosperm lineages) based on mapping of reads on single copy genes, that is linearily related to the sequencing depth.  

This link can be expressed like that :


$$ DEPTH = \beta *  average~coverage  * 1 000 000 $$

When the sequencing depth is known, the genome size can be estimated with:


$$ SIZE (~in~Mb) = \frac{total~of~nucleotides}{\beta * average~coverage * 1 000 000} $$


## Installation

### Conda environment 

It is recommended to use conda environment to use LocoGSE.

```bash
git clone https://github.com/institut-de-genomique/LocoGSE.git
cd LocoGSE/
conda env create -f environment.yml
conda activate LocoGSE
```


## Arguments

|  Option  |  Parameter(s)  |  Description  |  Requirement  |
|---   |:-:   |:-:   |--:  |
|  `--reads`  |  `fastq1 ....gz`  |  Input `fastq` file |  Required if there is no `--list_fastq`argument |
|  `--list_fastq`  |  `list_fastq.txt`  |  txt file with, on each line, the list of fastq files to be treated together (same sample, see [wiki](https://github.com/institut-de-genomique/LocoGSE/wiki/4.LocoGSE-tutorial) . The first column can be the name of the sample. |  Required if there is no `--reads`argument |
|  `--ref_prot`  |  `ref_prot_name`  | Path to a monocopy protein database to be used with [DIAMOND](https://github.com/bbuchfink/diamond). The prefix of these 2 files must be given.  |  Required. By Default : OneKP consensus obtained from https://github.com/smirarab/1kp/tree/master/alignments  | 
| `--recovery, -r`  | present or not |  Recovery option to continue the run started in the output directory provided   |  Optional  |
| `--threads, -t`  |  `int number`  |  Number of CPUs to be used during the mapping step  |  Optional  |
|  `--slope, -s `  |  `float number` |  Slope (regression factor) used to estimate sequencing depth from depth on monocopy proteins. It is specific to each plant lineage. Pre computed slopes are available for families listed in --list_families and lineages in --list_lineages : no need to provide a slope if your species is in the list, you can just provide either the family or the lineage.  |  Optional |
|  `--family, -f`  |  `name_of_family`  |  Specify the plant family in order to use a pre-computed slope  |  Optional if any slope are given  |
|  `--list_families`  |  present or not  |  Print all families with available pre-computed slope   |  Optional  |
|  `--lineage`  |  `name_of_lineage`  |  Specify the plant lineage in order to use a pre-computed slope  |  Optional if any slope or family are given  |
|  `--list_lineages`  |  present or not  |  Print all plant lineages with available pre-computed slope   |  Optional  |
| `--length_trim`  |  `int number`  |  Length (INT) to trim the sequences to  |  Optional (by default : 100 , highly recommended since the training step was performed with 100nt reads !) |
|  `--no_trim`  |  present or not  |  Desactivates the trimming step  |  Optional  | 
|  `--pegasus`  |  `yes`  |  Option to write a Pegasus script that can be manually launched to submit multiple mapping commands at the same time  |  Optional  |
|  `--picog` |   `yes`  |  To convert the default unit default size (Mb) to picograms   |  Optional  |
|  `--lgprot, -l`  |  `TSV file`  |  A tsv file with each protein (first column) and their length (second column)  |  Optional : if not provided, they will be computed  |
|  `--output, -o` |  `Name_of_output_directory`  |  A repository for the results  |  Optional. By default : results/   |
|  `--cleaning_output`  |  `yes`  |  Remove temporary files and only keep results : list of deviant genes, depth on monocopy gene set, and estimated genome size  |  Optional  |


## An example

`LocoGSE --list_fastq merge_fastq.txt --ref_prot database --output output_dir/ -l name_prot_length.tsv -f Asteraceae --threads 4 --picog yes --cleaning_output yes `

With this command, LocoGSE maps each sample in merge_fasts.txt on the protein database with 4 cpus. After a filtering step, it searches in the database the slope which is associated at Asteracea family. Moreover, it estimates the size of the genome sample in picogram and it cleans the output directory (output_dir) to keep only the list of deviant genes and the size.

## For information to calibrate the method, please read the [wiki](https://github.com/institut-de-genomique/LocoGSE/wiki/1.Home)


