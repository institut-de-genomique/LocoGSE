# LocoGSE : Low coverage based Genome Size Estimator

Many pipelines or algorithms can estimate the size of a genome. owever, they are either based on kmers, a method that requires high sequencing depth (>20X) or mapping on an assembly, which is not always alvailable.   
Moreover, assembling and annotating genomes just to get some statistics about their theoric size is not very time-efficient.  
Here we propose a novel method based on mapping of reads on single copy genes , that is linearily related to the sequencing depth.  
This link can be expressed like that :


$$ DEPTH = \beta *  averagecoverage  * 1 000 000 $$

When the sequencing depth is known, the genome size can be estimated with:


$$ SIZE (in Mb) = \frac{totalofnucleotides}{\beta * averagecoverage * 1 000 000} $$


## Installation

### Conda environment 

It is recommended to use conda environment to use LocoGSE.

```bash
git clone https://github.com/institut-de-genomique/LocoGSE.git
cd LocoGSE/
conda env create -f environment.yml
```


## Arguments

|  Option  |  Parameter(s)  |  Description  |  Requirement  |
|---   |:-:   |:-:   |--:  |
|  `--reads`  |  `fastq1 ....gz`  |  Input `fastq` file |  Required if there is no `--multi_files`argument |
|  `--multi_files`  |  `multi_fastq.txt`  |  TXT file with, in each line, the path of files which are processed together (example in [wiki](https://github.com/institut-de-genomique/LocoGSE/wiki/4.LocoGSE-tutorial) ) |  Required if there is no `--reads`argument |
|  `--ref_prot`  |  `ref_prot_name`  | the prefix of FA file with all reference proteins and a .dmmd  file for [DIAMOND](https://github.com/bbuchfink/diamond)  |  Required. By Default : `OneKpGenes database`  | 
| `--recovery, -r`  | present or not |  Recovery option to continue a previous run after main steps  |  Optional  |
| `--threads, -t`  |  `int number`  |  Number of cpus for the pipeline  |  Optional  |
|  `--slope, -s `  |  `float number` |  Specific slope of a family.  |  Optional |
|  `--family, -f`  |  `name_of_family`  |  The name of family of the sample if it is in the database  |  Optional if any slope are given  |
|  `--list_families`  |  present or not  |  Print all families which are treated by LocoGSE   |  Optional  |
|  `--lineage`  |  `name_of_lineage`  |  The name of lineage of the sample if it is in the database  |  Optional if any slope or family are given  |
|  `--list_lineages`  |  present or not  |  Print all lineages which are treated by LocoGSE   |  Optional  |
| `--length_trim`  |  `int number`  |  Length (INT) of the sequence for the trimming  |  Optional (by default : 100) but not recommended if it is not for a new calibration |
|  `--trim`  |  present or not  |  Option to not trim fastq. Any character after `--trim` argument can work  |  Optional  | 
|  `--name_samples`  |  `TSV file`  |  A TSV with, in each line, the name of sample in the same line that in multi_files  |  Optional  | 
|  `--pegasus`  |  `yes`  |  Option to write a pegasus script to process many files in parallel. It's only for the mapping step and the user must launch pegasus  |  Optional  |
|  `--picog` |   `yes`  |  To convert the default unit default size (Mb) in picogram   |  Optional  |
|  `--lgprot, -l`  |  `TSV file`  |  A tsv file with each protein (first column) and their length (second column)  |  Optional  |
|  `--output, -o` |  `Name_of_output_directory`  |  A repository for the results  |  Optional. By default : results/   |
|  `--cleaning_output`  |  `yes`  |  A option to have only the deviant genes and the size of genome of your sample  |  Optional  |


## An example

`LocoGSE --multi_files merge_fastq.txt --ref_prot database --output output_dir/ -l name_prot_length.tsv -f Asteraceae --threads 4 --picog yes --cleaning_output yes `

With this command, LocoGSE maps each sample in merge_fasts.txt on the protein database with 4 cpus. After a filtering step, it searchs in the database the slope which is associated at Asteracea family. Moreover, it estimates the size of the genome sample in picogram and it cleans the output directory (output_dir) to keep only the list of deviant genes and the size.

## For more information, please read the [wiki](https://github.com/institut-de-genomique/LocoGSE/wiki/1.Home)


