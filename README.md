# LocoGSE : Low coverage based Genome Size Estimator

Many pipelines or algorithms can estimate the size of a genome. However, they are either based on kmers, a method that requires high sequencing depth (>20X) or mapping on an assembly, which is not always available.   
Moreover, assembling and annotating genomes just to get some statistics about their theoric size is not very time-efficient.  
Here we propose a novel method (calibrated for most Angiosperm lineages) based on mapping of reads on single copy genes, that is linearily related to the sequencing depth.  

This link can be expressed like that :


$$ DEPTH = \beta *  average~Coverage  * 1 000 000 $$

When the sequencing depth is known, the genome size can be estimated with:


$$ SIZE (~in~Mb) = \frac{total~Nucleotides}{\beta * average~Coverage * 1 000 000} $$


<p align="center">
  <img width="400" alt="schema_A" src="https://github.com/institut-de-genomique/LocoGSE/blob/main/images/schema_A.png?raw=true">
</p>

Link to the article: [click me](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2024.1328966/full)

## Installation

### Conda environment 

It is recommended to use a conda environment to use LocoGSE.

```bash
git clone https://github.com/institut-de-genomique/LocoGSE.git
cd LocoGSE/
conda env create -f environment.yml
conda activate LocoGSE
```

### Optional dependency

If [Fastoche](https://github.com/institut-de-genomique/fastoche) is accessible in PATH when running LocoGSE, it will be used to significantly speed up fastq parsing.

## Input file validation

LocoGSE requires one or more fastq files in order to do its computations. During the mapping step, Diamond only reports read names up to the first space in the sequence header. This can cause incorrect results being reported by LocoGSE and we highly encourage users to first check if some reads have identical names when trimmed after the first space. To do so, we include a command `check_seq_names` that can be run like this:
```bash
check_seq_names fastq1 fastq2 ...
```
This command will issue warnings if some reads need to be renamed. To proceed with renaming the duplicated identifiers, the [Seqkit](https://bioinf.shenwei.me/seqkit/) software suite can be used, like this:
```bash
seqkit rename reads.fastq > reads_renamed.fastq
```

## Arguments
- `--reads FASTQ_PATH`: Input fastq file. **Required if there is no --list_fastq argument**
- `--list_fastq TXT_PATH`: Text file with, on each line, the name of the sample followed by the list of fastq files to be treated together (same sample), separator=space. **Required if there is no --reads argument**. Caution: the same slope will be applied for all the samples. To launch genome size predictions in different lineages, LocoGSE needs to be run several times
- `--ref_prot DB_PREFIX`: **Optional**. Path to a monocopy protein database to be used with DIAMOND (the two files REF.dmnd and REF.fa must exist). By Default if not provided: OneKP consensus obtained from https://github.com/smirarab/1kp/tree/master/alignments/alignments-FAA-masked.tar.bz . Alternatively the option --busco can be used to run LocoGSE on Busco Embryophyta odb10 ancestral sequences
- `--busco`: **Optional**. Use Busco Embryophyta odb10 dataset and associated slopes instead of OneKP
- `--slope NUMBER`: **Optional**. Slope (regression factor) used to estimate sequencing depth from depth on monocopy proteins. It is specific to each plant lineage. Pre computed slopes are available for families listed in --list_families and lineages in --list_lineages. There is no need to provide a slope if the lineage corresponding to the species of interest is in the list : you can just provide either the family (--family) or the lineage (--lineage). If none is provided, default slope is 1.62 (global slope for all plants in the calibration dataset)
- `--slope-file TSV_PATH`: **Optional**. If one wants to use their own custom slopes. Path to a three-column TSV file with the header (#Family\tPhylo_group\tslope)
- `--view-slope-file`: **Optional** Print OneKP default slopes to standard output
- `--family NAME`: **Optional**. Specify the plant family in order to use a pre-computed slope
- `--list_families`: **Optional**. Print all families with available pre-computed slope
- `--lineage NAME`: **Optional**. Specify the plant lineage in order to use a pre-computed slope
- `--list_lineages`: **Optional**. Print all plant lineages with available pre-computed slopes
- `--no_trim`: **Optional**. Deactivates the trimming step.
- `--picog`: **Optional**. Converts default units (Mb) to picograms.
- `--lgprot TSV_PATH`: **Optional**. To be used if a custom db is provided with --ref_prot. A TSV file (sep=\t) with protein names in the first column and their lengths in the second column. If not provided, LocoGSE will compute it.
- `--cleaning_output`: **Optional**. If present, remove temporary files and only keep main results (depth on monocopy gene set, number of nucleotides in the readset(s), and estimated genome size)
- `--threads NUMBER`: **Optional**. Number of CPUs to use during the mapping step. Default=1
- `--pegasus`: **Optional**. If present, write a Pegasus script "pegasus_script.txt" that can be manually launched to submit multiple mapping commands at the same time. Also creates directories input_dir_pegasus and output_dir_pegasus
- `--recovery`: **Optional**. If present, continue an interrupted run that was started in the output directory provided with --output

## An example

```
LocoGSE --list_fastq merge_fastq.txt \
  --ref_prot database --output output_dir/ \
   -f Asteraceae --threads 4 \
  --picog yes --cleaning_output yes
```

With this command, LocoGSE maps each sample in `merge_fasts.txt` on the protein database with 4 cpus. After a filtering step, it searches in the database the slope which is associated with the Asteracea family and estimates the size of the genome sample in picograms. Then, the the output directory (`output_dir`) is cleaned  to keep only the calculated depth and estimated genome size.

## Calibration

The first calibration was performed on a set of 430 plants as described in [Guenzi-Tiberi et al](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2024.1328966/full).

The calibration was updated in december 2024 on a set of 1474 samples among those released by the [PhyloAlps project](https://www.ebi.ac.uk/ena/browser/view/PRJEB85061) for which a known 1Cx genome size was obtained from [ Kew plant C-values database](https://cvalues.science.kew.org/) in order to cover more plant lineages, especially gymnosperms and pteridophytes.

The calibration data are available at this [link](https://github.com/user-attachments/files/18560451/DATA_CALIBRATION_V2.csv).
 
In case one wants to use the calibration slopes described in [Guenzi-Tiberi et al.](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2024.1328966/full) , they can be downloaded [here](https://github.com/user-attachments/files/18560444/LocoGSE_PlantFamilies.CoeffRegression.V1.txt) and provided to LocoGSE with the option `--slope-file` .


## For information on how to calibrate the method for other lineages, please read the [wiki](https://github.com/institut-de-genomique/LocoGSE/wiki/Home)


