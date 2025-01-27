#!/usr/bin/env python3

from .lib import prediction
from .lib import mapping
from .lib import filtering
from .lib import checking
from .lib import trimming
from .pegasus import writing_script
import argparse
import os
import sys
import time
import shutil


def run():
    parser = argparse.ArgumentParser(
        prog="LocoGSE : Low coverage Genome Size Estimator",
        description="\n\n A Genome Size Estimation program. It is based on a linear relation between the sequencing depth (linked to the genome size) and the depth of mapping short reads on a set of monocopy genes. \n The regression factor (slope) depends of the plant family/lineage : slopes are precomputed for a number of plant families. \n For questions : https://github.com/institut-de-genomique/LocoGSE/issues or pierre.guenzi.tiberi@gmail.com or fdenoeud@genoscope.cns.fr",
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=True,
    )

    args_grp = parser.add_argument_group("Arguments")

    # Mandatory arguments

    # Fastq file (short reads)
    args_grp.add_argument(
        "--reads",
        action="store",
        nargs="*",
        dest="reads",
        help="Input fastq file. Required if there is no --list_fastq argument",
        default="",
        required=False,
    )

    # Multiple files : listfastq = list of absolute paths to several fastq files (if you put several files on one line they will be treated together: for instance read1.fastq and read2.fastq)
    args_grp.add_argument(
        "--list_fastq",
        action="store",
        dest="list_fastq",
        help="Text file with, on each line, the name of the sample followed by the list of fastq files to be treated together (same sample), separator=space. Required if there is no --reads argument. Caution: the same slope will be applied for all the samples. To launch genome size predictions in different lineages, LocoGSE needs to be run several times",
        default="",
        required=False,
    )

    # FASTA file with all reference proteins
    args_grp.add_argument(
        "--ref_prot",
        action="store",
        dest="ref",
        help="Optional. Path to a monocopy protein database to be used with DIAMOND (the two files REF.dmnd and REF.fa must exist). By Default if not provided: OneKP consensus obtained from https://github.com/smirarab/1kp/tree/master/alignments/alignments-FAA-masked.tar.bz . Alternatively the option --busco can be used to run LocoGSE on Busco Embryophyta odb10 ancestral sequences",
        default="",
        required=False,
    )

    # Use Busco embryophyta proteins and slopes
    args_grp.add_argument(
        "--busco",
        action="store_true",
        dest="use_busco",
        help="Optional. Use Busco Embryophyta odb10 dataset and associated slopes instead of OneKP",
        required=False,
    )

    # Slope
    args_grp.add_argument(
        "--slope",
        "-s",
        action="store",
        dest="slope",
        help="Optional. Slope (regression factor) used to estimate sequencing depth from depth on monocopy proteins. It is specific to each plant lineage. Pre computed slopes are available for families listed in --list_families and lineages in --list_lineages. There is no need to provide a slope if the lineage corresponding to the species of interest is in the list : you can just provide either the family (--family) or the lineage (--lineage). If none is provided, default slope is 1.62 (global slope for all plants in the calibration dataset)",
        default=1.62,
        required=False,
    )
    args_grp.add_argument(
        "--slope-file",
        action="store",
        dest="slope_file",
        help="Optional. If one wants to use their own custom slopes. Path to a three-column TSV file with the header (#Family\\tPhylo_group\\tslope)",
        default=None,
        type=os.path.abspath,
        required=False,
    )
    args_grp.add_argument(
        "--view-slope-file",
        action="store_true",
        dest="view_slope_file",
        help="Prints the default OneKP slopes to standard output",
        default=False,
        required=False,
    )

    # Family
    args_grp.add_argument(
        "--family",
        "-f",
        action="store",
        dest="family",
        help="Optional. Specify the plant family in order to use a pre-computed slope",
        default="",
        required=False,
    )

    # List of families
    args_grp.add_argument(
        "--list_families",
        action="store_true",
        dest="list_families",
        help="Optional. Print all families with available pre-computed slope",
        default=False,
        required=False,
    )

    # Lineage
    args_grp.add_argument(
        "--lineage",
        action="store",
        dest="lineage",
        help="Optional. Specify the plant lineage in order to use a pre-computed slope",
        default="",
        required=False,
    )

    # List of lineages
    args_grp.add_argument(
        "--list_lineages",
        action="store_true",
        dest="list_lineages",
        help="Optional. Print all plant lineages with available pre-computed slopes",
        default=False,
        required=False,
    )

    # Trimming
    args_grp.add_argument(
        "--no_trim",
        action="store_true",
        dest="no_trim",
        help="Optional. Deactivates the trimming step",
        default=False,
        required=False,
    )

    # Size in pg
    args_grp.add_argument(
        "--picog",
        action="store_true",
        dest="picog",
        help="Optional. Converts default units (Mb) to picograms",
        default=False,
        required=False,
    )

    # Size of each reference protein, if the user provides another protein set (not OneKp)
    args_grp.add_argument(
        "--lgprot",
        "-l",
        action="store",
        dest="lgprot",
        help="Optional. To be used if a custom db is provided with --ref_prot. A TSV file (sep=\\t) with protein names in the first column and their lengths in the second column. If not provided, LocoGSE will compute it",
        default="",
        required=False,
    )

    # Clean output with only sample size
    args_grp.add_argument(
        "--cleaning_output",
        action="store_true",
        dest="clean",
        help="Optional. If present, remove temporary files and only keep main results (depth on monocopy gene set, number of nucleotides in the readset(s), and estimated genome size)",
        default=False,
        required=False,
    )

    # Threads
    args_grp.add_argument(
        "--threads",
        "-t",
        action="store",
        dest="threads",
        help="Optional. Number of CPUs to use during the mapping step. Default=1",
        default=1,
        required=False,
    )

    # Pegasus option
    args_grp.add_argument(
        "--pegasus",
        action="store_true",
        dest="pegasus",
        help="Optional. If present, write a Pegasus script 'pegasus_script.txt' that can be manually launched to submit multiple mapping commands at the same time. Also creates directories input_dir_pegasus and output_dir_pegasus",
        default=False,
        required=False,
    )

    # Recovery
    args_grp.add_argument(
        "--recovery",
        "-r",
        action="store_true",
        dest="recovery",
        help="Optional. If present, continue an interrupted run that was started in the output directory provided with --output",
        default=False,
        required=False,
    )

    # Sequence length
    args_grp.add_argument(
        "--length_trim",
        action="store",
        dest="length_sequence",
        help="Optional. Trim sequences to this size. Default: 100, highly recommended since the training step was performed with 100nt reads",
        default="100",
        required=False,
    )

    # Output
    args_grp.add_argument(
        "--output",
        "-o",
        action="store",
        dest="output_dir",
        help="Output directory name",
        default="results",
        required=False,
    )

    args = parser.parse_args()
    slope = args.slope

    if args.view_slope_file:
        path_module = os.path.abspath(__file__)
        path_database = path_module.replace(
            "LocoGSE.py", "slopes/PlantFamilies.CoeffRegression.V2.txt"
        )
        with open(path_database) as inf:
            for line in inf:
                print(line, end="")
            sys.exit(-1)

    if args.list_families:
        checking.list_families_print(args.use_busco)
        sys.exit(-1)

    if args.list_lineages:
        checking.list_lineages_print(args.use_busco)
        sys.exit(-1)

    # Loading samples
    name_reads = []
    number_nt_list = []
    list_fastq = []
    if args.reads == "":
        if args.list_fastq != "":
            number_nt_list, list_fastq, name_samples = checking.complete_multi_files(
                args.list_fastq
            )
        else:
            print("Please provide samples with either --reads or --list_fastq.")
            sys.exit(-1)
    else:
        number_nt_list, list_fastq, name_samples = checking.complete_single_files(
            args.reads
        )

    # Checking
    if args.ref == "":
        if args.use_busco:
            print(
                "No argument given to --ref_prot, defaulting to the Busco embryophyta database"
            )
            path_main = os.path.abspath(__file__)
            args.ref = path_main.replace("LocoGSE.py", "db/BUSCO.ancestral")
        else:
            print("No argument given to --ref_prot, defaulting to the OneKP database")
            path_main = os.path.abspath(__file__)
            args.ref = path_main.replace("LocoGSE.py", "db/OneKP.410genes.consensus")
    else:
        args.ref = os.path.abspath(args.ref)

    if not args.recovery:
        try:
            os.mkdir(args.output_dir)
        except:
            print(
                f"\n Output directory {args.output_dir} can not be created, please erase it before launching the programm !"
            )
            sys.exit(1)

    if args.lgprot != "":
        lgprot_path = os.path.abspath(args.lgprot)

    print(
        "Warning: when using multiple samples, they must belong to the same lineage as the slope used will be the same for every sample"
    )

    if args.family != "":
        slope = prediction.determine_slope_for_family(
            args.family, args.use_busco, args.slope_file
        )
    elif args.lineage != "":
        slope = prediction.determine_slope_for_lineage(
            args.lineage,
            args.use_busco,
            args.slope_file,
        )

    actual_path = os.getcwd()
    os.chdir(args.output_dir)
    global_start = time.perf_counter()

    # Checking steps
    mapping_step_finished = checking.checking_step_recovery_mapping()
    counting_nucleotides = checking.checking_step_recovery_computing_nbnt()
    checking.checking_step(args.ref, slope, mapping_step_finished)
    slope = float(slope)

    # Write pegasus script
    if args.pegasus:
        writing_script(
            list_fastq, args.threads, args.ref, name_samples, slope, actual_path
        )
        sys.exit(-1)

    # Trimming step
    if not args.no_trim:
        list_fastq = trimming.trimming_step(list_fastq, args.threads)
    else:
        # No trimming
        print("No trimming")

    # Map reads on protein
    diamond_df_list = []
    if mapping_step_finished != True:
        diamond_df_list = mapping.mapping_multi_files(
            list_fastq, args.ref, args.threads, name_samples
        )
    else:
        diamond_df_list = mapping.recovery_mapping_step(
            list_fastq, args.ref, args.threads, name_samples
        )

    # Find one best hit per read in diamond DataFrame
    besthit = []
    besthit = filtering.searching_besthit(diamond_df_list)

    # Compute the length of each protein

    if not args.lgprot:
        filtering.computing_length_prot_database(args.ref)
        lgprot_path = "ref_prot_length/ref_prot.tsv"

    # Find outlier genes and create a file after removing the outliers

    filtering.filter_sample(lgprot_path, besthit)

    # Compute the number of nucleotides in each FASTQ File (or list of files)
    dic_length_glob = {}
    if counting_nucleotides:
        dic_length_glob = prediction.recovery_number_nucleotides(
            diamond_df_list, list_fastq
        )
        prediction.writing_input_files_with_number_nucleotides(
            list_fastq, dic_length_glob
        )
    else:
        dic_length_glob = prediction.computing_number_nucleotides_multi_readsets(
            diamond_df_list, list_fastq, number_nt_list
        )
        prediction.writing_input_files_with_number_nucleotides(
            list_fastq, dic_length_glob
        )

    # Predicts sample size
    prediction.prediction_size_sample(
        dic_length_glob,
        slope,
        os.path.abspath("filtered_sample/df_with_sample_and_coverage.tsv"),
        args.picog,
    )

    # Cleaning the output dir
    if args.clean:
        shutil.rmtree("best_hit_per_read")
        shutil.rmtree("cmds")
        shutil.rmtree("deviant_genes")
        shutil.rmtree("ref_prot_length")
        shutil.rmtree("Sample_mapped", ignore_errors=True)
        shutil.rmtree("Trimmed_directory", ignore_errors=True)

    # END
    print(
        f"The depth can be found in {args.output_dir}/filtered_sample/df_with_sample_and_coverage.tsv"
    )
    print(
        f"\nGenome size prediction(s) can be found in {args.output_dir}/Sample_Size/samples_sizes.tsv"
    )
    print(
        f"\n Total running time : {float(time.perf_counter() - global_start)} seconds"
    )


if __name__ == "__main__":
    run()
