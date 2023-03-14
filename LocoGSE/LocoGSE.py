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


def run() :
    parser = argparse.ArgumentParser(
        prog="LocoGSE : Low coverage Genome Size Estimator",
        description="\n\n A Genome Size Estimation program. It is based on a linear relation between the sequencing depth (linked to the genome size) and the depth on a set of monocopy genes. \n The regression factor (slope) depends of the plant family/lineage : slopes are precomputed for a number of plant families. \n For questions : https://github.com/institut-de-genomique/LocoGSE/issues or pierre.guenzi.tiberi@gmail.com or fdenoeud@genoscope.cns.fr",
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=True,
    )

    mandatory_args = parser.add_argument_group("Mandatory arguments")

    # Mandatory arguments

    # Fastq file (short reads)
    mandatory_args.add_argument(
        "--reads",
        action="store",
        nargs="*",
        dest="reads",
        help="Directory with all fastq.gz files in single-end",
        default="",
        required=False,
    )

    # FASTA file with all reference proteins
    mandatory_args.add_argument(
        "--ref_prot",
        action="store",
        dest="ref",
        help="Path to a monocopy protein database to be used with diamond (the directory must contain a .fa file and a .dmnd file). The PREFIX (file name without the extension) of these 2 files must be given. By default : OneKP consensus obtained from https://github.com/smirarab/1kp/tree/master/alignments ! Caution: the option --lgprot needs to be provided too!",
        default="",
        required=False,
    )
    # Multiple files = list of absolute paths to several fastq files (if you put several files on one line they will be treated together: for instance read1.fastq and read2.fastq)
    mandatory_args.add_argument(
        "--list_fastq",
        action="store",
        dest="list_fastq",
        help="txt file with, on each line, the list of fastq files to be treated together (same sample)(for example: name_sample mypath/read1.fastq mypath/read2.fastq)" ,
        default="",
        required=False,
    )

    # Optional arguments

    optional_args = parser.add_argument_group("Optional arguments")

    # Recovery
    optional_args.add_argument(
        "--recovery",
        "-r",
        action="store_true",
        dest="recovery",
        help="Recovery option to continue the run started in the output directory provided !",
        default=None,
        required=False,
    )

    # Threads
    optional_args.add_argument(
        "--threads",
        "-t",
        action="store",
        dest="threads",
        help="Number of threads to be used during the mapping step",
        default=1,
        required=False,
    )

    # Slope
    optional_args.add_argument(
        "--slope",
        "-s",
        action="store",
        dest="slope",
        help="Slope (regression factor) used to estimate sequencing depth from depth on monocopy proteins. It is specific to each plant lineage. Pre computed slopes are available for families listed in --list_families and lineages in --list_lineagges : no need to provide a slope if your species is in the list, you can just provide either the family or the lineage", 
        default="",
        required=False,
    )

    # Family
    optional_args.add_argument(
        "--family",
        "-f",
        action="store",
        dest="family",
        help="Specify the plant family in order to use a pre-computed slope",
        default="",
        required=False,
    )

    # List of families
    optional_args.add_argument(
        "--list_families",
        action="store_true",
        dest="list_families",
        help="Print all families with available pre-computed slope",
        default=None,
        required=False,
    )

    # Lineage
    optional_args.add_argument(
        "--lineage",
        action="store",
        dest="lineage",
        help="Specify the plant lineage in order to use a pre-computed slope",
        default="",
        required=False,
    )

    # List of lineages
    optional_args.add_argument(
        "--list_lineages",
        action="store_true",
        dest="list_lineages",
        help="Print all plant lineages with available pre-computed slope",
        default=None,
        required=False,
    )

    # Sequence length
    optional_args.add_argument(
        "--length_trim",
        action="store",
        dest="length_sequence",
        help="Cumulative length of all sequences treated after trimming step (by default : 100 # Highly recommended since the training step was performed with 100nt reads !)",
        default="100",
        required=False,
    )

    # Trimming
    optional_args.add_argument(
        "--no_trim",
        action="store_true",
        dest="no_trim",
        help="Desactivates the trimming step (only if your reads are 100 nt long. Otherwise, highly recommended since the training step was performed with 100nt reads!)",
        default=None,
        required=False,
    )

    # Size in pg
    optional_args.add_argument(
        "--picog",
        action="store",
        dest="picog",
        help="An option to convert genome size to picograms (in MB by default)",
        default="n",
        required=False,
    )

    # Size of each reference protein, if the user provides another protein set (not OneKp)
    optional_args.add_argument(
        "--lgprot",
        "-l",
        action="store",
        dest="lgprot",
        help="A tabulated file with each protein name(1st column) and its length(2nd column) in aa, if the protein database specified is not OneKp (default)",
        default="",
        required=False,
    )

    # Pegasus option
    optional_args.add_argument(
        "--pegasus",
        action="store_true",
        dest="pegasus",
        help="Writes a Pegasus script that can be manually launched to submit multiple mapping commands at the same time",
        default=None,
        required=False,
    )

    # Output
    optional_args.add_argument(
        "--output",
        "-o",
        action="store",
        dest="output_dir",
        help="Output directory name",
        default="results",
        required=False,
    )

    # Clean output with only sample size
    optional_args.add_argument(
        "--cleaning_output",
        action="store_true",
        dest="clean",
        help="Remove temporary files and only keep results",
        default="n",
        required=False,
    )

    args = parser.parse_args()

    if args.list_families:
        checking.list_families_print()
        sys.exit(-1)

    if args.list_lineages:
        checking.list_lineages_print()
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
        number_nt_list, list_fastq, name_samples = checking.complete_single_files(args.reads)

    # Checking
    if args.ref == "":
        print("No argument given to --ref_prot, defaulting to the OneKP database")
        path_main = os.path.abspath(__file__)
        args.ref = path_main.replace("LocoGSE.py", "OneKP.410genes.consensus")
    else:
        args.ref = os.path.abspath(args.ref)

    if not args.recovery:
        try:
            os.mkdir(args.output_dir)
        except:
            print(
                f"\n Output directory {args.output_dir} can not be created, please erase it before launching the program !"
            )
            sys.exit(1)

    if args.lgprot != "":
        lgprot_path = os.path.abspath(args.lgprot)

    if args.slope == "":
        if args.family == "":
            if args.lineage == "":
                print(
                    "Warning: No family or lineage given to --family or --lineage, defaulting to slope=1. List of families and lineages pre-computed are available with --list_families and --list_lineages "
                )
                slope = 1
                no_slope = True
            else:
                slope = prediction.determine_slope_for_lineage(args.lineage)
                no_slope = False
        else:
            slope = prediction.determine_slope_for_family(args.family)
            no_slope = False
    else:
        slope = args.slope
        no_slope = False

    actual_path = os.getcwd()
    os.chdir(args.output_dir)
    global_start = time.perf_counter()

    # Checking steps
    mapping_step_finished = checking.checking_step_recovery_mapping()
    counting_nucleotides = checking.checking_step_recovery_computing_nbnt()
    checking.checking_step(args.ref, slope, mapping_step_finished, args.length_sequence)
    slope = float(slope)

    # Write pegasus script
    if args.pegasus:
        writing_script(
            list_fastq, args.threads, args.ref, name_samples, slope, actual_path
        )
        sys.exit(-1)

    # Trimming step
    if not args.no_trim:
        list_fastq = trimming.trimming_step(
            list_fastq, args.length_sequence, args.threads
        )
    else:
        # No trimming
        print("Any trimming")

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
    if counting_nucleotides == True:
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

    if no_slope == True:
        print(f"The depth can be found in {args.output_dir}filtered_sample/df_with_sample_and_coverage.tsv")
        print("\nA slope is needed to predict the genome sample size \n")
        print(
            "\n Please provide a family or lineage in order to use a pre-computed slope or calculate a slope as indicated in the wiki: https://github.com/institut-de-genomique/LocoGSE/wiki/2.Linear-regression \n"
        )
        sys.exit(-1)

    # Predicts sample size
    prediction.prediction_size_sample(
        dic_length_glob,
        slope,
        os.path.abspath("filtered_sample/df_with_sample_and_coverage.tsv"),
        args.picog,
    )

    # Cleaning the output dir
    if args.clean != "n":
        shutil.rmtree("best_hit_per_read")
        shutil.rmtree("cmds")
        shutil.rmtree("ref_prot_length")
        shutil.rmtree("Sample_mapped", ignore_errors=True)

    # END
    print(f"The depth can be found in {args.output_dir}filtered_sample/df_with_sample_and_coverage.tsv")
    print(f"\n Results can be found in {args.output_dir}Sample_Size/samples_sizes.tsv")
    print(
        f"\n Total running time : {float(time.perf_counter() - global_start)} seconds"
    )

if __name__ == '__main__':
    run()
