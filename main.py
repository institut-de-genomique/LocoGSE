#!/usr/bin/env python3

from lib import prediction
from lib import mapping
from lib import filtering
from lib import checking
from lib import trimming
import pegasus
import argparse
import os
import sys
import time
import shutil


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="LoCoGSe : Low Coverage Genome Size Estimation",
        description="\n\n A Genome Size Estimation program. It is based on a linear relation between the depth and the genome size. \n A correction which depends of the family is added at this depth for a best prediction. \n For a question : https://github.com/institut-de-genomique/LocoGSE/issues or pierre.guenzi.tiberi@gmail.com",
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=True,
    )

    mandatory_args = parser.add_argument_group("Mandatory arguments")

    # Mandatory arguments

    # Fastq files with reads
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
        help="Path to a protein database to be used with diamond (the directory must contain a .fa file and a .dmnd file). The prefix of these 2 files must be given. By default : OneKP consensus for the plants",
        default="",
        required=False,
    )
    # Multiple files
    mandatory_args.add_argument(
        "--multi_files",
        action="store",
        dest="multi_files",
        help="txt file with, in each line, each path file for the same sample",
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
        help="Recovery option to continue a previous run after main steps",
        default=None,
        required=False,
    )

    # Threads
    optional_args.add_argument(
        "--threads",
        "-t",
        action="store",
        dest="threads",
        help="Number of steps to be used during the mapping step",
        default=1,
        required=False,
    )

    # Slope
    optional_args.add_argument(
        "--slope",
        "-s",
        action="store",
        dest="slope",
        help="A slope which can be take to estimate the genome size. It is specific to each family",
        default="",
        required=False,
    )

    # Family
    optional_args.add_argument(
        "--family",
        "-f",
        action="store",
        dest="family",
        help="Some families are in the database with this program. Specify the family of the sample.",
        default="",
        required=False,
    )

    # List of families
    optional_args.add_argument(
        "--list_families",
        action="store_true",
        dest="list_families",
        help="Print all families which are treated by LocoGSE",
        default=None,
        required=False,
    )

    # Lineage
    optional_args.add_argument(
        "--lineage",
        action="store",
        dest="lineage",
        help="Some lineages are in the database with this program. Specify the lineage of the sample.",
        default="",
        required=False,
    )

    # List of lineages
    optional_args.add_argument(
        "--list_lineages",
        action="store_true",
        dest="list_lineages",
        help="Print all lineages which are treated by LocoGSE",
        default=None,
        required=False,
    )

    # Sequence length
    optional_args.add_argument(
        "--length_trim",
        action="store",
        dest="length_sequence",
        help="Length (INT) of the sequence after cutting step (by default : 100)",
        default="100",
        required=False,
    )

    # Trimming
    optional_args.add_argument(
        "--no_trim",
        action="store_true",
        dest="no_trim",
        help="Desactivates the trimming step",
        default=None,
        required=False,
    )

    # Name of sample
    optional_args.add_argument(
        "--name_samples",
        action="store",
        dest="name_samples",
        help="A TSV with, in each line, the name of sample in the same line that in multi_files",
        default="",
        required=False,
    )

    # Size in pg
    optional_args.add_argument(
        "--picog",
        action="store",
        dest="picog",
        help="An option to convert genome size in picogram (in MB by default)",
        default="n",
        required=False,
    )

    # Size of each reference protein
    optional_args.add_argument(
        "--lgprot",
        "-l",
        action="store",
        dest="lgprot",
        help="A TSV file with each protein name(1st column) and its length(2nd column)",
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
    multi_files = []
    if args.reads == "":
        if args.multi_files != "":
            number_nt_list, multi_files = checking.complete_multi_files(
                args.multi_files
            )
        else:
            print("Please provide samples with either --reads or --multi_files.")
            sys.exit(-1)
    else:
        number_nt_list, multi_files = checking.complete_single_files(args.reads)

    # Checking
    if args.ref == "":
        print("No argument given to --ref_prot, defaulting to the OneKP database")
        path_main = os.path.abspath(__file__)
        args.ref = path_main.replace("main.py", "lib/OneKP.410genes.consensus")
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

    if args.slope == "":
        if args.family == "":
            if args.lineage == "":
                print(
                    "No family or lineage given to --family or --lineage, defaulting to slope=1"
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

    if args.name_samples != "":
        file_name = open(os.path.abspath(args.name_samples), "r")
        lines_names = file_name.readlines()
        file_name.close()
        name_samples = []
        for line in range(0, len(lines_names)):
            name_samples.append(lines_names[line][0 : len(lines_names[line]) - 1])
    else:
        name_samples = ["Any names"]

    actual_path = os.getcwd()
    os.chdir(args.output_dir)
    global_start = time.perf_counter()

    # Checking steps
    mapping_step_finished = checking.checking_step_recovery_mapping()
    counting_nucleotides = checking.checking_step_recovery_computing_nbnt()
    checking.checking_step(args.ref, slope, mapping_step_finished, args.length_sequence)
    name_samples = checking.checking_number_name_samples(name_samples, multi_files)
    slope = float(slope)

    # Writing pegasus script
    if args.pegasus:
        pegasus.writing_script(
            multi_files, args.threads, args.ref, name_samples, slope, actual_path
        )
        sys.exit(-1)

    # Trimming step
    if not args.no_trim:
        multi_files = trimming.trimming_step(
            multi_files, args.length_sequence, args.threads
        )
    else:
        # No trimming
        print("Any trimming")

    # Mapping reads on protein
    diamond_df_list = []
    if mapping_step_finished != True:
        diamond_df_list = mapping.mapping_multi_files(
            multi_files, args.ref, args.threads, name_samples
        )
    else:
        diamond_df_list = mapping.recovery_mapping_step(
            multi_files, args.ref, args.threads, name_samples
        )

    # Find one best hit per read in diamond DataFrame
    besthit = []
    besthit = filtering.searching_besthit(diamond_df_list)

    # Compute the length of each protein

    if not args.lgprot:
        filtering.computing_length_prot_database(args.ref)
        lgprot_path = "ref_prot_length/ref_prot.tsv"

    # Find deviants genes and outliers and creating of a filtered file

    filtering.filter_sample(lgprot_path, besthit)

    # Compute the number of nucleotides in each FASTQ Files
    dic_length_glob = {}
    if counting_nucleotides == True:
        dic_length_glob = prediction.recovery_number_nucleotides(
            diamond_df_list, multi_files
        )
        prediction.writing_input_files_with_number_nucleotides(
            multi_files, dic_length_glob
        )
    else:
        dic_length_glob = prediction.computing_number_nucleotides_multi_readsets(
            diamond_df_list, multi_files, number_nt_list
        )
        prediction.writing_input_files_with_number_nucleotides(
            multi_files, dic_length_glob
        )

    if no_slope == True:
        print("\nA slope is needed to predict the genome sample size \n")
        print(
            "\n Please compute the slope as indicated in the wiki: https://github.com/institut-de-genomique/LocoGSE/wiki/2.Linear-regression \n"
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
        shutil.rmtree("filtered_sample")
        shutil.rmtree("ref_prot_length")
        shutil.rmtree("Sample_mapped", ignore_errors=True)

    # END
    print(f"\n Results can be found in {args.output_dir}")
    print(
        f"\n Total running time : {float(time.perf_counter() - global_start)} seconds"
    )
