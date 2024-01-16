import os
import sys
import subprocess
import pandas as pd

def check_dependencies() -> None:
    print("\n Checking dependencies ...", flush=True)
    FNULL = open(os.devnull, "w")

    # Looking for Diamond
    try:
        subprocess.call(["diamond", "view"], stdout=FNULL, stderr=FNULL)
        print("\t Found Diamond.", flush=True)
    except OSError:
        print("\t ERROR : Diamond not found. ", flush=True)
        sys.exit(-1)

    try:
        subprocess.call(["cutadapt"], stdout=FNULL, stderr=FNULL)
        print("\t Found cutadapt", flush=True)
    except OSError:
        print("\t ERROR : Cutadapt not found.", flush=True)
        sys.exit(-1)


def check_fasta_prot(prot: list) -> bool:
    authorized_chars = (
        "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ:@._"
    )
    for line in open(prot + ".fa"):
        if line.startswith(">"):
            header = line[1:].rstrip("\n")
            for char in header:
                if char not in authorized_chars:
                    return True
    return False


def check_creating_directory() -> None:
    try:
        os.mkdir("Sample_mapped")
        os.mkdir("deviant_genes")
        os.mkdir("best_hit_per_read")
        os.mkdir("cmds")
        os.mkdir("ref_prot_length")
        os.mkdir("filtered_sample")
        os.mkdir("Sample_Size")
        os.mkdir("Sample_number_nucleotides")
        os.mkdir("Trimmed_directory")
    except:
        print("ERROR : Could not create subdirectories")
        sys.exit(-1)


def check_slope_type(slope: str) -> None:
    try:
        val = float(slope)
    except ValueError:
        print("ERROR : The slope is not a number")
        sys.exit(-1)


def checking_step(
    ref: str, slope: str, mapping_step_finished: bool) -> None:
    check_dependencies()
    # Check the reference fasta file
    invalid_chars = check_fasta_prot(ref)
    if invalid_chars == True:
        print("Invalid characters were detected in fasta reference --ref_p")
        print(
            "Valid characters are : 0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ:@ "
        )
        sys.exit(-1)
    if mapping_step_finished != True:
        check_creating_directory()
    # Check the slope reference
    check_slope_type(slope)


def checking_step_recovery_mapping() -> bool:
    if os.path.exists("Sample_mapped/finished.txt"):
        print("Mapping step is already finished")
        mapping_step_finished = True
    else:
        mapping_step_finished = False

    return mapping_step_finished


def checking_step_recovery_computing_nbnt() -> bool:
    if os.path.exists("Sample_number_nucleotides/New_Input_File.txt"):
        print("Counting of nucleotides is already finished")
        counting_nucleotides = True
    else:
        counting_nucleotides = False
    return counting_nucleotides

def complete_multi_files(multi_files_a: str) -> tuple:
    final_multi_files = []
    final_number_nt_list = []
    final_name_samples_list = []
    multi = open(os.path.abspath(multi_files_a), "r")
    lines_multi = multi.readlines()
    multi.close()
    for line in range(0, len(lines_multi)):
        list_line = lines_multi[line].split()
        if not list_line[0].endswith(".gz") or not list_line[0].endswith(".fq") or not list_line[0].endswith(".fastq") :
             final_name_samples_list.append(list_line[0])
        else:
            final_name_samples_list.append("input_"+str(line))
        if list_line[-1].isnumeric():
            final_number_nt_list.append(int(list_line[-1]))
            final_multi_files.append(list_line[1 : len(list_line) - 1])
        else:
            final_number_nt_list.append("No_number")
            final_multi_files.append(list_line[1 : len(list_line)])
    return final_number_nt_list, final_multi_files, final_name_samples_list


def complete_single_files(reads: list) -> list:
    multi_files = []
    name_samples = []
    number_nt_list = []
    for read in range(0, len(reads)):
        multi_files.append(list(""))
        multi_files[read].append(reads[read])
        number_nt_list.append("No_number")
        name_samples.append("input_"+str(read))
    return number_nt_list, multi_files, name_samples

def list_families_print(use_busco):
    path_module = os.path.abspath(__file__)

    if not use_busco:
        path_database = path_module.replace(
            "checking.py", "PlantFamilies.CoeffRegression.txt"
        )
    else:
        path_database = path_module.replace(
            "checking.py", "PlantFamilies.CoeffRegression.BUSCO.txt"
        )

    database_organism = pd.read_csv(path_database, sep="\t", header=None)
    for family in database_organism[0]:
        print(family)

def list_lineages_print(use_busco):
    path_module = os.path.abspath(__file__)

    if not use_busco:
        path_database = path_module.replace(
            "checking.py", "PlantFamilies.CoeffRegression.txt"
        )
    else:
        path_database = path_module.replace(
            "checking.py", "PlantFamilies.CoeffRegression.BUSCO.txt"
        )
    database_organism = pd.read_csv(path_database, sep="\t", header=None)
    list_lineage_save=database_organism[1].unique()
    for lineage in list_lineage_save:
        print(str(lineage))
