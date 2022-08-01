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


def check_length_sequence(length_sequence: str) -> None:
    try:
        val = int(length_sequence)
    except ValueError:
        print("ERROR : the length of sequence is not a number")
        sys.exit(-1)


def checking_step(
    ref: str, slope: str, mapping_step_finished: bool, length_sequence: str
) -> None:
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
    check_length_sequence(length_sequence)


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


def checking_number_name_samples(name_samples: list, multi_files: list) -> list:
    if name_samples[0] == "Any names":
        name_samples.clear()
        for specie in range(0, len(multi_files)):
            name_samples.append("input_" + str(specie))
        return name_samples
    if len(name_samples) < len(multi_files):
        for x in range(0, len(multi_files) - len(name_samples)):
            name_samples.append("input_" + str(x))
        return name_samples
    elif len(name_samples) > len(multi_files):
        reducted_name_samples = name_samples[0 : len(multi_files) - len(name_samples)]
        return reducted_name_samples
    else:
        return name_samples


def complete_multi_files(multi_files_a: str) -> tuple:
    final_multi_files = []
    final_number_nt_list = []
    multi = open(os.path.abspath(multi_files_a), "r")
    lines_multi = multi.readlines()
    multi.close()
    for line in range(0, len(lines_multi)):
        list_line = lines_multi[line].split()
        if list_line[-1].isnumeric():
            final_number_nt_list.append(int(list_line[-1]))
            final_multi_files.append(list_line[0 : len(list_line) - 1])
        else:
            final_number_nt_list.append("No_number")
            final_multi_files.append(list_line)
    return final_number_nt_list, final_multi_files


def complete_single_files(reads: list) -> list:
    multi_files = []
    number_nt_list = []
    for read in range(0, len(reads)):
        multi_files.append(list(""))
        multi_files[read].append(reads[read])
        number_nt_list.append("No_number")
    return number_nt_list, multi_files

def list_families_print():
    path_module = os.path.abspath(__file__)
    path_database = path_module.replace(
        "checking.py", "PlantFamilies.CoeffRegression.txt"
    )
    database_organism = pd.read_csv(path_database, sep="\t", header=None)
    for family in database_organism[0]:
        print(family)

def list_lineages_print():
    path_module = os.path.abspath(__file__)
    path_database = path_module.replace(
        "checking.py", "PlantFamilies.CoeffRegression.txt"
    )
    database_organism = pd.read_csv(path_database, sep="\t", header=None)
    list_lineage_save=database_organism[1].unique()
    for lineage in list_lineage_save:
        print(str(lineage))
