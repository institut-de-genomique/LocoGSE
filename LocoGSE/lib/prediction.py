import os
import pandas as pd
from Bio import SeqIO
import gzip


def determine_slope_for_family(family: str, use_busco: bool) -> float:
    path_module = os.path.abspath(__file__)
    if not use_busco:
        path_database = path_module.replace(
            "prediction.py", "PlantFamilies.CoeffRegression.V2.txt"
        )
    else:
        path_database = path_module.replace(
            "prediction.py", "PlantFamilies.CoeffRegression.BUSCO.V2.txt"
        )
    database_organism = pd.read_csv(path_database, sep="\t", header=None)
    line_family = database_organism.loc[database_organism[0] == family].index.values
    slope = database_organism.iloc[line_family, 2]
    return list(slope)[0]


def determine_slope_for_lineage(lineage: str, use_busco: bool) -> float:
    path_module = os.path.abspath(__file__)
    if not use_busco:
        path_database = path_module.replace(
            "prediction.py", "PlantFamilies.CoeffRegression.V2.txt"
        )
    else:
        path_database = path_module.replace(
            "prediction.py", "PlantFamilies.CoeffRegression.BUSCO.V2.txt"
        )
    database_organism = pd.read_csv(path_database, sep="\t", header=None)
    line_family = database_organism.loc[database_organism[1] == lineage].index.values
    print(line_family)
    slope = database_organism.iloc[line_family[0], 2]
    print(slope)
    return list(slope)[0]


def computing_number_of_nucleotides_fastq(fastq_file: str) -> dict:
    total = 0
    dic_length = {}
    if fastq_file.endswith(".gz"):
        fastq_in = gzip.open(fastq_file, "rt")
    else:
        fastq_in = fastq_file
    for seq_record in SeqIO.parse(fastq_in, "fastq"):
        total += len(seq_record)
    dic_length["sorted_" + os.path.basename(fastq_file) + ".out"] = total
    return dic_length


def computing_number_nucleotides_multi_readsets(
    diamond_df_list: list, multi_reads: list, number_nt_list: list
) -> dict:
    dic_length_glob = {}
    for specie in range(0, len(diamond_df_list)):
        common_name = os.path.basename(diamond_df_list[specie])
        dic_length_specie = {}
        dic_length_specie["sorted_" + common_name] = 0
        if number_nt_list[specie] != "No_number":
            dic_length_glob["sorted_" + common_name] = number_nt_list[specie]
        else:
            for sample in range(0, len(multi_reads[specie])):
                dict_length_sample = {}
                dict_length_sample.update(
                    computing_number_of_nucleotides_fastq(multi_reads[specie][sample])
                )
                dic_length_specie["sorted_" + common_name] += dict_length_sample[
                    "sorted_" + os.path.basename(multi_reads[specie][sample]) + ".out"
                ]
            dic_length_glob["sorted_" + common_name] = dic_length_specie[
                "sorted_" + common_name
            ]
    return dic_length_glob


def number_nucleotides_in_df(ref_nucleotides: str, diamond_df_list: list) -> dict:
    dic_length_glob = {}
    ref_nucleotides_df = pd.read_csv(ref_nucleotides, sep="\t", header=None)
    for length in range(0, ref_nucleotides_df.shape[0]):
        for sample in range(0, len(diamond_df_list)):
            if (
                ref_nucleotides_df.iloc[length][0] in diamond_df_list[sample]
                and "sorted_" + os.path.basename(diamond_df_list[sample])
                in dic_length_glob
            ):
                dic_length_glob[
                    "sorted_" + os.path.basename(diamond_df_list[sample])
                ] += ref_nucleotides_df.iloc[length][1]
            if (
                ref_nucleotides_df.iloc[length][0] in diamond_df_list[sample]
                and "sorted_" + os.path.basename(diamond_df_list[sample])
                not in dic_length_glob
            ):
                dic_length_glob[
                    "sorted_" + os.path.basename(diamond_df_list[sample])
                ] = ref_nucleotides_df.iloc[length][1]
    return dic_length_glob


def prediction_size_sample(
    dic_length_glob: dict, slope: float, file_coverage: str, picog: str
) -> None:
    f_cov = pd.read_csv(file_coverage, header=None, sep="\t")
    dic_sample_size = {}
    for sample in list(dic_length_glob.keys()):
        line_sample = f_cov.loc[f_cov[0] == sample].index.values
        coverage = f_cov.loc[line_sample[0]][1]
        if picog != "n":
            size = dic_length_glob[sample] / (float(slope) * coverage * 1000000 * 978)
        else:
            size = dic_length_glob[sample] / (float(slope) * coverage * 1000000)
        dic_sample_size[os.path.basename(sample)] = size
    df_size = pd.DataFrame.from_dict(dic_sample_size, orient="index")
    df_size.to_csv(path_or_buf="Sample_Size/samples_sizes.tsv", sep="\t", header=False)


def writing_input_files_with_number_nucleotides(
    multi_files: list, dic_length_glob: dict
) -> None:
    new_input_file = open("Sample_number_nucleotides/New_Input_File.txt", "w+")

    for nbnt in range(0, len(dic_length_glob.values())):
        input_sample = " ".join(multi_files[nbnt])
        new_input_file.write(
            input_sample + " " + str(list(dic_length_glob.values())[nbnt]) + "\n"
        )
    new_input_file.close()


def recovery_number_nucleotides(diamond_df_list: list, multi_reads: list) -> dict:
    file_with_nbnt = open("Sample_number_nucleotides/New_Input_File.txt", "r")
    line_nbnt = file_with_nbnt.readlines()
    dic_length_glob = {}
    for specie in range(0, len(diamond_df_list)):
        common_name = os.path.basename(diamond_df_list[specie])
        dic_length_specie = {}
        dic_length_specie["sorted_" + common_name] = 0
        if specie < len(line_nbnt):
            list_line_nbnt = line_nbnt[specie].split()
            if list_line_nbnt[-1].isnumeric():
                dic_length_glob["sorted_" + common_name] = int(list_line_nbnt[-1])
            else:
                for sample in range(0, len(multi_reads[specie])):
                    dict_length_sample = {}
                    dict_length_sample.update(
                        computing_number_of_nucleotides_fastq(
                            multi_reads[specie][sample]
                        )
                    )
                    dic_length_specie["sorted_" + common_name] += dict_length_sample[
                        "sorted_"
                        + os.path.basename(multi_reads[specie][sample])
                        + ".out"
                    ]
                dic_length_glob["sorted_" + common_name] = dic_length_specie[
                    "sorted_" + common_name
                ]
        else:
            for sample in range(0, len(multi_reads[specie])):
                dict_length_sample = {}
                dict_length_sample.update(
                    computing_number_of_nucleotides_fastq(multi_reads[specie][sample])
                )
                dic_length_specie["sorted_" + common_name] += dict_length_sample[
                    "sorted_" + os.path.basename(multi_reads[specie][sample]) + ".out"
                ]
            dic_length_glob["sorted_" + common_name] = dic_length_specie[
                "sorted_" + common_name
            ]
    return dic_length_glob
