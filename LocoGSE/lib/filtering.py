import os
from Bio import SeqIO
import statistics
import pandas as pd


def search_best_hit_per_read(df_diamond: str) -> None:
    name_file = os.path.basename(df_diamond)
    df = pd.read_csv(df_diamond, sep="\t", header=None)
    sorted_df = df.sort_values(by=[10])
    sorted_df = sorted_df.astype(str)
    final_df = sorted_df.drop_duplicates(subset=0, keep="first")
    path_for_file = "best_hit_per_read/sorted_" + name_file
    final_df.to_csv(path_or_buf=path_for_file, sep="\t", header=False, index=False)


def computing_length_prot_database(ref_prot: str) -> None:
    file_length = open("ref_prot_length/ref_prot.tsv", "w")
    for seq_record in SeqIO.parse(ref_prot + ".fa", "fasta"):
        file_length.write(seq_record.id)
        file_length.write("\t")
        file_length.write(str(len(seq_record)))
        file_length.write("\n")
    file_length.close()


def writing_depth_for_each_protein(dict_depth_prot: dict, name_file: str) -> None:
    file_coverage_per_specie = open("filtered_sample/" + name_file, "w")
    file_coverage_per_specie.write(name_file+"\n")
    for name, depth in dict_depth_prot.items():
        file_coverage_per_specie.write(str(name) + "\t" + str(depth) + "\n")
    file_coverage_per_specie.close()


def filter_sample(prot_with_length: str, df_besthit: str) -> None:
    ref = pd.read_csv(prot_with_length, sep="\t", header=None)
    dic_sample = {}
    list_deviant_genes = []

    for file in range(0, len(df_besthit)):
        list_deviant_genes.append(os.path.basename(df_besthit[file]))
        df = pd.read_csv(df_besthit[file], sep="\t", header=None)
        lgprot_length_total = 0
        dic_length_prot = {}
        dic_depth_prot = {}
        for _, rows in ref.iterrows():
            line_sample = df.loc[df[1] == rows[0]].index.values
            length_per_prot = 0
            for x in line_sample:
                length_per_prot += df.iloc[x, 9] - df.iloc[x, 8] + 1
            lgprot_length_total += rows[1]
            dic_length_prot[rows[0]] = length_per_prot
            dic_depth_prot[rows[0]] = length_per_prot / rows[1]

        mean = sum(dic_length_prot.values()) / lgprot_length_total
        st_dev = statistics.pstdev(dic_length_prot.values())

        good_genes = []
        deviant_genes = []
        for z in dic_length_prot.keys():
            z_score = (dic_length_prot[z] - mean) / st_dev
            if z_score >= 1.96 or z_score <= (-1.96):
                deviant_genes.append(z)
            else:
                good_genes.append(z)
        for x in deviant_genes:
            dic_length_prot.pop(x)
            dic_depth_prot.pop(x)
            line_deviant = ref.loc[ref[0] == x].index.values
            lgprot_length_total -= ref.iloc[line_deviant[0], 1]
            list_deviant_genes.append(str(x))
        writing_depth_for_each_protein(
            dic_depth_prot, os.path.basename(df_besthit[file])
        )
        filtered_mean = sum(dic_length_prot.values()) / lgprot_length_total
        dic_sample[os.path.basename(df_besthit[file])] = str(filtered_mean)

    deviant_file = open("deviant_genes/deviant_file_liste.txt", "w")
    l = "\n".join(list_deviant_genes)
    deviant_file.write(l)
    deviant_file.close()
    final_df = pd.DataFrame.from_dict(dic_sample, orient="index")
    final_df.to_csv(
        path_or_buf="filtered_sample/df_with_sample_and_coverage.tsv",
        sep="\t",
        header=False,
    )


def searching_besthit(diamond_df_list: list) -> list:
    besthit = []
    for mapped_file in diamond_df_list:
        mapped_file_name = os.path.basename(mapped_file)
        search_best_hit_per_read(mapped_file)
        name_sorted_best_hit = ["best_hit_per_read/sorted_", mapped_file_name]
        besthit.append(os.path.abspath("".join(name_sorted_best_hit)))
    return besthit
