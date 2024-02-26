import subprocess
import os
import time
import pandas as pd


def launch_mapping(read: str, file_name: str, ref_prot: str, threads: int) -> None:
    ################DIAMOND MAPPING######################
    print(read)
    print("\n Generating diamond files ...", flush=True)
    cmd = [
        "diamond blastx -p",
        str(threads),
        "-q",
        read,
        "-o",
        "Sample_mapped/" + file_name + ".out",
        "-d",
        ref_prot,
        "-e 0.00001 -f 6 ",
    ]

    start = time.perf_counter()
    print(" ".join(cmd), flush=True, file=open("cmds/diamond.cmds", "w"))

    try:
        _ = subprocess.run(args=" ".join(cmd), shell=True, check=True)
    except Exception as error:
        print("\n Error : Couldn't mapping fastq", flush=True)
        print(error)
        exit(1)
    print(f"\n Done in {float(time.perf_counter() - start)} seconds", flush=True)


def concatenate_multi_reads_diamond_results(list_df: list, name_samples: str) -> str:
    result_df = pd.DataFrame(columns=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    for sample in list_df:
        sample_df = pd.read_csv(sample, sep="\t", header=None)
        # Rename the read so there is no collision between files
        sample_df[0] = sample_df[0].map(lambda x: str(x) + "_" + sample.split("/")[-1].split(".")[0])
        result_df = pd.concat([result_df, sample_df], ignore_index=True, sort=False)
    result_df.to_csv(
        path_or_buf="Sample_mapped/" + name_samples + "_merged.tsv",
        sep="\t",
        header=False,
        index=False,
    )
    return name_samples + "_merged.tsv"


def mapping_multi_files(
    multi_reads: list, ref: str, threads: int, name_samples: list
) -> list:
    diamond_df_list = []
    for sample in range(0, len(multi_reads)):
        sample_name = []
        for file in multi_reads[sample]:
            file_name = os.path.basename(file)
            sample_name.append(os.path.join("Sample_mapped", file_name + ".out"))
            launch_mapping(os.path.abspath(file), file_name, ref, threads)
        diamond_df_list.append(
            os.path.abspath(
                os.path.join(
                    "Sample_mapped",
                    concatenate_multi_reads_diamond_results(
                        sample_name, name_samples[sample]
                    ),
                )
            )
        )
    file_finished = open("Sample_mapped/finished.txt", "w+")
    file_finished.write("Mapping step is finished")
    file_finished.close()
    return diamond_df_list


def recovery_mapping_step(
    multi_reads: list, ref: str, threads: int, name_samples: list
) -> list:
    diamond_df_list = []
    for name in range(0, len(name_samples)):
        if os.path.exists("Sample_mapped/" + name_samples[name] + "_merged.tsv"):
            diamond_df_list.append(
                os.path.abspath("Sample_mapped/" + name_samples[name] + "_merged.tsv")
            )
        else:
            sample_name = []
            for file in multi_reads[name]:
                file_name = os.path.basename(file)
                sample_name.append(os.path.join("Sample_mapped", file_name + ".out"))
                launch_mapping(os.path.abspath(file), file_name, ref, threads)
            diamond_df_list.append(
                os.path.abspath(
                    os.path.join(
                        "Sample_mapped",
                        concatenate_multi_reads_diamond_results(
                            sample_name, name_samples[name]
                        ),
                    )
                )
            )

    return diamond_df_list
