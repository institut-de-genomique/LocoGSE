import subprocess
import os
import time


def trimming_cutadapt(fastq: str, threads: int) -> None:
    ##################CUTADAPT TRIMMING#################
    print(fastq)
    print("Generating trimmed fastq ...", flush=True)

    cmd = [
        "cutadapt -o Trimmed_directory/trimmed_" + os.path.basename(fastq),
        "-l",
        "100",
        fastq,
        "-j",
        str(threads),
    ]

    start = time.perf_counter()
    print(" ".join(cmd), flush=True, file=open("cmds/trimming.cmds", "w"))

    try:
        _ = subprocess.run(args=" ".join(cmd), shell=True, check=True)
    except Exception as error:
        print("\n Error : Couldn't trim fastq", flush=True)
        print(error)
        exit(1)

    print(f"\n Done in {float(time.perf_counter() - start)} secondes", flush=True)


def trimming_step(multi_files: list, threads: int) -> list:
    multi_trimmed_files = []
    for specie in range(0, len(multi_files)):
        multi_trimmed_files.append(list(""))
        for sample in range(0, len(multi_files[specie])):
            if (
                os.path.exists(
                    os.path.abspath(
                        "Trimmed_directory/trimmed_"
                        + os.path.basename(multi_files[specie][sample])
                    )
                )
                == True
            ):
                multi_trimmed_files[specie].append(
                    os.path.abspath(
                        "Trimmed_directory/trimmed_"
                        + os.path.basename(multi_files[specie][sample])
                    )
                )
            else:
                trimming_cutadapt(multi_files[specie][sample], threads)
                multi_trimmed_files[specie].append(
                    os.path.abspath(
                        "Trimmed_directory/trimmed_"
                        + os.path.basename(multi_files[specie][sample])
                    )
                )
    return multi_trimmed_files
