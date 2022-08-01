import os


def writing_script(multi_reads, threads, ref_prot, name_sample, slope, actual_path):
    os.chdir(actual_path)
    os.mkdir("input_dir_pegasus")
    os.mkdir("output_dir_pegasus")
    pegasus_script = open("pegasus_script.txt", "w")
    for specie in range(0, len(multi_reads)):
        input_file = open("input_dir_pegasus/input_" + str(specie) + ".txt", "w")
        print(" ".join(multi_reads[specie]) + "\n")
        input_file.write(" ".join(multi_reads[specie]) + "\n")
        input_file.close()
        name_input_file = open(
            "input_dir_pegasus/name_input_" + str(specie) + ".txt", "w"
        )
        name_input_file.write(name_sample[specie])
        name_input_file.close()
        pegasus_script.write(
            "TASK\tinput_"
            + str(specie)
            + " -c "
            + str(threads)
            + '\tbash -c " python3 '
            + os.path.dirname(__file__)
            + "/main.py "
            + " --multi_files "
            + os.path.abspath("input_dir_pegasus/input_" + str(specie) + ".txt")
            + " --ref_p "
            + str(ref_prot)
            + " -t "
            + str(threads)
            + " --slope "
            + str(slope)
            + " --output output_dir_pegasus/output_"
            + str(specie)
            + '" '
            + "--name_samples name_input_"
            + str(specie)
            + ".txt \n"
        )
    pegasus_script.close()
