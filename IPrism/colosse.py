import os
import sys

# To generate sh files for colosse...


def create_sh_file(source, destination):
    f_names = []
    listed = []
    for filename in os.listdir(source):
        path = "/".join([source, filename])
        f_names.append(path)

    nb_fois = len(f_names) // 8
    if len(f_names) - 8 * nb_fois != 0:
        nb_fois += 1

    for j in range(nb_fois):
        script_name = "script_" + str(j + 1) + ".sh"
        f_file = open(destination + "/" + script_name, 'w')
        f_file.write("\n".join(
            ["#!/bin/bash", "#PBS -N " + script_name, "#PBS -o " + script_name + ".stdout",
             "#PBS -e " + script_name + ".stderr",
             "#PBS -A nne-790-aa", "#PBS -l walltime=24:00:00", "#PBS -l nodes=1:ppn=8"
                , "#PBS -q default", "\n", "module load nne-790-ab/prodigal/2.6", "module load compilers/java/1.8",
             "module load compilers/intel/14.0", "module load apps/hmmer/3.1b2", "", "\n"]))

        if len(f_names) >= 8:
            for i in range(8):
                listed.append(f_names[i])
        else:
            listed = f_names

        for elt in listed:
            file = elt.split("/")[-1]
            grid_name = str(".".join(file.split(".")[0:-1])) + file[-1]
            command = "java -jar /rap/nne-790-ab/software2/prism-releases-master/prism.jar -a -p -f " + str(
                elt) + " -grid " + grid_name + " -reg -res -rib -sug -tt -w 10000 -r /rap/" \
                                               "nne-790-ab/software2/prism-releases-master/" \
                                               "prism/WebContent" + "& "
            if listed.index(elt) == len(listed) - 1:
                command = "java -jar /rap/nne-790-ab/software2/prism-releases-master/prism.jar -a -p -f " + str(
                    elt) + " -grid " + grid_name + " -reg -res -rib -sug -tt -w 10000 -r /" \
                                                   "rap/nne-790-ab/software2/prism-releases" \
                                                   "-master/prism/WebContent "
            f_file.write(command + "\n")

        f_file.write("wait")
        f_names = list(set(f_names) - set(listed))
        listed = []

        f_file.close()


if __name__ == '__main__':
    create_sh_file(sys.argv[1], sys.argv[2])
