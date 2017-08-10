import os
import glob
import sys


def all_pep_files(folder_path):
    args = []
    for path, dirs, files in os.walk(folder_path):
        if 'orfs_protein.fasta' in os.listdir(path):
            if 'cluster_1.json' in os.listdir(path):
                for file in glob.glob("/".join([path, 'orfs_protein.fasta'])):
                    header = path.split("/")[-1].split(".")[0]
                    modif_pep_file = "/".join([path, 'Orfs_protein.fasta'])
                    cmd = "sed s/\>/\>" + str(header) + "\|/ " + file + '> ' + str(modif_pep_file)
                    os.system(cmd)
                    args.append(modif_pep_file)

    # Concatenez les fichiers
    cmd = "cat  " + "\t".join(args) + "> cat_prodigal.fasta"
    os.system(cmd)

if __name__ == '__main__':
    all_pep_files(sys.argv[1])