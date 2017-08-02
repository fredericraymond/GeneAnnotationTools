import os
import sys
import collections
import numpy as np
from Sanalysis import cluster_finder, get_cluster_density


def parse_mibig_pathways(f_name):
    groups = collections.defaultdict(list)
    with open(f_name, 'r') as f:
        line = f.readline()
        current_entry = line.lstrip(">").rstrip("\n")
        current_cluster_no = current_entry.split("|")[0]
        groups[current_cluster_no].append(current_entry)
        for line in open(f_name):
            if line.startswith(">"):
                entry = line.lstrip(">").rstrip("\n")
                cluster_no = entry.split("|")[0]
                if cluster_no == current_cluster_no:
                    groups[current_cluster_no].append(entry)
                else:
                    current_entry = entry
                    current_cluster_no = current_entry.split("|")[0]
                    groups[current_cluster_no] .append(current_entry)

    return groups


def parse_mibig_cd_hit2d_results(f_name, genes):
    correspondences = collections.defaultdict(list)
    with open(f_name, "r") as f:
        line = f.readline()
        while line:
            if line.startswith(">"):
                ls = []
                cle = []
                newline = f.readline()
                while not newline.startswith(">") and newline:
                    current_name = newline.split(",")[-1].split("...")[0].split(">")[-1]
                    if current_name in genes:
                        cle.append(current_name)
                    else:
                        code = newline.split(">")[-1].split("...")[0]
                        ls.append(code)
                    newline = f.readline()

                else:
                    if cle:
                        if ls:
                            for i in cle:
                                correspondences[i] = np.unique(ls)
                    line = newline

    return correspondences


def parse_blast_results(f_name):
    blast_results = []
    with open(f_name, 'r') as f:
        line = f.readline()
    current_entry = line.split()[0]
    current_score = float(line.split()[-2])
    current_match = line.split()[1]
    for line in open(f_name):
        entry = line.split()[0]
        match = line.split()[1]
        score = float(line.split()[-2])
        if entry == current_entry:
            if score < current_score:
                current_entry = entry
                current_score = score
                current_match = match
        else:
            blast_results.append((current_entry, current_match, current_score))
            current_entry = entry
            current_score = score
            current_match = match
    blast_results.append((current_entry, current_match, current_score))
    return blast_results


# define pathways found considering input gene file
def find_correspondence(blast_r):
    groups_found = collections.defaultdict(list)
    for b in blast_r:
        query = b[0]
        match = b[1]
        gr = match.split("|")[0]
        groups_found[gr].append((query, match))
    return groups_found


def display_cluster_density(wg, f_name):

    contigs = []
    for cluster, values in wg.items():
        contigs = list(set(values.keys() + contigs))

    print contigs
    res = collections.defaultdict(dict)
    for c in contigs:
        res[c] = {}
        for pathways in wg.keys():
            if c in wg[pathways].keys():
                res[c][pathways] = []
                for x in wg[pathways][c]:
                    res[c][pathways].append(x)

    with open(f_name, "w") as f:
        f.write("\t".join(["Contigs ", " Pathways", "density"]))
        f.write("\n")
        for c, values in res.items():
            f.write(c)
            f.write("\n")
            f.write("=================================================================")
            f.write("\n")
            for p, liste in values.items():
                for clstr in liste:
                    mem1 = clstr[0]
                    mem2 = clstr[-1]
                    line = "\t".join([str(mem1), str(p), str(mem2)])
                    f.write(line)
                    f.write("\n")


if __name__ == '__main__':
    # all genomes...
    genome = sys.argv[1:]
    # Load pathways...
    group = parse_mibig_pathways("/home/rogia/Documents/MIBIG/MIBiG_prot_seqs_1.3.fasta")
    for index, g in enumerate(genome):
        print "prodigal is running on genome " + str(index + 1) + "...."
        f_pep = ".".join(g.split(".")[0:-1]) + ".pep"
        f_nuc = ".".join(g.split(".")[0:-1]) + ".nuc"
        f_out = ".".join(g.split(".")[0:-1]) + ".gff"
        command = "prodigal -i " + g + " -p meta -a " + f_pep + " -d " + f_nuc + " -f gff  -o " + f_out
        os.system(command)
        f_br = ".".join(g.split(".")[0:-1]) + ".blt"
        print "blat is running on genome " + str(index + 1) + "...."
        command = "blat /home/rogia/Documents/MIBIG/MIBiG_prot_seqs_1.3.fasta " + f_pep + " " + f_br +" -prot -out=blast8"
        os.system(command)
        print "parse blat results..."
        br = parse_blast_results(f_br)
        print "Finding pathways..."
        gf = find_correspondence(br)
        print "Finding clusters..."
        cf = cluster_finder(gf)
        get_cluster_density(cf, group)
        #
        f_clstr = ".".join(g.split(".")[0:-1]) + ".clstr"
        display_cluster_density(cf, f_clstr)
