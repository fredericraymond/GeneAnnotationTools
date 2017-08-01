import collections
import sys
import pickle
import numpy as np
from xlwt import Workbook
import os


def load_data(f_name):
    with open(f_name, 'rb') as f:
        unpickler = pickle.Unpickler(f)
        data = unpickler.load()

    return data


def parse_prodigal_results(f_name):
    genes = []
    for line in open(f_name):
        if line.startswith(">"):
            genes.append(line.split(">")[-1].split("#")[0].strip())

    return genes

# This function parses cd-hit results...


def find_correspondences(f_name,
                         genes):
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
                        code = newline.split(">EC_")[-1].split("|")[0]
                        ls.append(code)
                    newline = f.readline()

                else:
                    if cle:
                        if ls:
                            for i in cle:
                                correspondences[i] = np.unique(ls)
                    line = newline

    return correspondences


# define pathways found considering input gene file
def finding_pat(correspondences, group_by):
    groups_found = collections.defaultdict(list)
    for cle, values in correspondences.items():
        for i in values:
            for pathway, o_enzymes in group_by.items():
                if i in o_enzymes:
                    groups_found[pathway].append((cle, i))

    return groups_found


# Cette fonction cherche sur chacun des contigs des groupes de genes qui appartiennent au meme pathways....


def cluster_finder(groups_found):
    whole_genome_clusters = collections.defaultdict(dict)
    for pathway, orf in groups_found.items():
        clusters = collections.defaultdict(list)
        for indexed, gene in enumerate(orf):
            inp = gene[0]
            contig = "".join(inp.split("_")[0:-1])
            no = inp.split("_")[-1]
            if not no in clusters[contig]:
                clusters[contig].append(no)
        pseudo_clstr = collections.defaultdict(list)

        for cn, nos in clusters.items():
            for i in range(len(nos)):
                current = nos[i]
                tpl = (current,)
                for j in range(i + 1, len(nos)):
                    dif = abs(int(current) - int(nos[j]))
                    if 0 < dif <= 5:
                        tp = (nos[j],)
                        tpl = tpl + tp
                if len(tpl) > 1:
                    pseudo_clstr[cn].append(tpl)
            # Just to remove tuple overlaps...
            indexed = []
            if cn in pseudo_clstr:
                for k in range(len(pseudo_clstr[cn])):
                    current_tuple = pseudo_clstr[cn][k]
                    for m in range(len(current_tuple)):
                        for l in range(k + 1, len(pseudo_clstr[cn])):
                            cmp_tuple = pseudo_clstr[cn][l]
                            if current_tuple[m] in cmp_tuple:
                                if len(current_tuple) >= len(cmp_tuple):
                                    indexed.append(l)
                                else:
                                    indexed.append(k)

                if indexed:
                    b = list(set(indexed))
                    c = range(len(b))
                    for pos, ID in enumerate(b):
                        del pseudo_clstr[cn][ID - c[pos]]
        if len(pseudo_clstr) > 0:
            whole_genome_clusters[pathway] = pseudo_clstr

    return whole_genome_clusters


# Cette fonction evalue en termes de pourcentage tous le pathways retrouves dans le genome ...


def get_stats(correspondences, group_by):
    stats = {}
    ls = []
    results = {}
    for values in correspondences.values():
        ls = list(set(list(values) + ls))
    for sentier, contenu in group_by.items():
        stats[sentier] = 0
        for i in contenu:
            if i in ls:
                stats[sentier] += 1
    for p in group_by.keys():
        percent = (float(stats[p]) / float(len(group_by[p]))) * 100
        if percent > 0:
            results[p] = int(round(percent, 0))

    return results


# Cette fonction evalue en termes de pourcentage tous le pathways retrouves dans le genome par contig...
def get_cluster_density(wg, groups):

    for pathways in wg:
        no = len(groups[pathways])
        for contig in wg[pathways].keys():
            liste_des_clusters = wg[pathways][contig]
            for index, i in enumerate(liste_des_clusters):
                percent = round((float(len(i)) / float(no)) * 100, 0)
                wg[pathways][contig][index] = (i, percent)

# Get brenda' s enzymes found in our genome....


def abundance_stats(correspondences):
    ls = []
    results = {}
    for values in correspondences.values():
        ls = list(set(list(values) + ls))

    for ecc in ls:
        cpt = 0
        for val in correspondences.values():
            if ecc in val:
                cpt += 1

        if cpt > 0:
            results[ecc] = cpt
    return results


# Build excel table...
def draw_table(test_result, f_name):
    book = Workbook()
    feuil1 = book.add_sheet('feuille 1')

    keys_list = []

    for cle in test_result.values():
        keys_list = list(set(cle.keys() + keys_list))

    for fam in keys_list:
        pos = keys_list.index(fam)
        line, col = pos + 1, 0
        feuil1.write(line, col, fam)
        for cle in test_result.keys():
            no_col = test_result.keys().index(cle) + 1
            if fam in test_result[cle]:
                content = test_result[cle][fam]

            else:
                content = 0

            feuil1.write(line, no_col, content)
    for cle in test_result.keys():
        line = feuil1.row(0)
        col = test_result.keys().index(cle) + 1
        line.write(col, cle)

    book.save(f_name)


if __name__ == '__main__':
    # create big data....
    bgd = collections.defaultdict(dict)
    abS = collections.defaultdict(dict)
    # all genomes...
    genome = sys.argv[1:]
    # Load pathways...
    group = load_data('/home/rogia/Documents/SCRIPTS/folder_2/Data/datas/group_by')
    for index, g in enumerate(genome):
        print "prodigal is running on genome " + str(index + 1) + "...."
        f_pep = ".".join(g.split(".")[0:-1]) + ".pep"
        f_nuc = ".".join(g.split(".")[0:-1]) + ".nuc"
        f_out = ".".join(g.split(".")[0:-1]) + ".gff"
        command = "prodigal -i " + g + " -p meta -a " + f_pep + " -d " + f_nuc + " -f gff  -o " + f_out
        os.system(command)
        print " parse prodigal results...."
        # Parse prodigal results...
        gn = parse_prodigal_results(f_pep)
        print " running cd-hit-2d to compare against brenda' s sequence.... "
        f_csltr = f_pep + "_data"
        # Database i2 must be replaced...
        cml = "cdhit-2d -i2 /home/rogia/Documents/SCRIPTS/folder_2/Data/datas/header_modif+BRENDA_90 -i " + f_pep + " -o " \
              + f_csltr + "  -c 0.9 -n 5  -M 1600 -d 0"
        os.system(cml)
        print "Parse cd-hit results..."
        cf = find_correspondences(f_csltr + ".clstr", gn)
        ab = abundance_stats(cf)
        abS[g] = ab
        # get stats...
        rs = get_stats(cf, group)
        print rs.items()
        # add genome information to big data...
        bgd[g] = rs
    print " writing stats file ...."
    draw_table(bgd, "Pathways.xls")
    draw_table(abS, "enzymes.xls")
