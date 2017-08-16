import collections
import pickle
import numpy as np
from xlwt import Workbook
import pandas as pd
import os
import xlrd


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
        # Etape 1: On identifie tous les genes du pathways qui appartiennent au meme contigs
        for indexed, gene in enumerate(orf):
            inp = gene[0]
            corr = gene[-1]
            contig = "_".join(inp.split("_")[0:-1])
            no = inp.split("_")[-1]
            if not no in clusters[contig]:
                clusters[contig].append(no)
        pseudo_clstr = collections.defaultdict(list)
        # Etape 2: On retrouve parmi les genes du contig ceux qui sont contigus
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
            # on calcule la densite du cluster...
            for contig, valeurs in pseudo_clstr.items():
                for pos, bcg in enumerate(valeurs):
                    cpt = []
                    for i in bcg:
                        gene = contig + "_" + str(i)
                        for gr in orf:
                            if gene in gr:
                                for corr in gr[1:]:
                                    cpt.append(corr)
                    cpt = np.unique(cpt)
                    pseudo_clstr[contig][pos] = (pseudo_clstr[contig][pos], len(cpt), cpt)

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

        if stats[p] > 0:
            results[p] = stats[p] #str(len(group_by[p]))  #int(round(percent, 0))

    return results


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


def display_cluster_density(wg, f_name, group_by):

    contigs = []
    for cluster, values in wg.items():
        contigs = list(set(values.keys() + contigs))

    res = collections.defaultdict(dict)
    for c in contigs:
        res[c] = {}
        for pathways in wg.keys():
            if c in wg[pathways].keys():
                res[c][pathways] = []
                for x in wg[pathways][c]:
                    res[c][pathways].append(x)

    with open(f_name, "w") as f:
        f.write("\t".join(["Contigs ", "genes", "Pathways", "Pathways_length ", "cluster density", "Corresponding EC"]))
        f.write("\n")
        for c, values in res.items():
            for p, liste in values.items():
                for clstr in liste:
                    series = []
                    for i in clstr[0]:
                        series.append(str(i))
                    mem1 = ",".join(series)
                    mem2 = clstr[1]
                    mem3 = ",".join(clstr[-1])
                    line = "\t".join([str(c), str(mem1), str(p), str(len(group_by[p])), str(mem2), str(mem3)])
                    f.write(line)
                    f.write("\n")


def add_pathways_length(f_name, group_by):
    data = []
    wb = xlrd.open_workbook(f_name)
    sh = wb.sheet_by_index(0)
    for rownum in range(sh.nrows):
        data.append(sh.row_values(rownum)[0])

    data = data[1:]
    no = []
    df = pd.read_excel(f_name)
    for pathways in data:
        length = len(group_by[pathways])
        no.append(length)

    df.insert(0, "Pathways length", no)
    os.system("rm " + f_name)
    writer = pd.ExcelWriter(f_name)
    df.to_excel(writer, 'Sheet1')
    writer.save()


def results_in_tsv_file(test_result, f_name):
    with open(f_name, "w") as f:
        sample = test_result.keys()
        sample.insert(0, "FEATURES")
        print sample
        f.write("\t".join(sample))
        f.write("\n")
        keys_list = []
        for cle in test_result.values():
            keys_list = list(set(cle.keys() + keys_list))
        for fea in keys_list:
            cpt = [fea]
            for cle in test_result.keys():
                if fea in test_result[cle]:
                    content = test_result[cle][fea]

                else:
                    content = 0
                cpt.append(str(content))
            f.write("\t".join(cpt))
            f.write("\n")