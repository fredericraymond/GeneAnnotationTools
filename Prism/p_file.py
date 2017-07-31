import json
import numpy as np
import glob
from xlwt import Workbook
import os
import os.path
from lxml import html
import shutil
import collections
from Bio import SeqIO
import re

nodes = collections.defaultdict(str)


class Cluster:
    def __init__(self):
        self.family = []
        self.orfs = []
        self.frame = "Unfound"
        self.start = 0
        self.end = 0
        self.contig_name = []
        self.homolog_clusters = []
        self.type = []


class Orf:
    def __init__(self):
        self.name = "Unfound"
        self.frame = "Unfound"
        self.stop = 0
        self.start = 0
        self.domains = []


class Motif:
    def __init__(self, start, stop, score, name):
        self.start = start
        self.stop = stop
        self.score = score
        self.name = name


class Organism:
    def __init__(self, name, path):
        self.f_name = name
        self.clusters = []
        self.path = path

    def display_organism(self):

        print self.f_name
        print "\t".join(["Results : ", "Clusters", str(len(self.clusters))])
        print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        for cluster in self.clusters:
            print "\t".join(["cluster", str(self.clusters.index(cluster) + 1)])
            print "------------"
            print "Family : ",
            for name in cluster.family:
                print name,
            print "\n"

            print "Orfs in cluster : " + "\n"
            print "\t".join(["orfs", "start", "stop", "frame", "domains"])
            for orf in cluster.orfs:
                print "\t".join([orf.name, str(orf.start), "", str(orf.stop), orf.frame, " "]),
                for domain in orf.domains:
                    print domain.name,
                print "\n"
            print "\n"

    # For each cluster and for each gene, we'll identify it's sequence.....
    def find_gene_sequence(self, no_cluster):

        sequences_by_orf_by_cluster = {}
        orf_names_by_cluster = []
        for clusters in self.clusters:

            no = self.clusters.index(clusters)
            if no_cluster == no + 1:
                filename = "/".join([self.path, "cluster_" + str(no + 1) + ".html"])

                for orf in self.clusters[no].orfs:
                    orf_names_by_cluster.append(orf.name)

                with open(filename, "r") as f:
                    page = f.read()
                    tree = html.fromstring(page)
                    sequence = tree.find_class("sequence")
                    for i in sequence:
                        name_orf = orf_names_by_cluster[sequence.index(i)]

                        sequences_by_orf_by_cluster[name_orf] = i.text_content().lstrip('\n')

        return sequences_by_orf_by_cluster

    # We'll find (using prodigal file) each gene considering it's sequence.....
    def find_gene_prodigal_identifier(self):

        peptides = list(SeqIO.parse(self.path + "/" + self.f_name, "fasta"))
        for index, item in enumerate(self.clusters):
            dict = self.find_gene_sequence(index + 1)
            for orf, pep in dict.items():
                subject = pep
                for pos, record in enumerate(peptides):
                    pattern = str(record.seq)
                    match = re.match(pattern, subject)
                    if match is not None:
                        if match.end() == len(subject):
                            item.contig_name.append(str(record.name))
                            break

                        elif match.end() == len(pattern) or pattern == subject[0:match.end() + 1]:
                            for i in range(1, len(peptides) - pos + 1):
                                pattern = "".join([str(p.seq) for p in peptides[pos:pos + i]])
                                if pattern[-1] == "*":
                                    pattern = pattern[0:-1]

                                if pattern == pep:
                                    for x in peptides[pos:pos + i]:
                                        item.contig_name.append(x.name)
                                    break

    # This function reads cd-hit results file and parse it....
    # Only cd-hit clusters with a gene detected by prism will be annotated ....

    def nodes_attributes(self, f1_name):

        with open(f1_name, 'r') as f:
            line = f.readline()
            while len(line.split()) != 0 and line[0] == ">":
                node = line.split(">")[-1].rstrip('\n')
                if node not in nodes:
                    nodes[node] = {}
                ss_line = f.readline()
                while len(ss_line.split()) != 0 and ss_line[0] != ">":
                    sample_gene = ss_line.split(">")[1].split("...")[0]
                    if self.f_name == sample_gene.split("|")[0]:
                        gene = sample_gene.split("|")[-1]
                        for clusters in self.clusters:
                            if gene in clusters.contig_name:
                                fam = "/".join(clusters.type)
                                if fam not in nodes[node]:
                                    nodes[node][fam] = 1
                                else:
                                    nodes[node][fam] += 1
                    ss_line = f.readline()

                else:
                    line = ss_line

        return nodes

   # We move all prism prodigal results to a new diretory....
   # These files will be used to run cd-hit...
  # This is useful to be sure that we'll we work with the same gene whatever the cd-hit clusters won' t be correctly anotated


def move_all_progigal_output(folder_path):
    os.makedirs("/home/rogia/FR_genome/Data/prodigal")
    for path, dirs, files in os.walk(folder_path):
        for filename in files:
            if filename.endswith('.fasta.json') or filename.endswith('.fna.json'):
                for file in glob.glob("/".join([path, 'orfs_protein.fasta'])):
                    new_name = os.path.join(os.path.sep, path, ".".join(filename.split(".")[0:-1]))
                    os.rename(file, new_name)
                    shutil.copy2(new_name, os.path.join(os.sep, os.getcwd(), '/home/rogia/FR_genome/Data/prodigal'))

# The folder path must contains all prism results directory...
# Useful to parse more than one genome prism reults...


def do_it_for_all(folder_path):
    result = []
    f_names = []

    for path, dirs, files in os.walk(folder_path):
        for filename in files:
            if filename.endswith('.json') and 'cluster_' not in filename:
                where = "/".join([path, filename])
                f_names.append(where)

    for files in f_names:
        org = read_prism_file(files)
        result.append(org)
    return result

# This function parse prism results file....


def read_prism_file(f_name):
    # Load prism file
    obj = json.loads(open(f_name).read())
    path = "/".join(f_name.split("/")[0:-1])
    # Preparing object i will use
    # 1..
    cluster = Cluster()
    # 2..
    name = obj["prism_results"]["input"]["filename"]
    org = Organism(name, path)
    # 3..
    gene = Orf()
    # 2 Cases...
    # no clusters detected...
    if not obj["prism_results"]["clusters"]:
        return org
    else:
        # iterates ...
        # Step 1
        for clusters in obj["prism_results"]["clusters"]:
            # Seting cluster object parameters
            cluster.end = clusters["end"]
            cluster.family = np.unique(clusters["family"])
            cluster.frame = clusters["frame"]
            cluster.start = clusters["start"]
            cluster.type = clusters["type"]
            cluster.homolog_clusters = clusters["homolog_clusters"]
            # Step 2
            # Fill orf's list
            for orf in clusters["orfs"]:
                gene.start = orf["start"]
                gene.stop = orf["stop"]
                gene.frame = orf["frame"]
                # gene.sequence = orf["sequence"]
                gene.name = orf["name"]

                # Finding domains...
                for domain in orf["domains"]:
                    motif = Motif(domain["start"], domain["stop"], domain["score"], domain["full_name"])
                    gene.domains.append(motif)

                    # Add orf to cluster orfs
                cluster.orfs.append(gene)
                # Reinitialize...
                gene = Orf()
            # Add cluster to organism clusters
            org.clusters.append(cluster)
            # Reinitialize...
            cluster = Cluster()

    return org

# Count for each genome how many clusters of each type...


def compares_organism_considering_clusters_families(result):
    test_result = {}
    for genome in result:
        test_result[genome.f_name] = {}
        for i in range(len(genome.clusters)):
            fam = "/".join(genome.clusters[i].type)
            if fam not in test_result[genome.f_name]:
                test_result[genome.f_name][fam] = 1
            else:
                test_result[genome.f_name][fam] += 1

    return test_result

# Produces excel table wih all genome and their stats...


def display_test(test_result, f_name):
    book = Workbook()
    feuil1 = book.add_sheet('feuille 1')

    keys_list = []

    for cle in test_result.values():
        keys_list = list(set(cle.keys() + keys_list))

    for fam in keys_list:
        pos = keys_list.index(fam)
        ligne, col = 0, pos + 1
        feuil1.write(ligne, col, fam)
        for cle in test_result.keys():
            no_ligne = test_result.keys().index(cle) + 1
            ligne = feuil1.row(no_ligne)
            if fam in test_result[cle]:
                contenu = test_result[cle][fam]

            else:
                contenu = 0

            ligne.write(col, contenu)
    for cle in test_result.keys():
        ligne, col = test_result.keys().index(cle) + 1, 0
        feuil1.write(ligne, col, cle)

    book.save(f_name)
