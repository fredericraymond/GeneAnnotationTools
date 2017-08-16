from p_file import *
import sys
import docopt

usage = """IPrism For the interpretation and representation of the prism results.

Usage:
  IPrism.py [options] <inputfolder>


Options:
  -h --help     Show this screen
  --version     Version 0.0
  -S FILE1 --table FILE1    prefix of Name of Statistical table(output)(prefix will be use to generate prefix.tsv and prefix.xlsx)
  -c CLSTRFILE --chr CLSTRFILE     cd-hit results file
  -a FILE2 --attr FILE2    prefix of Name of the attributes file. This file is used for cytoscape.(output)(prefix will be use to generate prefix.tsv and prefix.xlsx)
"""

if __name__ == '__main__':
    args = docopt.docopt(doc=usage, argv=sys.argv[1:])
    #print args
    if args['<inputfolder>']:
        all_org = do_it_for_all(args['<inputfolder>'])
        result = compares_organism_considering_clusters_families(all_org)

        if args['--table']:
            display_test(result, args['--table'] + ".xlsx")
            results_in_tsv_file(result, args['--table'] + ".tsv")

        if args['--chr']:
            for org in all_org:
                if org.clusters:
                    print org.f_name
                    org.write_cluster_fasta_files()
                    org.find_ids()
                    org.nodes_attributes(args[
                                             '--chr'])  # )  # "/home/rogia/FR_genome/Data/cdhit/20170707_OrthologComparison_Nares+Modern_Nostoc+GCF_20025-70.clstr"

    no_correspondences = []

    for cle, valeur in nodes.items():
        if len(valeur) == 0:
            no_correspondences.append(cle)

    for cle in no_correspondences:
        del (nodes[cle])

    if args['--attr']:
        display_test(nodes, args['--attr'] +".xlsx")  # 20170707_OrthologComparison_Nares+Modern_Nostoc+GCF_20025-70%_cluster_mapping_Statistiques.xls'
        results_in_tsv_file(nodes, args['--attr'] +".tsv")