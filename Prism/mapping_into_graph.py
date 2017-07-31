from p_file import *
import sys
import time

# annotate cd-hit clusters...to generate graph...
if __name__ == '__main__':
    all_org = do_it_for_all(
        sys.argv[1])
    start_time = time.time()
    for org in all_org:
        if org.clusters:
            print org.f_name
            org.find_gene_prodigal_identifier()
            org.nodes_attributes(sys.argv[2]
                                 )  # "/home/rogia/FR_genome/Data/cdhit/20170707_OrthologComparison_Nares+Modern_Nostoc+GCF_20025-70.clstr"
    print time.time() - start_time

    empty_keys = []

    for cle, valeur in nodes.items():
        if len(valeur) == 0:
            empty_keys.append(cle)

    for cle in empty_keys:
        del (nodes[cle])

    display_test(nodes,
                 sys.argv[
                     3])  # '20170707_OrthologComparison_Nares+Modern_Nostoc+GCF_20025-70%_cluster_mapping_Statistiques.xls'
