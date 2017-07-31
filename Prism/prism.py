from p_file import *
import sys
# This file generate only stats for genomes tested...

if __name__ == '__main__':
    all_org = do_it_for_all(
        sys.argv[1])
    result = compares_organism_considering_clusters_families(all_org)
    display_test(result, sys.argv[2])

