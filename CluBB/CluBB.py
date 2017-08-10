import docopt
from Sanalysis import *
import glob
import sys
import os

usage = """CluBB Clustering based on brenda's enzymes.

Usage:
  CluBB.py [options] <inputfolder>


Options:
  -h --help     Show this screen
  --version     Version 0.0
  -p, --boolean   Runs prodigal if -p is specified in command line.By default, p is set to False

"""

if __name__ == '__main__':
    args = docopt.docopt(doc=usage, argv=sys.argv[1:])
    path = args['<inputfolder>'] + "/*"
    input_files = []
    # Load pathways...
    group = load_data('group_by')
    # create big data....
    bgd = collections.defaultdict(dict)
    abS = collections.defaultdict(dict)
    if args['--boolean']:

        for index, g in enumerate(glob.glob(path)):
            go = g.split("/")[-1]
            print "prodigal is running on genome " + str(index + 1) + "...."
            f_pep = ".".join(go.split(".")[0:-1]) + ".pep"
            f_nuc = ".".join(go.split(".")[0:-1]) + ".nuc"
            f_out = ".".join(go.split(".")[0:-1]) + ".gff"
            command = "prodigal -i " + g + " -p meta -a " + f_pep + " -d " + f_nuc + " -f gff  -o " + f_out
            os.system(command)
            input_files.append(f_pep)

    else:
        input_files = glob.glob(path)

    for doc in input_files:
        f_pep = ".".join(doc.split("/")[-1].split(".")[0:-1])
        print " parse prodigal results...."
        # Parse prodigal results...
        gn = parse_prodigal_results(doc)
        print " running cd-hit-2d to compare against brenda' s sequence.... "
        f_csltr = f_pep + "_data"
        cml = "cd-hit-2d -i2 header_modif+NDn_BRENDA_sequences.fasta -i " + doc + " -o " \
           + f_csltr + "  -c 0.9 -n 5  -M 4000 -d 0"
        os.system(cml)
        print "Parse cd-hit results..."
        cf = find_correspondences(f_csltr + ".clstr", gn)
        ab = abundance_stats(cf)
        abS[f_pep] = ab
        # get stats...
        rs = get_stats(cf, group)
        # print rs.items()
        #add genome information to big data...
        bgd[f_pep] = rs
        # Finding pathways....
        pf = finding_pat(cf, group)
        print "Finding clusters..."
        wg = cluster_finder(pf)
        # get clusters density...
        f_clstr = f_pep + ".clstr"
        display_cluster_density(wg, f_clstr, group)
    print " writing stats file ...."
    draw_table(bgd, "Pathways.xls")
    draw_table(abS, "enzymes.xls")


