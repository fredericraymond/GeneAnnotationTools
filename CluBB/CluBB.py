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
  -p, --boolean   Runs prodigal and cd-hit-2d  if -p is specified in command line.By default, p is set to False.If not specified, we consider that your input folder
                  contains a prodigal pep file or cd-hit results
  -per POURCENTAGE --percent POURCENTAGE  Specify identity percent
  -c, --cdhitr   Runs only cd-hit-2d  if -c is specified. It suppose that the input folder contains prodigal results(pep files)
  
  NB: iF not -p  or not -c specified, we' ll consider that the input folder contains cd-hit-2d results
                

"""

if __name__ == '__main__':
    args = docopt.docopt(doc=usage, argv=sys.argv[1:])
    path = args['<inputfolder>'] + "/*"
    input_files = []
    per = 0.9
    # Load pathways...
    group = load_data('Data/group_by')  # Veuillez a mettre le bon path
    # create big data....
    bgd = collections.defaultdict(dict)
    abS = collections.defaultdict(dict)
    if args['--boolean']:
        pep_files = []
        for index, g in enumerate(glob.glob(path)):
            go = g.split("/")[-1]
            print "prodigal is running on genome " + str(index + 1) + "...."
            f_pep = ".".join(go.split(".")[0:-1]) + ".pep"
            f_nuc = ".".join(go.split(".")[0:-1]) + ".nuc"
            f_out = ".".join(go.split(".")[0:-1]) + ".gff"
            command = "prodigal -i " + g + " -p meta -a " + f_pep + " -d " + f_nuc + " -f gff  -o " + f_out
            os.system(command)
            pep_files.append(f_pep)

        for doc in pep_files:
            f_pep = ".".join(doc.split("/")[-1].split(".")[0:-1])
            print " running cd-hit-2d to compare against brenda' s sequence.... "
            if args['--percent']:
                per = args['--percent']

            cml = "cd-hit-2d -i2 Data/header_modif+NDn_BRENDA_sequences.fasta -i " + doc + " -o " \
                  + f_pep + "  -c " + str(per) + " - n 5 - M 4000 - d 0 "
            os.system(cml)

            input_files.append(f_pep + ".clstr")

    else:
        if args['--cdhitr']:
            if args['--percent']:
                per = args['--percent']
            for doc in glob.glob(path):
                f_pep = ".".join(doc.split("/")[-1].split(".")[0:-1])
                print " running cd-hit-2d to compare against brenda' s sequence.... "
                cml = "cd-hit-2d -i2 Data/header_modif+NDn_BRENDA_sequences.fasta -i " + doc + " -o " \
                      + f_pep + "  -c  " + str(per) + "  -n 5  -M 4000 -d 0"
                os.system(cml)
        else:

            input_files = glob.glob(path)

    for doc in input_files:
        if doc.split("/")[-1].endswith("clstr"):
            f_pep = ".".join(doc.split("/")[-1].split(".")[0:-1])
            print "Parse cd-hit results..."
            cf = find_correspondences(doc)
            ab = abundance_stats(cf)
            abS[f_pep] = ab
            # get stats...
            rs = get_stats(cf, group)
            # print rs.items()
            # add genome information to big data...
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
    add_pathways_length("Pathways.xls", group)
    from_excel_to_tsv("Pathways.xls")
    draw_table(abS, "enzymes.xls")
    results_in_tsv_file(abS, "enzymes.tsv")
