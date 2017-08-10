import urllib2
from bs4 import BeautifulSoup
import collections
import pickle

# to open the link...


def open_link(link):
    try:
        page = urllib2.urlopen(link)
        soup = BeautifulSoup(page, "lxml")
    except:
        soup = None

    return soup

# to get all brenda's enzymes and their common names...


def get_enzymes_from_brenda():
    enzymes_common_names = {}
    wiki = "http://www.brenda-enzymes.org/all_enzymes.php"
    soup = open_link(wiki)
    # assert soup.title.string == 'BRENDA - All enzymes'
    table = soup.find('div', class_="equal")
    # print table
    for row in table.findAll('div', {'class': 'row'}):
        line = row(text=True)  # just keep the index 0, 1 ecc_id common_name
        ecc_code = line[0]
        enzyme_common_name = line[1]
        enzymes_common_names[str(ecc_code)] = str(enzyme_common_name)

    return enzymes_common_names

# this returns brenda's pathways and link associated...


def get_brenda_pathways():
    # Just make list of metabolic link
    pathways = collections.defaultdict(list)
    wiki = " http://www.brenda-enzymes.org/search_result.php?a=137"
    soup = open_link(wiki)
    tab = soup.findAll('a')
    for a in tab:
        if 'href' in a.attrs:
            link = a['href']
            if 'Search=Search' in str(link):
                pathway_name = a(text=True)[0]
                print pathway_name
                pathways[str(pathway_name)].append(str(link))

    return pathways

# this returns for each enzyme, the set of pathways ...


def get_enzymes_pathways():
    enzymes_pathways = collections.defaultdict(list)
    wiki = 'http://www.brenda-enzymes.org/result_download.php?a=137&RN=&RNV=1&os=&pt=&FNV=&tt=' \
           '&SYN=&Textmining=&T[0]=2&W[1]=&T[1]=2&nolimit=1'
    soup = open_link(wiki)

    for line in soup.getText().split("\n"):
        if len(line.split("\t")) >= 4:
            enzyme = line.split("\t")[0]
            pathway = line.split("\t")[-2]
            enzymes_pathways[str(enzyme)].append(str(pathway))

    return enzymes_pathways

# this returns for each pathway, the set of enzymes...


def sort_by_metabolic_pathways(b_pathways, b_enzymes_pathways):
    group_by_pathways = collections.defaultdict(list)

    for pathway in b_pathways.keys():
        for ec, ec_pathways in b_enzymes_pathways.items():
            if pathway in ec_pathways:
                group_by_pathways[pathway].append(ec)

    return group_by_pathways

# i download all sequences from brenda's database because they ' re well identified...
# i have noticed that brenda' s database is redundant...


def get_sequences_from_brenda(enzymes_common_names):

    with open("BRENDA_sequences.fasta", "w") as f:
        for ecc in enzymes_common_names.keys():
            print ecc, enzymes_common_names.keys().index(ecc)
            wiki = 'http://www.brenda-enzymes.org/sequences.php?download=allfasta&ec=' + str(ecc)
            soup = open_link(wiki)
            if soup is not None and len(soup) == 1:
                for line in soup.getText().split("\r\n"):
                    if "|" in line:
                        f.write(">" + line)

                    else:
                        f.write(line)
                    f.write("\n")
            f.write("\n")


# This code proves that there's 200 000 repeated sequences in brenda's database....


def eliminate_duplicated_lines(f_name):
    memory = []
    f_output = "NDn_" + f_name
    fiche = open(f_output, 'w')
    with open(f_name, 'r') as f:
        line = f.readline()
        while line:
            if line.startswith(">"):
                up = line.split("|")[0].split(">")[-1]
                if up not in memory:
                    print up
                    memory.append(up)
                    fiche.write(line)
                    line_suv = f.readline()
                    while not line_suv.startswith(">"):
                        fiche.write(line_suv)
                        line_suv = f.readline()
                    else:
                        line = line_suv
                else:
                    line = f.readline()
            else:
                line = f.readline()
    fiche.close()

# This function modify headers to facilitate cd-hit results parsing...


def adjust_rp_file_header(f_name):
        fiche = open("header_modif+" + f_name, 'w')
        for line in open(f_name):
            if line.startswith(">"):
                ec_code = line.split("|")[2]
                b1 = line.lstrip(">").split(ec_code)[0]
                b2 = line.split("|")[-1]
                l = ">" + "_".join(ec_code.split()) + "|" + b1 + b2
                fiche.write(l)
            else:
                fiche.write(line)


# the aims of this script is to build a database, so we' ll saved all data' s...


def save_objects(f_name, data):
    with open(f_name, 'wb') as f:
        pickler = pickle.Pickler(f)
        pickler.dump(data)


if __name__ == '__main__':
    # Get enzymes and their common names...
    print "Get enzymes and their common names..."
    enzymes = get_enzymes_from_brenda()
    for cle, values in enzymes.items():
        print cle, values
    # Load enzymes sequences from BRENDA database...
    print "Load enzymes sequences from BRENDA database..."
    get_sequences_from_brenda(enzymes)
    # Get brenda' s pathways
    print "Get brenda' s pathways"
    metabo_pathway = get_brenda_pathways()
    for cle, values in metabo_pathway.items():
        print cle, values
    # Get enzymes and their pathways...
    print "Get enzymes and their pathways..."
    path = get_enzymes_pathways()
    for cle, values in path.items():
        print cle, values
    # Sort enzymes_pathways by metabolic pathways...
    group_by = sort_by_metabolic_pathways(metabo_pathway, path)
    # Let' save our data's...
    save_objects('enzymes', enzymes)
    print "enzymes saved...."
    save_objects('pathways', metabo_pathway)
    print "brenda pathways saved..."
    save_objects('enzymes+pathways', path)
    print "enzymes and pathways saved..."
    save_objects('group_by', group_by)
    print "group_by metabolic pathways saved..."
    # Let' s do a little modifications..
    eliminate_duplicated_lines("BRENDA_sequences.fasta")
    adjust_rp_file_header("NDn_BRENDA_sequences.fasta")