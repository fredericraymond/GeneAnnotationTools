import sys
import csv
import collections
from xlwt import Workbook
import glob


def load_cognames2003_2004_csv(f_name):
    cr = csv.reader(open(f_name, "rb"))
    kogs = collections.defaultdict(list)
    for index, row in enumerate(cr):
        if len(row) >= 7:
            cog_id = row[6]
            protein_id = row[0]
            kogs[cog_id].append(protein_id)

    return kogs


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


def find_correspondence(kogs, blast_r):
    data = []
    for b in blast_r:
        query = b[0]
        match = b[1]
        for cog_id, protein_id in kogs.items():
            if match.split("|")[1] in protein_id:
                data.append((query, match, str(b[2]), cog_id))

    return data


def get_cog_count(data):
    cog_list = {}
    for c in data:
        cog = c[-1]
        if cog in cog_list:
            cog_list[cog] += 1
        else:
            cog_list[cog] = 1
    return cog_list


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
                cpt .append(str(content))
            f.write("\t".join(cpt))
            f.write("\n")


if __name__ == '__main__':
    lars = collections.defaultdict(dict)
    # Load cogs ...
    tab = load_cognames2003_2004_csv("cog2003-2014.csv")
    for b in glob.glob(sys.argv[1] + "/*"):
        print "parse blast results...."
        blast_r = parse_blast_results(b)
        # For each blast entry, find the corresponding kogg...
        print "Finding correspondences..."
        df = find_correspondence(tab, blast_r)
        # edit result file...
        print "editing results...."
        f_ile = str(b.split("/")[-1]) + ".corr"
        with open(f_ile, "w") as f:
            f.write("\t".join(["QUERY", "MATCH", "SCORE", "COG\n"]))
            for line in df:
                l = "\t".join([str(line[0]), str(line[1]), str(line[2]), str(line[3])])
                f.write(l)
                f.write("\n")
        # Count...
        print "Counting..."
        c_g = get_cog_count(df)
        # add to dictionary...
        lars[str(b.split("/")[-1])] = c_g
    print "draw table...."
    draw_table(lars, "Cog_Count.xls")
    results_in_tsv_file(lars, "cog_count.tsv")
