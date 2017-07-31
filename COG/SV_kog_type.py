import sys
import csv
import collections


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
    f_name = open("annotation.txt", 'w')
    f_name.write("\t".join(["Query", "Best COG match", "E score", "Cog Code"]))
    for b in blast_r:
        query = b[0]
        match = b[1]
        for cog_id, protein_id in kogs.items():
            if match.split("|")[1] in protein_id:
                f_name.write("\n" + "\t".join([query, match, str(b[2]), cog_id]))

    f_name.close()


if __name__ == '__main__':

    # Load cogs ...
    tab = load_cognames2003_2004_csv("cog2003-2014.csv")
    # Load and parse blast results...
    blast_r = parse_blast_results(sys.argv[0])

    # For each blast entry, find the corresponding kogg...
    find_correspondence(tab, blast_r)

