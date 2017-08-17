import numpy as np
import csv
import collections
import glob
import sys


def load_cognames2003_2004_csv(f_name):
    cr = csv.reader(open(f_name, "rb"))
    kogs = collections.defaultdict(list)
    for index, row in enumerate(cr):
        if len(row) >= 7:
            cog_id = row[6]
            protein_id = row[0]
            kogs[cog_id].append(protein_id)

    return kogs


def parse_cd_hit2d_results(f_name):
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
                    if not current_name.startswith("gi"):
                        cle.append(current_name)
                    else:
                        code = newline.split(">")[-1].strip("\n")
                        ls.append(code)
                    newline = f.readline()

                else:
                    if cle:
                        if ls:
                            for i in cle:
                                correspondences[i] = np.unique(ls)
                    line = newline

    return correspondences


def find_correspondence(kogs, correspondences):
    data = []
    for gene, values in correspondences.items():
        corr = {}
        for match in values:
            id = match.split("...")[0].split("|")[1]
            per = match.split("...")[-1]
            for cog_id, protein_id in kogs.items():
                if id in protein_id:
                    if cog_id not in corr:
                        corr[cog_id] = (per, match.split("...")[0])
                    else:
                        if corr[cog_id][0] < per:
                            match = corr[cog_id][1]
                            corr[cog_id] = (per, match)

        for cog in corr:
            data.append((gene, corr[cog][-1], corr[cog][0], cog))

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


def results_in_tsv_file(test_result, f_name):
    with open(f_name, "w") as f:
        sample = test_result.keys()
        sample.insert(0, "FEATURES")
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


if __name__ == '__main__':
    lars = collections.defaultdict(dict)
    # Load cogs ...
    tab = load_cognames2003_2004_csv("cog2003-2014.csv")
    for b in glob.glob(sys.argv[1] + "/*"):
        if b.split("/")[-1].endswith("clstr"):
            print "parse cd-hit-2d results..."
            corr = parse_cd_hit2d_results(b)
            print corr
            print "find correspondences...."
            df = find_correspondence(tab, corr)
            print "editing results...."
            f_ile = str(".".join(b.split("/")[-1].split(".")[0:-1])) + ".corr"
            with open(f_ile, "w") as f:
                f.write("\t".join(["QUERY", "MATCH", "SIMILARITE", "COG\n"]))
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
    results_in_tsv_file(lars, "cog_count.tsv")
