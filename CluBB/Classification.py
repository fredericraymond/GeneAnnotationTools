import pickle


def load_data(f_name):
    with open(f_name, 'rb') as f:
        unpickler = pickle.Unpickler(f)
        data = unpickler.load()

    return data


def visualize_brenda_pathways(pathways, ec_number):
    line = " EC_NUMBER" + "\t" + "\t".join(pathways.keys()) + "\n"
    with open("visualize_pathways.tsv", "w") as f:
        f.write(line)
        for ec_code in ec_number.keys():
            summarize = []
            for x in pathways.keys():
                if ec_code.rstrip("\n") in pathways[x]:
                    boolean = 1
                else:
                    boolean = 0
                summarize.append(str(boolean))
            line = ec_code + "\t" + "\t".join(summarize) + "\n"
            if str("1") in summarize:
                f.write(line)


def get_similar_pathways(pathways):

    with open("pathways_with_same_enzymes.txt", "w") as f:
        for x in pathways:
            corr = []
            poi = [y for y in pathways if y != x]
            for y in poi:
                if set(pathways[x]) == set(pathways[y]):
                    corr.append(y)
                    poi.remove(y)
            if corr:
                corr.append(x)
                f.write(";".join(corr) + "\n")

if __name__ == '__main__':
    ecc = load_data("Data/enzymes")
    path = load_data("Data/group_by")
    #get_similar_pathways(path)
    visualize_brenda_pathways(path, ecc)
