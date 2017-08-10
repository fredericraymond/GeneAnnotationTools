import sys


def divide(f_name):
    cpt = 0
    base = f_name + "_1"
    fiche = open(base, "w")
    for line in open(f_name):
        if line.startswith(">"):
            cpt += 1
        if cpt <= 1000:
            fiche.write(line)

        else:
            cpt = 0
            fiche.close()
            fiche = open(f_name + "_" + str(int(base.split("_")[-1]) + 1), "w")
            base = f_name + "_" + str(int(base.split("_")[-1]) + 1)

if __name__ == '__main__':

    divide(sys.argv[1])