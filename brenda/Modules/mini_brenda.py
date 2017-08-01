import os


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


if __name__ == '__main__':
    # Run cd-hit on  brenda's sequences....
    cmd = "cdhit -i NDn_BRENDA_sequences.fasta -o BRENDA_90 -c 0.90 -n 5 -d 0 -M 3728"
    os.system(cmd)
    # modify headers......
    adjust_rp_file_header("BRENDA_90")