import xlrd
import sys


def from_excel_to_tsv(f_name):
    output_file = ".".join(f_name.split("/")[-1].split(".")[0:-1]) + ".tsv"
    # ouverture du fichier Excel
    wb = xlrd.open_workbook(f_name)
    with open(output_file, "w") as f:
        for index in range(len(wb.sheet_names())):
            sheet = wb.sheet_by_index(index)
            for rownum in range(sheet.nrows):
                row = sheet.row_values(rownum)
                line = []
                for val in row:
                    line.append(str(val))

                f.write("\t".join(line))
                f.write("\n")

if __name__ == '__main__':
    from_excel_to_tsv(sys.argv[1])