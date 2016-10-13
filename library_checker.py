import sys, os,csv, getopt
import csv
import ntpath

#use : python library_checker.py <library folder to check>

def main(argv):


    lib_folder=sys.argv[1]
    for ffile in os.listdir(lib_folder):
        if not ffile.startswith("."):
            libchecker(os.path.join(lib_folder,ffile))



def libchecker(library):
    name_dict=dict()
#    IPython.embed()
    with open(library, "r") as tsv:
        for line in csv.reader(tsv, dialect="excel-tab"):
            if line[0].strip() in name_dict:

                name_dict[line[0].strip()]+=1
            else:
                name_dict[line[0].strip()]=1
    tsv.close()
    with open(library, "r") as tsv:
        lines=tsv.readlines()
    tsv.close()
    with open(library, "w") as tsv:
        for line in lines:
            if name_dict[line.split()[0]]==10:
                tsv.writelines(line)
    tsv.close()

if __name__ == "__main__":

    main(sys.argv[1:])
