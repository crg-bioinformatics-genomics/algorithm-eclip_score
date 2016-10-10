import numpy as np
#import IPython
import ast, sys, os,csv, getopt
import csv
import ntpath
from multiprocessing import Pool

current=os.path.dirname(os.path.realpath(__file__))

def main(argv):

    if len(sys.argv)==2:
        cores = sys.argv[1]
        if not RepresentsInt(cores):
            cores=4
        else:
            cores=int(cores)
    else:
        cores=4
    print "Cores used: ", cores
    rna_lib_folder=sys.argv[3]
    protein_lib_folder=sys.argv[2]
    print "RNA library folder: ", os.path.join(rna_lib_folder)
    print "Protein library folder: ", os.path.join(protein_lib_folder)

#    IPython.embed()
    pool = Pool(processes=cores)
    for prot_file in os.listdir(os.path.join(current,protein_lib_folder)):
        if not prot_file.startswith("."):
            libchecker(os.path.join(current,protein_lib_folder,prot_file))
    for rna_file in os.listdir(os.path.join(current,rna_lib_folder)):
        if not rna_file.startswith("."):
            libchecker(os.path.join(current,rna_lib_folder,rna_file))

    for prot_file in os.listdir(os.path.join(current,protein_lib_folder)):
        if not prot_file.startswith("."):
            for rna_file in os.listdir(os.path.join(current,rna_lib_folder)):
                if not rna_file.startswith("."):
                    #IPython.embed()
                    pool.apply_async(f, args=(os.path.join(current,rna_lib_folder,rna_file),os.path.join(current,protein_lib_folder,prot_file)))

    pool.close()
    pool.join()
#    IPython.embed()

def f(rna_file,prot_file):
    prot_library=[]
    prot_frag_names=[]
    with open(prot_file, "r") as tsv:
        for line in csv.reader(tsv, dialect="excel-tab"):
            prot_frag_names.append(line[0].strip())
            line=[float(k) for k in line[1:51]]
            prot_library.append(line)
    prot_frag_names=prot_frag_names[0::10]
    prot_library=np.array(prot_library)


    rna_frag_names=[]
    rna_library=[]
    with open(rna_file, "r") as tsv:
        for line in csv.reader(tsv, dialect="excel-tab"):
            rna_frag_names.append(line[0].strip())
            line=[float(k) for k in line[1:51]]
            rna_library.append(line)
    rna_frag_names=rna_frag_names[0::10]
    rna_library=np.array(rna_library)


    with open(os.path.join(current,'template/params/coefficients.pro.txt'), "r") as f:
        prot_coef=[float(k.strip()) for k in f.readlines()]
    prot_coef=np.array(prot_coef)
    prot_coef=prot_coef/float(max(np.abs(prot_coef)))
    f.close()

    with open(os.path.join(current,'template/params/coefficients.rna.txt'), "r") as f:
        rna_coef=[float(k.strip()) for k in f.readlines()]
    rna_coef=np.array(rna_coef)
    rna_coef=rna_coef/float(max(np.abs(rna_coef)))
    f.close()

    coef=[]
    with open(os.path.join(current,'template/params/coefficients.txt'), "r") as f:
        lines=f.readlines()
    for line in lines:
            coef.append([float(k) for k in line.strip().split()])
    coef=np.array(coef)
    f.close()


    prot_fragments_summed=[]
    for i in xrange(0,len(prot_library),10):
        prot_fragments_summed.append(np.dot(prot_library[i:i+10].T,prot_coef))
    prot_fragments_summed=np.array(prot_fragments_summed)


    rna_fragments_summed=[]
    for i in xrange(0,len(rna_library),10):
        rna_fragments_summed.append(np.dot(rna_library[i:i+10].T,rna_coef))
    rna_fragments_summed=np.array(rna_fragments_summed)

#    IPython.embed()
    with open(os.path.join(current,'pre-compiled/out.merged.'+ntpath.basename(prot_file).split('.')[0]+'-'+ntpath.basename(rna_file).split('.')[0]+'.txt'), "w") as writef:
        for pf_index in range(0,len(prot_fragments_summed)):
            table=[]
            for i in range(0,len(prot_fragments_summed[pf_index])):
                table.append(coef[i]*prot_fragments_summed[pf_index][i])
            table=np.array(table)
            for rf_index in range(0,len(rna_fragments_summed)):
                writef.writelines(prot_frag_names[pf_index]+" "+rna_frag_names[rf_index]+" "+ str(round(sum(np.dot(table,rna_fragments_summed[rf_index])),2))+" -1\n")


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

def RepresentsInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

if __name__ == "__main__":

    main(sys.argv[1:])
