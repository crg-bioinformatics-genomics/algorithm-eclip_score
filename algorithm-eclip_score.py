#!/usr/bin/env python
import argparse
import yaml
import os
import subprocess
import shutil, errno
import sys

# we want to be agnostic to where the script is ran
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
WORKER_PATH = os.path.realpath(os.curdir)

# Function to copy the output directory content
def copyfolder(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc: # python >2.5
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise

# read the task definition yaml file
with open(os.path.join(SCRIPT_PATH, "catrapid_fragments_ultra.yaml"), "r") as task_f:
	task_definition = yaml.load(task_f)

parser = argparse.ArgumentParser(
   description='Launches catRAPID fragments with properly parset parameters')

# parser.add_argument(
#    '-fileA', type=str, default=["none"], nargs=1, help='Dataset A')
#
# parser.add_argument(
#    '-fileB', type=str, default=["none"], nargs=1, help='Dataset B')

parser.add_argument(
   '-output_dir', type=str, nargs=1,
   help='Directory where the output is going to be stored')

# accept form fields
for item in task_definition['form_fields']:
   nargs = 1 if item['required'] else "?"
   parser.add_argument(
       '-FORM%s'%item['name'], type=str, default=["none"], nargs=nargs,
       help='Form argument: %s' % item)

# this parse_args stops if any unexpected arguments is passed
args = parser.parse_args()

OUTPUT_PATH = os.path.join(WORKER_PATH, args.output_dir[0])

# import code; code.interact(local=locals())
# import IPython
# IPython.embed()

import re
import StringIO
from Bio import SeqIO

Ppat = re.compile('>.*?\n[ARNDCQEGHILKMFPSTWYV]+', re.IGNORECASE)
if Ppat.match(args.FORMprotein_seq[0]) == None:
	args.FORMprotein_seq[0] = ">input_protein\n"+args.FORMprotein_seq[0]
protSeq = []
for record in SeqIO.parse(StringIO.StringIO(args.FORMprotein_seq[0]), "fasta"):
	protSeq.append(record)

protFile = os.path.join(OUTPUT_PATH.replace("output/", ""),"protein.fasta")
output_handle = open(protFile, "w")
SeqIO.write(protSeq, output_handle, "fasta")
output_handle.close()

Rpat = re.compile('>.*?\n[GATCU]+', re.IGNORECASE)
if Rpat.match(args.FORMrna_seq[0]) == None:
	args.FORMrna_seq[0] = ">input_rna\n"+args.FORMrna_seq[0]
rnaSeq = []
for record in SeqIO.parse(StringIO.StringIO(args.FORMrna_seq[0]), "fasta"):
	rnaSeq.append(record)

rnaFile = os.path.join(OUTPUT_PATH.replace("output/", ""),"rna.fasta")
output_handle = open(rnaFile, "w")
SeqIO.write(rnaSeq, output_handle, "fasta")
output_handle.close()

os.chdir(SCRIPT_PATH)
# print(WORKER_PATH)

random_number = (""" "{}" """.format(args.output_dir[0])).split("/")[3]

args.FORMtitle = "".join([t.replace(' ', '_') for t in args.FORMtitle])
# os.rename(os.path.join(WORKER_PATH,args.fileA[0]), os.path.join(WORKER_PATH,args.fileA[0].replace(' ', '-')))
# os.rename(os.path.join(WORKER_PATH,args.fileB[0]), os.path.join(WORKER_PATH,args.fileB[0].replace(' ', '-')))

if type(args.FORMtitle)==list:
	title = "none"
else:
	title = args.FORMtitle.replace(" ", "_")

logfile = open("pylog."+str(random_number)+".txt","w")

# cmd = """bash runcatrapid.fragments.sh "{}" "{}" "{}" "{}" "{}" "{}" "{}" """.format(random_number, args.FORMemail, title, args.FORMprotFrag[0], "150", os.path.join(WORKER_PATH,args.fileA[0].replace(' ', '-')), os.path.join(WORKER_PATH,args.fileB[0].replace(' ', '-')) )
cmd = """bash algorithm-eclip_score.sh "{}" "{}" "{}" "{}" "{}" "{}" "{}" "{}" """.format(random_number, args.FORMemail[0], title, args.FORMprotFrag[0], "150", protFile, rnaFile, args.FORMextras[0])

p = subprocess.Popen(cmd, cwd=SCRIPT_PATH, shell=True)

logfile.write(str(p.returncode)+"\n")

p.wait()

logfile.write(str(p.returncode)+"\n")

if p.returncode == 0:

	TMP_PATH = "./tmp/{}/outputs/".format(random_number)
	logfile.write(TMP_PATH+"\n")
	dirList=os.listdir(TMP_PATH)
	shutil.copyfile(os.path.join(OUTPUT_PATH.replace("output/", ""),"protein.fasta"), OUTPUT_PATH+"protein.fasta")
	shutil.copyfile(os.path.join(OUTPUT_PATH.replace("output/", ""),"rna.fasta"), OUTPUT_PATH+"rna.fasta")
	logfile.write("copied fastas\n")
	for file in dirList:
		if os.path.isfile(TMP_PATH+file): #ignore directories
			shutil.copyfile(TMP_PATH+file, OUTPUT_PATH+file)
			logfile.write("copied all files\n")
	if os.path.exists(OUTPUT_PATH+"images") == False :
		copyfolder(SCRIPT_PATH+"/images", OUTPUT_PATH+"images")
		logfile.write("copied images folder\n")

	from django.template import Template
	from django.template import Context
	from django.conf import settings
	from django.template import Template

	settings.configure(TEMPLATE_DIRS=(os.path.join(SCRIPT_PATH,'./')), DEBUG=True, TEMPLATE_DEBUG=True)

	logfile.write(args.FORMprotFrag[0]+"\n")
	if args.FORMprotFrag[0] == "uniform":
	#	# read the template file into a variable
		with open(os.path.join(OUTPUT_PATH, "index.mat.tmp.html"), "r") as template_file:
			template_string = "".join(template_file.readlines())
		logfile.write("copied index.mat.tmp.html\n")

	if args.FORMprotFrag[0] == "weighted":
		with open(os.path.join(OUTPUT_PATH, "index.Z.tmp.html"), "r") as template_file:
			template_string = "".join(template_file.readlines())
		logfile.write("copied index.Z.tmp.html\n")

	import datetime

	# create template from the string
	t = Template(template_string)

	# context contains variables to be replaced
	c = Context(
		{
			"title": title,
			"rnaFrag" : "150",
			"randoms" : random_number,
# 			"fileA" : protFile,
# 			"fileB" : rnaFile,
			"generated" : str(datetime.datetime.now()),
		}
	)

	# and this bit outputs it all into index.html
	with open(os.path.join(OUTPUT_PATH, "index.html"), "w") as output:
		output.write(t.render(c))
	logfile.write("created index.html\n")

else:
	sys.exit("The execution of the bash script failed.")
	logfile.write("bash script failed\n")

logfile.write("that's it!\n")
logfile.close()
os.remove("pylog."+str(random_number)+".txt")
