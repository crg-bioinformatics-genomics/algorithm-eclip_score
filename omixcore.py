#!/usr/bin/env python
import argparse
import yaml
import os
import subprocess
import shutil, errno
import sys
import json
from django.template import Template
from django.template import Context
from django.conf import settings
from django.template import Template
import datetime
import IPython
# we want to be agnostic to where the script is ran
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
WORKER_PATH = os.path.realpath(os.curdir)

# Function to copy the output directory content
def copyfolder(src, dst):
	try:
		if os.path.exists(dst):
			shutil.rmtree(dst)
		shutil.copytree(src, dst)
	except OSError as exc: # python >2.5
		if exc.errno == errno.ENOTDIR:
			shutil.copy(src, dst)
		else: raise

def error_page(OUTPUT_PATH,SCRIPT_PATH,title,random_number,error_message):
	settings.configure(TEMPLATE_DIRS=(os.path.join(SCRIPT_PATH,'./')), DEBUG=True, TEMPLATE_DEBUG=True)
	with open(os.path.join(SCRIPT_PATH, "index.error.html"), "r") as template_file:
		   template_string = "".join(template_file.readlines())

	# create template from the string
	t = Template(template_string)

	c = Context(
	{"title": title,
	 "randoms" : random_number,
	 "generated" : str(datetime.datetime.now()),
	 "error_message" : error_message,
	 })

	with open(os.path.join(OUTPUT_PATH, "index.html"), "w") as output:
		output.write(t.render(c))


# read the task definition yaml file
with open(os.path.join(SCRIPT_PATH, "omixcore.yaml"), "r") as task_f:
	task_definition = yaml.load(task_f)

parser = argparse.ArgumentParser(
   description='Launches')

# parser.add_argument(
#	'-fileA', type=str, default=["none"], nargs=1, help='Dataset A')
#
# parser.add_argument(
#	'-fileB', type=str, default=["none"], nargs=1, help='Dataset B')

parser.add_argument(
   '-output_dir', type=str, nargs=1,
   help='Directory where the output is going to be stored')


# accept form fields
for item in task_definition['form_fields']:
	nargs = 1 if item['required'] else "?"
	'''if item['name']=="rna_seq":
		nargs='*'
	'''
   	parser.add_argument(
	   '-FORM%s'%item['name'], type=str, default="none", nargs=nargs,
	   help='Form argument: %s' % item)

# this parse_args stops if any unexpected arguments is passed
args = parser.parse_args()

OUTPUT_PATH = os.path.join(WORKER_PATH, args.output_dir[0])


import re
import StringIO
from Bio import SeqIO

Ppat = re.compile('>.*?\n[ARNDCQEGHILKMFPSTWYV]+', re.IGNORECASE)
if Ppat.match(args.FORMprotein_seq[0]) == None:
	args.FORMprotein_seq[0] = ">input_protein\n"+args.FORMprotein_seq[0]
protSeq = []
for record in SeqIO.parse(StringIO.StringIO(args.FORMprotein_seq[0]), "fasta"):
	protSeq.append(record)
	break


args.FORMrna_seq=[args.FORMrna_seq]
Rpat = re.compile('>.*?\n[GATCU]+', re.IGNORECASE)
rnaSeq = []
records=SeqIO.parse(StringIO.StringIO(args.FORMrna_seq[0]), "fasta")

for record in records:
	if record.id=="":
		record.id="input_rna"
	check=">"+record.id+"\n"+str(record.seq)
	if Rpat.match(check) != None:
		rnaSeq.append(record)


rnafolder=os.path.join(OUTPUT_PATH.replace("output/", ""),"rnas")
if not os.path.exists(rnafolder):
    os.makedirs(rnafolder)


valid_entries=0
rnaAllFile=os.path.join(OUTPUT_PATH,"rna.fasta")
output_all_handle = open(rnaAllFile, "w")
for entry in rnaSeq:
	if len(entry.seq)>=500:
		valid_entries+=1
		entry.id=re.sub('[^0-9a-zA-Z]+', '_', entry.id)
		entry.description=""
		rnaFile = os.path.join(rnafolder,entry.id)
		output_handle = open(rnaFile, "w")
		SeqIO.write(entry, output_handle, "fasta")
		SeqIO.write(entry, output_all_handle, "fasta")
		output_handle.close()

protFile = os.path.join(OUTPUT_PATH.replace("output/", ""),"protein.fasta")
output_handle = open(protFile, "w")
for protein_record in SeqIO.parse(StringIO.StringIO(args.FORMprotein_seq[0]), "fasta"):
	output_handle.write(str(">input_protein\n"+protein_record.seq))
	break

# SeqIO.write(protSeq, output_handle, "fasta")
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

modeFile=open(os.path.join(OUTPUT_PATH,"mode"),'w')
modeFile.writelines(args.FORMmode[0])
modeFile.close()

if args.FORMmode[0]=="custom" and valid_entries==0:
	error_message="Custom option selected, with 0 valid entries. Please check the lengths of the custom transcript sequences to be at least 500nt and re-submit!"
	error_page(OUTPUT_PATH,SCRIPT_PATH,title,random_number,error_message)
	logfile.write("custom option and valid entries. Created error index.html\n")
	sys.exit()
	#sys.exit("Wrong Submission. Custom option selected, with 0 valid entries. The execution of the bash script failed.")

if len(protein_record.seq)<=150:

	error_message="Protein sequence too short! Please check the length of the protein sequence to be at least 150aa and re-submit!"
	error_page(OUTPUT_PATH,SCRIPT_PATH,title,random_number,error_message)
	logfile.write("Protein sequence too short. Created error index.html\n")
	sys.exit()
	#sys.exit("Wrong Submission. Protein sequence too short! The execution of the bash script failed.")


cmd = """bash omixcore.sh "{}" "{}" "{}" "{}" "{}" "{}" "{}"  """.format(random_number, args.FORMemail[0], title, "150", protFile, args.FORMmode[0],rnafolder)

p = subprocess.Popen(cmd, cwd=SCRIPT_PATH, shell=True)

logfile.write(str(p.returncode)+"\n")

p.wait()

logfile.write(str(p.returncode)+"\n")

if p.returncode == 0:

	TMP_PATH = "./tmp/{}/outputs/".format(random_number)
	logfile.write(TMP_PATH+"\n")
	dirList=os.listdir(TMP_PATH)
	shutil.copyfile(os.path.join(OUTPUT_PATH.replace("output/", ""),"protein.fasta"), OUTPUT_PATH+"protein.fasta")

	logfile.write("copied fastas\n")
	for file in dirList:
		if os.path.isfile(TMP_PATH+file): #ignore directories
			shutil.copyfile(TMP_PATH+file, OUTPUT_PATH+file)
			logfile.write("copied all files\n")

	if os.path.exists(TMP_PATH+"rna_libs") == True :
		copyfolder(TMP_PATH+"rna_libs", OUTPUT_PATH+"rna_libs")
		logfile.write("copied rna_libs folder\n")

	if os.path.exists(OUTPUT_PATH+"images") == False :
		copyfolder(SCRIPT_PATH+"/images", OUTPUT_PATH+"images")
		logfile.write("copied images folder\n")



	settings.configure(TEMPLATE_DIRS=(os.path.join(SCRIPT_PATH,'./')), DEBUG=True, TEMPLATE_DEBUG=True)




	with open(os.path.join(OUTPUT_PATH,"Signature_prediction.txt"), "r") as sign_result:
		PredictionScore=float(sign_result.readlines()[0])
	if PredictionScore<0.5:
		with open(os.path.join(SCRIPT_PATH, "index.nrbp.html"), "r") as template_file:
			   template_string = "".join(template_file.readlines())

		# create template from the string
		t = Template(template_string)
		c = Context(
		{"title": title,

		 "PredictionScore" : PredictionScore,
		 "randoms" : random_number,
		 "generated" : str(datetime.datetime.now()), })

		with open(os.path.join(OUTPUT_PATH, "index.html"), "w") as output:
 			output.write(t.render(c))
 		logfile.write("created nrbp index.html\n")

	else:
		with open(os.path.join(SCRIPT_PATH, "index.mat.html"), "r") as template_file:
			template_string = "".join(template_file.readlines())
		logfile.write("copied index.mat.tmp.html\n")

		# create template from the string
		t = Template(template_string)

		# context contains variables to be replaced
		c = Context(
			{
				"title": title,
				"rnaFrag" : "150",
				"randoms" : random_number,
				 "fileA" : protFile,
	#			 "fileB" : rnaFile,
				"generated" : str(datetime.datetime.now()),
			}
		)

		# and this bit outputs it all into index.html
		rnas=[]
		values=[]
		with open(os.path.join(OUTPUT_PATH,"filter.processed.txt"), "r") as fltr:
			for line in fltr:
				if len(line.split())==2:
					rnas.append(line.split()[0].strip())
					values.append(line.split()[1].strip())
		fltr.close()
		with open(os.path.join(OUTPUT_PATH,"transcript.rows"),'w') as fobj:
			json.dump(zip(rnas,values),fobj)
		fobj.close()
		with open(os.path.join(OUTPUT_PATH, "index.html"), "w") as output:
			output.write(t.render(c))
		logfile.write("created index.html\n")

else:
	error_message="Webserver failed to run prediction. Please try again!"
	error_page(OUTPUT_PATH,SCRIPT_PATH,title,random_number, error_message)

	logfile.write("bash script failed\n")
	sys.exit()


logfile.write("that's it!\n")
logfile.close()
os.remove("pylog."+str(random_number)+".txt")
