from __future__ import division
import os
import sys
import random
import subprocess
import pickle
import math
import numpy
import scipy
import collections

# include CM modules
import libprofiles as prof

import sklearn
from sklearn.ensemble import RandomForestClassifier
from sklearn import preprocessing
from sklearn import cross_validation
from sklearn import neighbors
from sklearn import metrics
from sklearn import svm
from sklearn import ensemble
from sklearn.externals import joblib


def cleverMachine_petr(test_sequence,property_table,AAs):
	"""
	profile calculations from Petr
	"""
	def process_sequence(sequence, window_size=7):
		"""Single sequence evaluation using scales provided"""
		records = [("seq_submit", sequence)]  # + sequences
		scored = prof.score_proteins(records, AAs, property_table)
		smoothed = map(lambda scored: prof.smooth(scored, window_size=window_size), scored)
		cur_prot = smoothed[0]
		return cur_prot,scored
	# derives scale information
	NO_OF_PROPERTIES = property_table.shape[0]
	# property_groups, ungrouped_properties = lc.property_group_init(NO_OF_PROPERTIES)
	# contains 80x379 entries - 80 scales with 379 entries in each row
	result_smoothed ,prob_resultsnorm= process_sequence(test_sequence)
	return result_smoothed,prob_resultsnorm


def create_dic_fromFasta(inputfilename):
	#log.write(inputfilename)
	#print ("Reading sequences...")
	infile=open(inputfilename,"r")
	name=""
	string=""
	start=0
	dic={}
	num=0
	for line in infile:
		#print line
		if line[0]==">":
			if start==1:
				if string!="Sequence unavailable" and string!="No UUR is annotated for this transcript":
					dic[name]=string
					num+=1
				#print name, num
				string=""
				name=""
				new=""
				name=line[0:len(line)]

				for a in name:
					if a!=">" and a!="\n":
						new+=a
						#let=raw_input(": ")
				name=new
			else:
				start=1
				new=""
				name=line[0:len(line)]
				for a in name:
					if a!=">" and a!="\n":
						new+=a
						#let=raw_input(": ")
				name=new
		else:
			for a in line:
				if a!="\n":
					if a=="T":
						string+="T"
					else:
						string+=a

	if string!="Sequence unavailable" and string!="No UUR is annotated for this transcript":
		dic[name]=string
		num+=1

	return dic


def accum(currentRRMprof, currentRBPprof,protLen):
##calculates the accum of the profile
	correlationField={}
	for a in range(21): #-1...0...1):
		correlationField[a]=0
	outstring=""
	#print correlationField, currentRRMprof, currentRBPprof
	for key in correlationField: ##init dict
		correlationField[key]=0
	for sw in range(0,len(currentRBPprof)-len(currentRRMprof)):
		correlation=numpy.corrcoef(currentRRMprof,currentRBPprof[sw:sw+len(currentRRMprof)] )[0, 1]
		try:
			enterinField=round(correlation,1)
			correlationField[abs(enterinField*10+10)]+=1
		except:
			enterinField=round(numpy.nan_to_num(correlation),1)
			correlationField[(enterinField*10+10)]+=1
		od = collections.OrderedDict(sorted(correlationField.items()))
		dictentry=""
		for k, v in od.iteritems():
			dictentry+=str(v)+"\t"
		accumula=int(dictentry.split()[0])
		outstring=""
		liste=[]
		for i in range(21):
			accumula+=float(dictentry.split()[i])
			liste.append(accumula/int(protLen))
		newliste=[]
		for i in range(len(liste)):
			if (liste[19]!=0):
				newliste.append(liste[i]*100/(liste[19]*100))
			else:
				newliste.append(0)

		for i in newliste:
			outstring+=str(i)+" "

	return outstring

def features(RBP,RBPlen,RBP_profiles,scalet,motifprofiles):

	onlyBestCorr=numpy.zeros([1, len(scalet)])
	p=tp=0
	for key in (scalet):
		RRM=str(key.split()[0])
		scale=int(key.split()[1])
		bestcorr=float(scalet[key].split()[1])
		#print key,RRM,scale,bestcorr,RBP
		rbpprof=RBP_profiles[RBP][scale]
		#print RBP_profiles[RBP][scale]
		#print motifprofiles
		rrmprof=motifprofiles[RRM+" "+str(scale)]

		#calculate accumulation
		accumprofile= accum(rrmprof, rbpprof,RBPlen)
		if accumprofile!="":
			onlyBestCorr[0,tp]=accumprofile.split()[int(bestcorr*10+10)]
			tp+=1
		else:
			onlyBestCorr[0,tp]=0.0
	#print data,onlyBestCorr
	#prediction=BIGmodel.predict_proba(data)
	return onlyBestCorr#,prediction,BIGmodel.classes_



fastafile="./data/inseq.fasta"
print "read fastafile ",fastafile

infile=open(fastafile,"r")
posRBP={}
name=""
for line in infile:
	if line.startswith(">"):
		posRBP["RBP"]=line



AAs = "ARNDCQEGHILKMFPSTWYVX"
titles, property_table  = prof.parse_scales(  os.path.join("", "./scales/scales_generated.txt"), AAs )
#print titles

print "load scan files"
print "load C_scales"
selectedMotSca_short={}
infilename = "./data/scaletta_scan/C_scaletta_0.75.txt"
print infilename
ff=open(infilename,"r") ######16 0.3 ENSP00000304350_1208-1343 17 0.664549863466898
for line in ff:
	if line[0]!="#":
		#print line
		motif=line.split()[2]
		scale=line.split()[3]
		corr=line.split()[1]
		selectedMotSca_short["C "+motif+" "+scale]=corr

print "load NC_scales"
infilename = "./data/scaletta_scan/NC_scaletta_0.75.txt"
print infilename
ff=open(infilename,"r") ######16 0.3 ENSP00000304350_1208-1343 17 0.664549863466898
for line in ff:
	if line[0]!="#":
		#print line
		motif=line.split()[2]
		scale=line.split()[3]
		corr=line.split()[1]
		selectedMotSca_short["NC "+motif+" "+scale]=corr

print "load P_scales"
infilename ="./data/scaletta_scan/P_scaletta_0.75.txt"
print infilename
ff=open(infilename,"r") ######16 0.3 ENSP00000304350_1208-1343 17 0.664549863466898
for line in ff:
	if line[0]!="#":
		#print line
		motif=line.split()[2]
		scale=line.split()[3]
		corr=line.split()[1]
		selectedMotSca_short["P "+motif+" "+scale]=corr


print "load meta files"
print "load C_scales"
selectedMotSca_long={}
infilename = "./data/scaletta_meta/C_scaletta_0.75.txt"
print infilename
ff=open(infilename,"r") ######16 0.3 ENSP00000304350_1208-1343 17 0.664549863466898
for line in ff:
	if line[0]!="#":
		#print line
		motif=line.split()[2]
		scale=line.split()[3]
		corr=line.split()[1]
		selectedMotSca_long["C "+motif+" "+scale]=corr

print "load NC_scales"
infilename = "./data/scaletta_meta/NC_scaletta_0.75.txt"
print infilename
ff=open(infilename,"r") ######16 0.3 ENSP00000304350_1208-1343 17 0.664549863466898
for line in ff:
	if line[0]!="#":
		#print line
		motif=line.split()[2]
		scale=line.split()[3]
		corr=line.split()[1]
		selectedMotSca_long["NC "+motif+" "+scale]=corr

print "load P_scales"
infilename ="./data/scaletta_meta/P_scaletta_0.75.txt"
print infilename
ff=open(infilename,"r") ######16 0.3 ENSP00000304350_1208-1343 17 0.664549863466898
for line in ff:
	if line[0]!="#":
		#print line
		motif=line.split()[2]
		scale=line.split()[3]
		corr=line.split()[1]
		selectedMotSca_long["P "+motif+" "+scale]=corr




print "load pfamModels-pfamID mapping"
infile=open("./data/PfamIDsRnaBinding.txt","r")
cpfam={}
for line in infile:#PFAMID=1
	cpfam[line.split()[0].split(".")[0]]=1


print "load pfamModels-nonclassical mapping"
infile=open("./data/listOfNonClassicalBaltzCastelloKwon.Pfam","r")
ncpfam={}
for line in infile:#PFAMID = 1
	ncpfam[line.split()[0]]=1


print "load pfamModels-putative mapping"
infile=open("./data/listOfunknownPutativeBaltzCastelloKwon.Pfam","r")
putpfam={}
for line in infile:#PFAMID = 1
	putpfam[line.split()[0]]=1



print "load motifs scan"
motifsshort={}
motifsshort=create_dic_fromFasta("./data/motifs.txt")
motifprofiles_short={}
for key in selectedMotSca_short:
	mm=key.split()[1]
	ss=int(key.split()[2])
	smoothed, norm= cleverMachine_petr(motifsshort[mm],property_table,AAs)
	motifprofiles_short[mm+" "+str(ss)]=smoothed[ss]


print "load motifs metapredictor"
motifslong={}
motifslong=create_dic_fromFasta("./data/motifsANDncMotifs30cdhit.txt")
motifprofiles_long={}
for key in selectedMotSca_long:
	mm=key.split()[1]
	ss=int(key.split()[2])
	smoothed, norm= cleverMachine_petr(motifslong[mm],property_table,AAs)
	motifprofiles_long[mm+" "+str(ss)]=smoothed[ss]



super_model = joblib.load("./data/SuperModel/FinalRBFModel.model")
subprocess.call("python Filter.py",shell=True) #creates feature vector for prediction
test=open("./outputs/FinalFeatures.txt","r").readline()

featvec=numpy.zeros([1, 60])
for i in range(1,len(test.split())):
	featvec[0,i-1]=test.split()[i]
#In [168]: pred_model.classes_
#Out[168]: array([-1,  1])
decision=super_model.predict_proba(featvec)[0][1]
print decision


HMMERpresent=0
for RBP in posRBP:

	print "check HMMER"
	o=open("./outputs/HMMER.fasta","w")
	o.write(">"+str(RBP)+"\n")
	o.write(posRBP[RBP]+"\n")
	o.close()
	subprocess.call("hmmscan --cut_ga --noali --domtblout ./outputs/fragmentation.log ../../../Pfam/Pfam-A.hmm ./outputs/HMMER.fasta > ./outputs/hmm_log.txt",shell=True)
	infile=open("./outputs/fragmentation.log","r")
	small={}
	for line in infile:
		if line[0]!="#":
			pfamid=line.split()[1].split(".")[0]
			small[pfamid+" "+line.split()[17]+" "+line.split()[18]]=line.split()[17]+" "+line.split()[18]+" "+line.split()[0] #id, position from to, domname
	infile.close()
	print small
	o=open("./outputs/HMMER.results","w")
	final={}
	if small!=[]:
		for a in small:
			if (a.split()[0] in cpfam) or(a.split()[0] in ncpfam):#or(a.split()[0] in putpfam):
				o.write(a.split()[0]+" "+a.split()[1]+" "+a.split()[2]+" "+small[a].split()[2]+"\n")
				HMMERpresent=1
				final[a]=small[a]
	o.close()
	if os.stat("./outputs/HMMER.results").st_size == 0:
		o=open("./outputs/HMMER.results","w")
		o.write("Nodomain")
		o.close()


	o=open("./outputs/Signature_prediction.txt","w")
	o.write(str(max(round(decision,2),HMMERpresent)))
	o.close()
	print max(round(decision,2),HMMERpresent)


#pred_model.classes_
#Out[74]: array([-1,  1])
