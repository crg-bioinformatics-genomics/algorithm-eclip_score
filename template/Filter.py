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
from matplotlib import pyplot as plt
# include CM modules
import libprofiles as prof
import libcombinations as lc  


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




    
def features(RBP,RBPlen,RBP_profiles,selectedMotSca,motifprofiles):

	onlyBestCorr=numpy.zeros([1, len(selectedMotSca)])
	p=tp=0	
	for key in (selectedMotSca):
		
		RRM=str(key.split()[1])
		scale=int(key.split()[2])
		bestcorr=float(selectedMotSca[key].split()[1])
		rbpprof=RBP_profiles[RBP][scale]
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

	

AAs = "ARNDCQEGHILKMFPSTWYVX"    
titles, property_table  = prof.parse_scales(  os.path.join("", "./scales/scales_generated.txt"), AAs )    
#print titles

log=open("./log.log","w")
selectedMotSca={}



infile=open("./data/scaletta_meta/C_scaletta_0.75.txt","r") #24 -0.6 ENSP00000378461_58-109 52 0.602412505731066
#####gest best motif+scales
for line in infile:
	if line[0]!="#":
		#print line
		motif=line.split()[2]
		scale=line.split()[3]
		selectedMotSca["C "+motif+" "+scale]=line.split()[4]+" "+line.split()[1] #[motif scale]=cov corr
#print selectedMotSca
log.write(str(len(selectedMotSca))+"\n")
infile=open("./data/scaletta_meta/NC_scaletta_0.75.txt","r") #24 -0.6 ENSP00000378461_58-109 52 0.602412505731066
#####gest best motif+scales
for line in infile:
	if line[0]!="#":
		#print line
		motif=line.split()[2]
		scale=line.split()[3]
		selectedMotSca["NC "+motif+" "+scale]=line.split()[4]+" "+line.split()[1] #[motif scale]=cov corr
#print selectedMotSca
infile=open("./data/scaletta_meta/P_scaletta_0.75.txt","r") #24 -0.6 ENSP00000378461_58-109 52 0.602412505731066
#####gest best motif+scales
for line in infile:
	if line[0]!="#":
		#print line
		motif=line.split()[2]
		scale=line.split()[3]
		selectedMotSca["P "+motif+" "+scale]=line.split()[4]+" "+line.split()[1] #[motif scale]=cov corr
#print selectedMotSca



######calculate motif profile





motifs={}
motifs=create_dic_fromFasta("./data/motifsANDncMotifs30cdhit.txt")
motifprofiles={}

for key in selectedMotSca:
	mm=key.split()[1]
	ss=int(key.split()[2])
	smoothed, norm= cleverMachine_petr(motifs[mm],property_table,AAs)
	motifprofiles[mm+" "+str(ss)]=smoothed[ss]
#print motifprofiles

######calculate pos


posRBP={}
infile=open("./data/inseq.fasta","r")

name=""
for line in infile:
	if line!=">input_protein":
		posRBP["RBP"]=line



posRBP_profiles={}
for key in posRBP:
    #print key,posRBP[key]
    smoothed, norm= cleverMachine_petr(posRBP[key],property_table,AAs)
    posRBP_profiles[key]=smoothed



outfile=open("./outputs/FinalFeatures.txt","w")
#outfile=open("POStoJoinNewFeatures.txt","w")
for RBP in posRBP: 
	print RBP
	outstring="1 "
	RBPlen=len(posRBP[RBP])
	totaldata = numpy.zeros([1, len(selectedMotSca)])


	totaldata= features(RBP,RBPlen,posRBP_profiles,selectedMotSca,motifprofiles)

	for i in range(totaldata.shape[1]):
		outstring+=str(round(totaldata[0,i],2))+" "		
	outfile.write(outstring+"\n")

outfile.close()

