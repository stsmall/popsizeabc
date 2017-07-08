# This file includes routines to read genotype data in vcf format and store them in the class data (see data.py)
# It can also write data from class hapdata into diverse formats (only MSMC here).

import sys
import data
import re
import numpy as np
import numpy.ma as ma
from bisect import insort
from gwas import missing
from gwas import summary_stat as ss
from scipy.stats import bernoulli
import gzip
import os


defaultParams= {
    "MissPop":'0',
    "MissGeno":'0',
    "MissPheno":'-9',
    "MissParent":'0',
    "MissSex":'0'
    }

key="kshuwewekiuvf"

#################### vcf files ######################

def parsePedFile_nogeno(fileName,dataset,params=defaultParams):
    '''
    Read a ped file without genotype columns to update individual information in an existing dataset
    '''
    fic=open(fileName)
    print 'Reading',fileName
    for ligne in fic:
	buf=ligne.split()
	indiv_name=buf[1]
	if indiv_name in dataset.indiv:
	    pop=((buf[0]!=params['MissPop']) and buf[0]) or None
	    fatherID=((buf[2]!=params['MissParent']) and buf[2]) or None
            motherID=((buf[3]!=params['MissParent']) and buf[3]) or None
            sex=(buf[4]!=params['MissSex'] and buf[4]) or None
	    pheno=((buf[5]!=params['MissPheno']) and buf[5]) or None
	    dataset.updateIndividual(indiv_name,pop=pop,sex=sex,fatherID=fatherID,motherID=motherID,phenotype=pheno)           

def read_list(fileName):
    '''
    Read a one column text file and stores it in a list (for instance of individuals)
    '''
    fic=open(fileName)
    print 'Reading',fileName
    res=[]
    for ligne in fic:
	buf=ligne.split()
	if len(buf) > 1:
	    print 'warning : only the first column will be read'
	res.append(buf[0])
    return res

def parseVcfFile(fileName,includeInd=[],gz=True):
    '''
    Read a genotype file at vcf format
    options are:
      - includeInd: store info for these individuals
    '''
    if gz:
	fic=gzip.open(fileName,'rb')
    else:
	fic=open(fileName,'r')
    print 'Reading',fileName
    # count snps
    nsnp=0
    for ligne in fic:
	nsnp+=1
    fic.seek(0)
    ligne=fic.readline()
    while ligne[0:2]=='##':
	nsnp-=1
	ligne=fic.readline()
    if ligne[0] == '#':
	nsnp-=1
    else:
	print 'wrong format : missing header line'
	raise ValueError
    print nsnp,' SNP found'
    # count / create a hash of / create a sorted list of indices for - valid individuals
    buf=ligne.split()
    nindiv=0
    allIndiv={}
    individx=[]
    for i in range(9,len(buf)):
	indiv=buf[i]	
	if indiv in includeInd:
	    nindiv+=1
            allIndiv[indiv]=i # not used for the moment
            allIndiv[key+str(i)]=indiv
	    individx.append(i)
    print nindiv, 'individuals stored'
    # create individuals
    dataset=data.Dataset(fileName,nsnp=nsnp,nindiv=nindiv)
    for i in individx:
	dataset.addIndividual(ID=buf[i])
    #print ligne
    # create snps and genotypes
    ligne=fic.readline()
    Map=data.Map()
    while ligne != "":
	buf=ligne.split()
	if buf[2]=='.':
	    snp_name=buf[0]+":"+buf[1]
	else:
	    snp_name=buf[2]	
	dataset.addSnp(snp_name)
	alt=buf[4].split(',')
	dataset.snp[snp_name].initAlleles(buf[3],alt[0]) # if multiple alternative alleles, only the first one is considered
	Map.addMarker(M=snp_name,C=buf[0],posP=float(buf[1]))
	for i in individx:
	    try:
		s=buf[i]
	    except IndexError:
		print ligne
	    try:
		all1=dataset.snp[snp_name].alleles[int(s[0])]
	    	all2=dataset.snp[snp_name].alleles[int(s[2])]
		dataset.SetGenotype(allIndiv[key+str(i)],snp_name,all1+all2)
	    except:
	        dataset.SetGenotype(allIndiv[key+str(i)],snp_name,'toto') # will be considered as missing value
	ligne=fic.readline()
    return [dataset,Map]

