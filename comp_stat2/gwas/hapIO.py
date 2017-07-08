# This file includes routines to read haplotype data with diverse formats (here beagle and ms) and store them in the class hapdata (see hapdata.py)
# It can also write data from class hapdata into diverse formats (only MSMC here).

import sys
import hapdata
import re
import numpy as np
import numpy.ma as ma
from gwas import missing
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

#################### beagle files ######################

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

def parseBeagleFiles_allSNP(DatasetName,MapName,PhaseName,includeInd=[],indels=False):
    '''
    Creates an haplotype dataset based on :
      - a map file with p lines and 7 columns (p is the number of polymorphic sites), similar to example_chr1.map.gz, at gz format
      - a phased Beagle file with (p+1) lines and (2+2n) columns (n is the number of individuals), similar to example_chr1.phase.gz, at gz format
    Options:
      - includeInd: store info only for these individuals.
      - indels : if True, indels will be included.
    '''
    # count the snps
    nsnp=0
    for line in os.popen('zcat '+MapName):
	buf=line.split()
	if (indels) or (not indels and buf[6] == 'S'): 
	    nsnp+=1
    print nsnp,' SNP found'
    ficphase=os.popen('zcat '+PhaseName)
    # count / create a hash of / create a sorted list of indices for - valid individuals
    line=ficphase.readline()
    buf=line.split()
    if buf[0] == 'I':
	nindiv=0
        allIndiv={}
        individx=[] # provides the column number for the first haplotype of each indiv
        for i in range(2,len(buf)):
	    indiv=buf[i]	
	    if indiv in includeInd:
		try:
        	    myindiv=allIndiv[indiv]
    		except KeyError: 
	            nindiv+=1
                    allIndiv[indiv]=i # not used for the moment
                    allIndiv[key+str(i)]=indiv
	            individx.append(i)
        print nindiv, 'individuals stored'
    else:
	print 'Wrong format : missing header line'
	raise ValueError
    # create individuals
    dataset=hapdata.HapDataset(DatasetName,nsnp=nsnp,nhaplo=2*nindiv)
    for i in individx:
	dataset.addIndividual(ID=buf[i])
    # create SNP
    Map=hapdata.Map()
    for line in os.popen('zcat '+MapName):
	buf=line.split()
	if (indels) or (not indels and buf[6] == 'S'):
	    snp_name=buf[0]
	    dataset.addSnp(snp_name)
	    dataset.snp[snp_name].initAlleles(buf[2],buf[3])
	    buf2=snp_name.split(":")
	    Map.addMarker(M=snp_name,C=buf2[0],posP=float(buf[1]))
    # create haplotypes
    line=ficphase.readline()
    ficmap=os.popen('zcat '+MapName)
    line2=ficmap.readline()
    while line != "":
	buf=line.split()	
	buf2=line2.split()
	snp_name=buf[1]
	if (indels) or (not indels and buf2[6] == 'S'):
	    for i in individx:
	        all1=dataset.snp[snp_name].alleles[int(buf[i])-1]
	        all2=dataset.snp[snp_name].alleles[int(buf[i+1])-1]
		dataset.SetIndivHaplotypes(allIndiv[key+str(i)],snp_name,all1,all2)
	line=ficphase.readline()
	line2=ficmap.readline()
    return [dataset,Map]

#################### ms files ######################

def create_pos(pos1,L):
    '''
    creates an array of positions betwen 1 and L, from an array of positions between 0 and 1
    '''
    pos2=np.round(L*pos1)
    pos2[0]=max(1,pos2[0])
    for i in range(1,len(pos2)):
	if pos2[i]<=pos2[i-1]:
	    pos2[i]=pos2[i-1]+1
    return pos2    

def parsems_onepop_nofile_fast(DatasetName,output,n2=np.inf):
    '''
    Creates an haplotype dataset from an ms output string with one single pop.
    Command line includes at least flags -t and -r
    Assumes a single chromosome
    n2 is the target sample size. If not specified, the actual sample size is used.
    '''
    lines=output.splitlines()
    # first line - general parameters
    buf=lines[0].split()
    n=int(buf[1])
    if n2>n:
	n2=n
    rep=int(buf[2])
    L=int(buf[7])
    # skip next line
    k=1
    # skip one blank line and the line with \\ and read the one with segsites
    k+=3
    buf=lines[k].split()
    nb_snp=int(buf[1])
    if nb_snp>0:
	# creates dataset
    	dataset=hapdata.HapDataset(DatasetName,nsnp=nb_snp,nhaplo=n2)
	# reads snp positions
	k+=1
	buf=lines[k].split()
	dataset.snp_pos=create_pos(np.array(buf[1:],dtype='float'),L)
	# read haplotypes
	k+=1
	k2=0
	#sel_hap=np.random.choice(n,size=n2,replace=False)
	sel_hap=np.random.permutation(n)
	sel_hap=sel_hap[0:n2]
	for i in sel_hap:
	    dataset.UpdatePop(k2,'pop1')
	    dataset.Data[k2,:]=np.array(list(lines[k+i])[:nb_snp],dtype='int16')
	    k2+=1
    else:
	# creates empty dataset
    	dataset=hapdata.HapDataset(DatasetName,nsnp=0,nhaplo=n2)
    return dataset

def write_msmc_fast(mydata,output):
    '''
    writes the MSMC input file for a dataset extracted from ms - fast version
    '''
    u=hapdata.Haplos_and_counts_fast(mydata,'pop1')
    pos=u[0]
    hap=u[2]
    fic=open(output,'w')
    pos0=0
    code=['A','C']
    for i in range(len(pos)):
	fic.write('1\t'+str(int(pos[i]))+'\t'+str(int(pos[i]-pos0))+'\t')
	for j in hap[:,i]:
	    fic.write(code[j])
	fic.write('\n')
	pos0=pos[i]
    fic.close()


