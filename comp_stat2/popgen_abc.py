#!/usr/bin/python
import sys
import os
from gwas import hapIO
from gwas import summary_stat as ss
from gwas import hapdata
import numpy as np
import subprocess
import tarfile
import msprime

def simul_stats_one_rep_macld(outfile_name,rep,nb_seg,L,n,times,params,interval_list,mac=1,mac_ld=10,save_msp=False,prob_err=0):
    '''
    Simulates one sample with msprime and returns a vector of summary statistics
    stats : nb SNP, folded AFS, LD
    Each allele is switched ( 0 becomes 1 and 1 becomes 0) with probability prob_err (default is 0).
    '''
    # creates demographic events
    mymut=params[0]
    myrec=params[1]
    N=params[2]
    mydemo=[]
    for i in range(1,len(times)):
	mydemo.append(msprime.PopulationParametersChange(time=times[i],initial_size=params[i+2]))
    # simulate and summarize each segment
    count_list=[]
    pos_list=[] # for ld
    hap_list=[] # for ld
    nb_snp=0
    for i in range(nb_seg):
	tree_sequence = msprime.simulate(sample_size=n, Ne=N, length=L, recombination_rate=myrec, mutation_rate=mymut, demographic_events=mydemo)
    	p=tree_sequence.get_num_mutations()
    	if p>0:
	    # creates dataset
    	    mydata=hapdata.HapDataset('ms_rep',nsnp=p,nhaplo=n)
	    # reads snp positions
	    pos=[]
	    for mut in tree_sequence.mutations():
		pos.append(mut[0])
	    mydata.snp_pos=np.array(pos,dtype='float')
	    # read haplotypes
	    j=0
	    for h in tree_sequence.haplotypes():
	    	mydata.UpdatePop(j,'pop1')
	    	mydata.Data[j,:]=np.array(list(h)[:p],dtype='int16')	    
	    	j+=1
	    # stores snp positions and haplotypes if required
	    if save_msp==True:
		outfile=open(outfile_name+'_rep'+str(rep+1)+'_seg'+str(i+1)+'.msp','w')
		for mut in tree_sequence.mutations():
		    outfile.write(str(mut[0])+' ')
		outfile.write('\n')
		for h in tree_sequence.haplotypes():
	    	    outfile.write(h+'\n')	        
    	    # updates lists
	    if prob_err>0:
		mydata.introduce_errors(prob_err)
            u=hapdata.counts_fast(mydata,'pop1',mac=mac)       
            count_list.append(u)
	    nb_snp+=len(u)
	    u=hapdata.Haplos_and_counts_fast(mydata,'pop1',mac=mac_ld)
	    pos_list.append(u[0])
            hap_list.append(u[2])
	else:
	    print 'no segsite in segment', i+1, 'of replicate', rep+1
    geno_list=ss.hap_to_geno(hap_list)
    # computation of summary statistics
    res_afs=ss.histo(count_list,n/2)
    res_ld_zyg=ss.distrib_zyg_r2(pos2_list,geno2_list,interval_list)
    return np.concatenate((np.array([np.float(nb_snp)/(L*nb_seg)],dtype='float'),res_afs,res_ld_zyg[0]))

def comp_stats_one_rep_macld(outfile_name,rep,nb_seg,L,n2,interval_list,mac=1,mac_ld=10,prob_err=0):
    '''
    Reads haplotype files and returns a vector of summary statistics
    stats : nb SNP, folded AFS, LD
    The file for segment j is called outfile_name_repi_segj.msp (where i is given by argument rep)
    Summary statistics are computed for n2 randomly sampled sequences and the first nb_seg segments
    Each allele is switched (0 becomes 1 and 1 becomes 0) with probability prob_err (default is 0).
    '''
    # open tar archive
    mytar=tarfile.open(outfile_name+'.msp.tar.bz2','r:bz2')
    # simulate and summarize each segment
    count_list=[]
    pos_list=[] # for ld
    hap_list=[] # for ld
    nb_snp=0
    for i in range(nb_seg):
	try:
	    # computes initial sample size
	    output=mytar.extractfile(outfile_name+'_rep'+str(rep+1)+'_seg'+str(i+1)+'.msp').read()
	    lines=output.splitlines()
	    n=len(lines)-1
	    sel_hap=np.random.choice(n,size=n2,replace=False)
	    # number if snps
	    buf=lines[0].split()
	    p=len(buf)	    
	    # creates dataset
    	    mydata=hapdata.HapDataset('ms_rep',nsnp=p,nhaplo=n2)
	    # reads snp positions
	    mydata.snp_pos=np.array(buf,dtype='float')
	    # read haplotypes
	    for j in range(n2):
	    	mydata.UpdatePop(j,'pop1')
	    	mydata.Data[j,:]=np.array(list(lines[sel_hap[j]+1])[:p],dtype='int16')
	    if prob_err>0:
		mydata.introduce_errors(prob_err)
	    # updates lists
            if prob_err>0:
		mydata.introduce_errors(prob_err)
            u=hapdata.counts_fast(mydata,'pop1',mac=mac)       
            count_list.append(u)
	    nb_snp+=len(u)
	    u=hapdata.Haplos_and_counts_fast(mydata,'pop1',mac=mac_ld)
	    pos_list.append(u[0])
            hap_list.append(u[2])
	except IOError:
	    print 'File '+outfile_name+'_rep'+str(rep+1)+'_seg'+str(i+1)+'.msp does not exist'
    mytar.close()
    geno_list=ss.hap_to_geno(hap_list)
    # computation of summary statistics
    res_afs=ss.histo(count_list,n2/2)
    res_ld_zyg=ss.distrib_zyg_r2(pos2_list,geno2_list,interval_list)
    return np.concatenate((np.array([np.float(nb_snp)/(L*nb_seg)],dtype='float'),res_afs,res_ld_zyg[0]))

