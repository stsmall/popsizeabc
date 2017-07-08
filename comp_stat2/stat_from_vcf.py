# This script computes the AFS, LD and IBS summary statistics described in the manuscript, from a sample of genomes at vcf format.
# The list of individuals used for the computation (a subset of the individuals included in the vcf file) must be provided in a separate file.
# The population of each individual in the vcf must be provided in the first column of a file at ped format.

# Summary statistics are stored in a file with suffix .stat, which has a single column with one statistic per line.
# Statistics are the same as those produced by simul_data.py.

# For this example the vcf file is Chr1.vcf.gz, the list of animals is list_indiv_Jersey.txt, the ped file is indiv.ped and the result file is example_vcf_Jersey_n30_mac6_macld6.stat.
# All these names can easily be changed by modifying the script below.

################################################
############### dependencies ###################
################################################

#!/usr/bin/python

import sys
from gwas import IO
from gwas import data
from gwas import summary_stat as ss
import numpy as np
import time

################################################
############### paramaters #####################
################################################

pop='Jersey'
list_ani=IO.read_list('../cattle_data/list_indiv_Jersey.txt') # list of diploid animals used for computing the summary statistics
n=len(list_ani)*2 # haploid sample size
mac=6 # minor allele count threshold for AFS and IBS statistics computation
mac_ld=6 # minor allele count threshold for LD statistics computation
L=2000000 # size of each segment, in bp.

# time windows
nb_times=21 # number of time window
Tmax=130000 # the oldest time window will start at Tmax
a=0.06 # the length of time windows increases when time increases, at a speed that is determined by this coefficient  
# computation of time windows based on the above parameters
times=-np.ones(shape=nb_times,dtype='float')
for i in range(nb_times):
    times[i]=(np.exp(np.log(1+a*Tmax)*i/(nb_times-1))-1)/a
print "Population size changes at the following times (in generations):"
print times
print ""

# LD statistics parameters 
per_err=5 # the length of each interval, as a percentage of the target distance
r=10**(-8) # the per generation per bp recombination rate (an approximation is sufficient value)
# creation of the bins of physical distance for which the average LD will be computed, based on the time windows defined above.
interval_list=[]
for i in range(nb_times-1):
    t=(times[i+1]+times[i])/2
    d=1/(2*r*t)
    if d <= L:
        interval_list.append([d-per_err*d/100,d+per_err*d/100])
t=Tmax+times[nb_times-1]-times[nb_times-2]
d=10**8/(2*t)
#d=1/(2*r*t)
interval_list.append([d-per_err*d/100,d+per_err*d/100])
print "Average LD will be computed for the following distance bins (in bp) :"
print interval_list
print ""

################################################
############### the programm ###################
################################################

print('Number of diploid individuals : '+str(n/2)+'\n')
nb_dist=len(interval_list)
print('Total number of statistics : '+str(nb_dist+n/2+1)+'\n')

print 'Started the analysis, current time :'+time.ctime()
sys.stdout.flush()

# summarize each chromosome
count_list=[]
pos_list=[] # for ld
geno_list=[] # for ld
nb_snp=0
Lchr=0
for chro in range(1,2):
    # store data and pre-analyse it
    print 'Processing chromosome ',chro
    infile_vcf='../cattle_data/Chr'+str(chro)+'.vcf.gz'
    [mydata,mymap]=IO.parseVcfFile(infile_vcf,includeInd=list_ani)
    IO.parsePedFile_nogeno('../cattle_data/indiv.ped',mydata)  
    u=data.Genos_and_counts(mydata,mymap,pop,mac=mac)
    count_list.append(u[1][0])
    p=len(u[0][0])
    nb_snp+=p
    Lchr+=u[0][0][p-1]
    u=data.Genos_and_counts(mydata,mymap,pop,mac=mac_ld)
    pos_list.append(u[0][0])
    geno_list.append(u[2][0])
    sys.stdout.flush()

# compute summary statistics
res_afs=ss.histo(count_list,n/2)
print 'AFS computed, current time :'+time.ctime()
sys.stdout.flush()
u=ss.break_chr(pos_list,geno_list,L)
pos_list=u[0]
geno_list=u[1]
print 'Broke chromosomes in independent segments of size '+str(L)+', current time :'+time.ctime()
sys.stdout.flush()
res_ld_zyg=ss.distrib_zyg_r2(pos_list,geno_list,interval_list)
print 'LD statistics computed, current time :'+time.ctime()
sys.stdout.flush()

# print the result
np.savetxt('res/example_vcf_'+pop+'_n'+str(n)+'_mac'+str(mac)+'_macld'+str(mac_ld)+'.stat',np.concatenate((np.array([np.float(nb_snp)/np.float(Lchr)],dtype='float'),res_afs[0],res_ld_zyg[0])),fmt='%.3e')
print "Finished job, ","current time :", time.ctime()
