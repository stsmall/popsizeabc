# This script computes the AFS, LD and IBS summary statistics described in the manuscript, from a set of genomic samples that have already been produced by the script simul_data.py (with option save_ms=True). It can be used for instance to compute summary statistics using different MAF thresholds, or a different sample size, than in the first analysis.
# To call this script from a terminal, go into the 'abc_code' directory and run the command :
# python simul_stat.py

# This script produces no .param files, but a .stat file that is similar to that produced by simul_data.py

# With the current parameter values, the only difference with simul_data.py is that we use a MAF threshol of 20%, similar to what is done in the manuscript for the cattle data analysis.
# This can easily be changed by modifying variables in the parameters section below.
# In practice, it is strongly recommended to excute this script using many processors in parallel (with different output names!) and to merge all output files at the end.
 
################################################
############### dependencies ###################
################################################

#!/usr/bin/python

import sys
import time
import datetime
import numpy as np
import popgen_abc

start_time=time.time()

################################################
############### paramaters #####################
################################################

# general parameters
outfile_name='res/example_simu' # root of the output files. This must be the same as that used in simul_data.py.
nb_rep=100 # number of simulated datasets. This must be the same as that used in simul_data.py.
nb_seg=100 # number of independent segments in each dataset. This must be the same as that used in simul_data.py.
L=2000000 # length of one segment, in bp. This must be the same as that used in simul_data.py.
n=50 # haploid sample size. This must be the same as that used in simul_data.py.

n2=50 # number of haploid genomes used to compute summary statistics. must be lower than n and even.
nb_seg2=100 # number of independent segments used to compute summary statistics. must be lower than nb_seg.
mac=10 # minor allele count threshold for AFS and IBS statistics computation
mac_ld=10 # minor allele count threshold for LD statistics computation

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
# creation of the bins of physical distance for which the average LD will be computed, based on the time windows defined above
interval_list=[]
for i in range(nb_times-1):
    t=(times[i+1]+times[i])/2
    d=1/(2*r*t)
    if d <= L:
        interval_list.append([d-per_err*d/100,d+per_err*d/100])
t=Tmax+times[nb_times-1]-times[nb_times-2]
d=10**8/(2*t)
interval_list.append([d-per_err*d/100,d+per_err*d/100])
print "Average LD will be computed for the following distance bins (in bp) :"
print interval_list
print ""

################################################
############### the programm ###################
################################################

# create the matrices where results (parameters and statistics) are stored
nb_dist=len(interval_list)
results=-np.ones(shape=[nb_rep,nb_dist+n2/2+1],dtype='float') # missing statistic values will be set to 1. 
									      # this mostly concerns LD statistics, because it is sometimes impossible to find SNP pairs in a given distance bin.
print('Total numberb of statistics : '+str(nb_dist+n2/2+1)+'\n')

# compute summary statistics
outfile_name_2_read = outfile_name + "_n" + str(n) + "_s" + str(nb_seg)
outfile_name_2_write = outfile_name + "_n" + str(n2) + "_s" + str(nb_seg2) + "_mac" + str(mac) + "_macld" + str(mac_ld)
print 'Started the computation'
for i in range(nb_rep):
    elapsed_time=time.time()-start_time
    print 'Computing statistics for replicate', i+1,", current time :", time.ctime(),", elapsed time :", datetime.timedelta(seconds=elapsed_time)
    sys.stdout.flush()
    try:
        results[i,:]=popgen_abc.comp_stats_one_rep_macld(outfile_name_2_read,i,nb_seg2,L,n2,interval_list,mac=mac,mac_ld=mac_ld)
    except:
        print 'Problem with replicate', i+1
        pass
# print the result
np.savetxt(outfile_name_2_write+'.stat',results[0:nb_rep,:],fmt='%.3e')
print "Printed the results"

# the end
elapsed_time=time.time()-start_time
print "Finished job, ","current time :", time.ctime(), "elapsed time :", datetime.timedelta(seconds=elapsed_time)

