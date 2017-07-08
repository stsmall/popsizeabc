# This function loads the .stat files produced by stat_from_beagle.py, and performs some basic operations.
# The difference with load_obs.R is that two stat files are loaded, one for AFS and IBS statistics, the other for LD statistics
# This allows to combine AFS and LD statistics with a different MAF threshold

# Several parameters need to be defined before calling this function :
# - infile_obs
# - infile_obs_ld
# - n
# - mac

obs=read.table(infile_obs)
obs=obs[,1]
obs_ld=read.table(infile_obs_ld)
obs_ld=obs_ld[,1]

# replaces LD statistics in obs bu those in obs_ld
obs[n/2+1+(1:nb_dist)]=obs_ld[n/2+1+(1:nb_dist)]
rm(obs_ld)

# remove stats for allele counts below mac
if (mac>1){
	ind_rem=2:mac
	obs=obs[-ind_rem]
}

print(paste("stat vector dimension :",length(obs),sep=" "))

