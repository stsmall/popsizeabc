# This function loads the .stat files produced by stat_from_beagle.py, and performs some basic operations.

# Several parameters need to be defined before calling this function :
# - infile_obs
# - n
# - mac

obs=read.table(infile_obs)
obs=obs[,1]

# remove stats for allele counts below mac
if (mac>1){
	ind_rem=2:mac
	obs=obs[-ind_rem]
}

print(paste("stat vector dimension :",length(obs),sep=" "))

