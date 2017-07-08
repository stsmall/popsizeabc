# This function loads the .param and .stat files produced by simul_data.py or simul_stat.py, and performs some basic operations.
# The difference with load_simu.R is that two stat files are loaded, one for AFS and IBS statistics, the other for LD statistics
# This allows to combine AFS and LD statistics with a different MAF threshold
 
# Several parameters need to be defined before calling this function :
# - infile_params
# - infile_stat
# - infile_stat_ld
# - n
# - mac

nrows <- -1
params_all <- read.table(infile_params,nrows=nrows)
stat_all <- read.table(infile_stat,nrows=nrows)
stat_all_ld=read.table(infile_stat_ld,nrows=nrows)

# Parameter names
colnames(params_all) <- c("m","r",paste("N",round(generations_numbers),sep=""))

# Statistic names
d=rep("?",nb_dist)
ind=which(ld_distances<10000)
d[ind]=paste(round(ld_distances[ind]),'b',sep="")
d[-ind]=paste(round(ld_distances[-ind]/1000),'kb',sep="")

var_names1=c("AFS0",paste("AFS_",1:(n/2),sep=""),paste("LD_zyg_",d,sep=""))
probs=c(0.0001,0.001,0.01,0.1,0.25,0.5,0.75,0.9,0.99,0.999,0.9999)
nb_prob=length(probs)
m=c(1)
nb_m=length(m)
var_names2=paste("IBS",m[1],"_",probs,sep="")
colnames(stat_all) <- c(var_names1,var_names2)

# replaces LD statistics in stat_all bu those in infile_stat_ld
stat_all[,n/2+1+(1:nb_dist)]=stat_all_ld[,n/2+1+(1:nb_dist)]
rm(stat_all_ld)

# remove stats for allele counts below mac
if (mac>1){
	ind_rem=2:mac
	stat_all=stat_all[,-ind_rem]
}

# remove missing values
gr <- good.rows(stat_all)
stat <- stat_all[gr,]
params <- params_all[gr,]

print(paste("stat matrix dimensions :",dim(stat)[1],dim(stat)[2],sep=" "))
print(paste("params matrix dimensions :",dim(params)[1],dim(params)[2],sep=" "))

