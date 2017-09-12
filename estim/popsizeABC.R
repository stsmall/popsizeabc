#!/usr/bin/env Rscript
# load libs
library("optparse")
library("abc")
rm(list=ls())
source('myfunctions.R')

# load options
option_list = list(make_option(c("-s", "--stats"), type="character", default=NULL, help="filename of simulated stats", metavar="character"),
                   make_option(c("-o", "--obs"), type="character", default=NULL, help="filename of observed stats", metavar="character"),
                   make_option(c("-p", "--params"), type="character", default=NULL, help="filename of parameters stats", metavar="character"),
                   make_option(c("-n", "--haps"), type="character", default=NULL, help="number of haploid inds", metavar="character"),
                   make_option(c("-m", "--mac"), type="integer", default=3, help="minor allele count cutoff, [default= %default]", metavar="character"),
                   make_option(c("-g", "--gens"), type="integer", default=1, help="generation time, [default= %default]", metavar="character"),
                   make_option(c("-t", "--tmax"), type="integer", default=13000, help="time maximum,[default= %default]", metavar="character"),
                   make_option(c("-u","--mu"), type="double", default=2.9E-9, help="mutation rate per bp"),
                   make_option("--outfile", type="character", default="out.pdf"),
                   make_option("--scaled", action="store_true"),
                   make_option("--crossval", action="store_true"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Defines time windows and distance bins for LD statistics, similar to that in simul_data.py
nb_times=21
Tmax=opt$tmax
generations_numbers=rep(-1, nb_times)
a=0.06
for (i in 0:(nb_times-1)){
  generations_numbers[i+1]=(exp(log(1+a*Tmax)*i/(nb_times-1))-1)/a
}
gen_time=opt$gens
generations_centers <- rowMeans(cbind(generations_numbers,c(generations_numbers[-1],Tmax+2*Tmax-2*generations_numbers[nb_times-1])))
generations_log_centers <- rowMeans(log10(cbind(c(1,generations_numbers[-1]),c(generations_numbers[-1],Tmax+2*Tmax-2*generations_numbers[nb_times-1]))))
generations <- c(0,log10(generations_numbers[-1]))
years=c(0,log10(gen_time*generations_numbers[-1]))
years_centers=c((years[1:(nb_times-1)]+years[2:nb_times])/2,6)
ld_distances=10^8/(2*generations_centers)  # 10^8 is default to correspond with 1E-8
ind=which(ld_distances<=2000000)  # fragments size
ld_distances=ld_distances[ind]
nb_dist=length(ld_distances)

# loads simulated sample
infile_params=opt$params
infile_stat=opt$stats
n=opt$haps
mac=opt$mac
nrows <- -1
params_all <- read.table(infile_params, nrows=nrows)
stat_all <- read.table(infile_stat, nrows=nrows)
colnames(params_all) <- c("m", "r", paste("N", round(generations_numbers), sep=""))
# Statistic names
d=rep("?", nb_dist)
ind=which(ld_distances<10000)
d[ind]=paste(round(ld_distances[ind]), 'b', sep="")
d[-ind]=paste(round(ld_distances[-ind]/1000), 'kb', sep="")
var_names1=c("AFS0", paste("AFS_",1:(n/2),sep=""), paste("LD_zyg_", d, sep=""))
## IBS
#probs=c(0.0001,0.001,0.01,0.1,0.25,0.5,0.75,0.9,0.99,0.999,0.9999)
#nb_prob=length(probs)
#m=c(1)
#nb_m=length(m)
#var_names2=paste("IBS",m[1],"_",probs,sep="")
#colnames(stat_all) <- c(var_names1,var_names2)
colnames(stat_all) <- c(var_names1)
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

# loads observed samples
infile_obs=opt$obs
obs=read.table(infile_obs)
obs=obs[,1]
# remove stats for allele counts below mac
if (mac>1){
  ind_rem=2:mac
  obs=obs[-ind_rem]
}
print(paste("stat vector dimension :",length(obs),sep=" "))

# choice of summary stats
nb_afs=n/2-mac+2
ind_afs=1:nb_afs
ind_ld=nb_afs+(1:(nb_dist-1)) # the LD statistic corresponding to the shortest distance is removed
#ind_ibs=nb_afs+nb_dist+(1:(nb_m*nb_prob))
ind_stat=c(ind_afs, ind_ld) # only AFS and LD statistics
# abc
abc_res=abc(obs[ind_stat], log10_pop(params)[,-1], stat[,ind_stat], tol=0.005, method="neuralnet") 
# a tolerance of order 0.001 would give much better results, but this would require more simulated samples.
abc_estim=summary(abc_res, print=FALSE)[3,]  # median
#abc_estim=summary(abc_res,print=FALSE)[5,]  # mode
# plot
pdf(paste(opt$outfile, ".pdf"), height=4, width=4)
par(mar=c(3,3,1,2),cex=0.7,mgp=c(1.5,0.5,0))
plot(NA, xlim=c(0,5), ylim=c(2,5), xlab="years before present (log scale)", ylab="effective population size (log scale)", axes=F)
axis(1, at=0:5, labels=c("0", "10", "100", "1,000", "10,000", "100,000"))
axis(2, at=c(2, log10(300), 3, log10(3000), 4, log10(30000), 5), labels=c("100", "300", "1,000", "3,000", "10,000", "30,000", "100,000"))
lines(years, abc_estim[-1], type="s", lwd=2) # the first term of the vector is removed because it corresponds to the recombination rate
abc_q=summary(abc_res, print=FALSE, intvl=0.9)[2,] # 5% quantile
print(abc_q)
lines(years, abc_q[-1], type="s", lwd=2, lty=3)
abc_q=summary(abc_res, print=FALSE, intvl=0.9)[6,] # 95% quantile
print(abc_q)
lines(years, abc_q[-1], type="s", lwd=2, lty=3)
dev.off()
# plot detailed information about the estimation of all parameters
plot(abc_res, param=log10_pop(params)[,-1], file=paste(opt$outfile, "-detailed.pdf"), subsample=20)

# cross validation
if (opt$crossval == TRUE){
  nval=100 # number of tested histories
  # abc
  res=cv4abc(param=log10_pop(params)[,-1], sumstat=stat[,ind_stat], nval=nval, tols=c(0.005), statistic="median", method="neuralnet") 
  par.estim=(res$estim)[[1]]
  par.true=res$true
  errors=colSums((par.estim-par.true)**2)/(length(par.true[,1])*var(par.true[,dim(par.true)[2]]))
  # plot
  pdf("cv_errors.pdf", height=4, width=4)
  par(mar=c(3,3,1,2), cex=0.7, mgp=c(1.5,0.5,0))
  plot(NA, xlim=c(0,5.5), ylim=c(0,1), xlab="generations before present (log scale)", ylab="population size prediction error ", xaxt="n")
  axis(1, at=0:5, labels=c("0", "10", "100", "1,000", "10,000", "100,000"))
  y=errors[-1]
  points(c(generations, 5.5), c(y, y[nb_times]), type="s", lwd=2)
  dev.off()
}


