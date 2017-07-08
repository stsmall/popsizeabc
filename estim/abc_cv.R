rm(list=ls())
library(abc)
source('generations.R')
source('myfunctions.R')

nval=100 # number of tested histories

# loads simulated samples, all SNP for AFS stats and MAF 20% for LD stats
infile_params="../simu_stat/simu_n50_s100.params"
infile_stat="../simu_stat/simu_n50_s100_mac1_macld1.stat"
infile_stat_ld="../simu_stat2/simu_n50_s100_mac10_macld10.stat"
n=50 # haploid sample size
mac=1 # minor allele count for AFS and LD statistics
source("load_simu_combine.R")

# choice of summary stats
nb_afs=n/2-mac+2
ind_afs=1:nb_afs
ind_ld=nb_afs+(1:(nb_dist-1)) # the LD statistic corresponding to the shortest distance is removed
ind_ibs=nb_afs+nb_dist+(1:(nb_m*nb_prob))
ind_stat=c(ind_afs,ind_ld) # only AFS and LD statistics

# abc
res=cv4abc(param=log10_pop(params)[,-1],sumstat=stat[,ind_stat],nval=nval,tols=c(0.005),statistic="median",method="ridge") # ridge is used here in order to reduce computation time
par.estim=(res$estim)[[1]]
par.true=res$true
errors=colSums((par.estim-par.true)**2)/(length(par.true[,1])*var(par.true[,dim(par.true)[2]]))

# plot
pdf("cv_errors.pdf",height=4,width=4)
par(mar=c(3,3,1,2),cex=0.7,mgp=c(1.5,0.5,0))
plot(NA,xlim=c(0,5.5),ylim=c(0,1),xlab="generations before present (log scale)",ylab="population size prediction error ",xaxt="n")
axis(1,at=0:5,labels=c("0","10","100","1,000","10,000","100,000"))
y=errors[-1]
points(c(generations,5.5),c(y,y[nb_times]),type="s",lwd=2)
dev.off()


