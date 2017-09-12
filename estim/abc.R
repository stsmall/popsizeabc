rm(list=ls())
library(abc)
source('/home/ssmall2/Desktop/popsizeSTS/generations.R')
source('/home/ssmall2/Desktop/popsizeSTS/myfunctions.R')

# loads simulated samples, MAF 20% for both AFS and LD stats
infile_params="/home/ssmall2/Desktop/popsizeSTS/wbpng.n40_s40.params"
infile_stat="/home/ssmall2/Desktop/popsizeSTS/wbpng.n40_s40_mac3_macld3.stat"
n=40 # haploid sample size
mac=3 # minor allele count
source("/home/ssmall2/Desktop/popsizeSTS/load_simu.R")

# loads observed samples
infile_obs="/home/ssmall2/Desktop/popsizeSTS/popsizeSTS_obs.stats"
source("/home/ssmall2/Desktop/popsizeSTS/load_obs.R")

# choice of summary stats
nb_afs=n/2-mac+2
ind_afs=1:nb_afs
ind_ld=nb_afs+(1:(nb_dist-1)) # the LD statistic corresponding to the shortest distance is removed
ind_stat=c(ind_afs, ind_ld) # only AFS and LD statistics

# abc
abc_res=abc(obs[ind_stat],log10_pop(params)[,-1],stat[,ind_stat],tol=0.1,method="neuralnet") # a tolerance of order 0.001 would give much better results, but this would require more simulated samples.
abc_estim=summary(abc_res,print=FALSE)[3,] # median
#abc_estim=summary(abc_res,print=FALSE)[5,] # mode

# plot
pdf('/home/ssmall2/Desktop/ex1_Jersey_mac6_macld6_estim.pdf',height=4,width=4)
par(mar=c(3,3,1,2),cex=0.7,mgp=c(1.5,0.5,0))
plot(NA,xlim=c(0,5),ylim=c(2,5),xlab="years before present (log scale)",ylab="effective population size (log scale)",axes=F)
axis(1,at=0:5,labels=c("0","10","100","1,000","10,000","100,000"))
axis(2,at=c(2,log10(300),3,log10(3000),4,log10(30000),5),labels=c("100","300","1,000","3,000","10,000","30,000","100,000"))
lines(years,abc_estim[-1],type="s",lwd=2) # the first term of the vector is removed because it corresponds to the recombination rate
abc_q=summary(abc_res, print=FALSE, intvl=0.9)[2,] # 5% quantile
print(abc_q)
lines(years,abc_q[-1],type="s",lwd=2,lty=3)
abc_q=summary(abc_res,print=FALSE,intvl=0.9)[6,] # 95% quantile
print(abc_q)
lines(years,abc_q[-1],type="s",lwd=2,lty=3)
dev.off()

# plot detailed information about the estimation of all parameters
plot(abc_res,param=log10_pop(params)[,-1],file='/home/ssmall2/Desktop/ex1_Jersey_mac6_macld6_detailed',subsample=20)
