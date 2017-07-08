rm(list=ls())
library(abc)
source('generations.R')
source('myfunctions.R')

# loads simulated samples, MAF 20% for both AFS and LD stats
infile_params="../simu_stat/simu_n50_s100.params"
infile_stat="../simu_stat2/simu_n50_s100_mac10_macld10.stat"
n=50 # haploid sample size
mac=10 # minor allele count
source("load_simu.R")

# loads observed samples
infile_obs="../cattle_stat/Holstein_n50_mac10_macld10.stat"
source("load_obs.R")

# choice of summary stats
nb_afs=n/2-mac+2
ind_afs=1:nb_afs
ind_ld=nb_afs+(1:(nb_dist-1)) # the LD statistic corresponding to the shortest distance is removed
ind_ibs=nb_afs+nb_dist+(1:(nb_m*nb_prob))
ind_stat=c(ind_afs,ind_ld) # only AFS and LD statistics

# abc
abc_res=abc(obs[ind_stat],log10_pop(params)[,-1],stat[,ind_stat],tol=0.005,method="neuralnet")
abc_estim=summary(abc_res,print=FALSE)[3,] # median
#abc_estim=summary(abc_res,print=FALSE)[5,] # mode

# plot
pdf('estim_Holstein_MAF20.pdf',height=4,width=4)
par(mar=c(3,3,1,2),cex=0.7,mgp=c(1.5,0.5,0))
plot(NA,xlim=c(0,5),ylim=c(2,5),xlab="years before present (log scale)",ylab="effective population size (log scale)",axes=F)
axis(1,at=0:5,labels=c("0","10","100","1,000","10,000","100,000"))
axis(2,at=c(2,log10(300),3,log10(3000),4,log10(30000),5),labels=c("100","300","1,000","3,000","10,000","30,000","100,000"))
lines(years,abc_estim[-1],type="s",lwd=2) # the first term of the vector is removed because it corresponds to the recombination rate 
lines(c(4,4),c(2.1,5),lwd=3,lty=3)
text(x=4,y=2,labels="domestication")
dev.off()

# plot detailed information about the estimation of all parameters
plot(abc_res,param=log10_pop(params)[,-1],file='estim_Holstein_MAF20_detailed')

# sample parameter values from the posterior distribution
nb_rep=10
h=hist(abc_res,breaks=50,plot=FALSE)
post_params=array(data=0,dim=c(nb_rep,dim(params)[2]))
post_params[,1]=10^(-8) # fixed mutation rate
for (rep in 1:nb_rep){
	post_params[rep,2]=sample(h[[1]]$mids,1,prob=h[[1]]$density) # recombination rate
	post_params[rep,3]=sample(h[[2]]$mids,1,prob=h[[2]]$density) # population size in the most recent time window
	for (j in 4:dim(params)[2]){ # population size in other time windows
		post_params[rep,j]=sample(h[[j-1]]$mids,1,prob=h[[j-1]]$density)
		d=abs(post_params[rep,j]-post_params[rep,j-1])
		while (d>1){
			post_params[rep,j]=sample(h[[j-1]]$mids,1,prob=h[[j-1]]$density)
			d=abs(post_params[rep,j]-post_params[rep,j-1])
		}
	}
}
post_params[,-(1:2)]=10**post_params[,-(1:2)]
write.table(post_params,file='estim_Holstein_MAF20.post_params',quote=F,row.names=F,col.names=F)




