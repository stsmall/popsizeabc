# Defines time windows and distance bins for LD statistics, similar to that in simul_data.py

nb_times=21
Tmax=13000
generations_numbers=rep(-1,nb_times)
a=0.06
for (i in 0:(nb_times-1)){
    generations_numbers[i+1]=(exp(log(1+a*Tmax)*i/(nb_times-1))-1)/a
}
gen_time=1

generations_centers <- rowMeans(cbind(generations_numbers,c(generations_numbers[-1],Tmax+2*Tmax-2*generations_numbers[nb_times-1])))
generations_log_centers <- rowMeans(log10(cbind(c(1,generations_numbers[-1]),c(generations_numbers[-1],Tmax+2*Tmax-2*generations_numbers[nb_times-1]))))
generations <- c(0,log10(generations_numbers[-1]))

years=c(0,log10(gen_time*generations_numbers[-1]))
years_centers=c((years[1:(nb_times-1)]+years[2:nb_times])/2,6)

ld_distances=10^8/(2*generations_centers)
ind=which(ld_distances<=2000000)
ld_distances=ld_distances[ind]
nb_dist=length(ld_distances)
