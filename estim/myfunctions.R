# This file includes auxiliary functions

good.rows <- function(data1){
index.positive <- data1>=0
return(!as.logical(rowSums(!index.positive)))
}

filter.negative <- function(data1, data2){
good.rows1 <- good.rows(data1)
return(data.frame(stat=I(data1[good.rows1,]), params=I(data2[good.rows1,])))
}

square_distance_var <- function(data1, data2){
return(colMeans((data1-data2)**2)/apply(data1,2,var))
}

log10_pop <- function(parameters){
return(cbind(parameters[,1:2],log10(parameters[,-1:-2])))
}
