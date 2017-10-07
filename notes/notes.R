library(MASS)
library(lme4)
library(spgwr)
library(dummies)
# library(Hmisc)
 # library(spatstat)
# require(Rcpp)
# require(RcppArmadillo)
 require(Formula)

source("R/Functions.R")


#here i try different things
# library(plyr)
load(file="notes\\notesSpNstNoOut.RData")
l=1
sampleData = Pop[[l]][S[[l]],c("long","lat","clusterid", "x", "y")]
popoutData = Pop[[l]][-S[[l]],c("long","lat","clusterid", "x")]
pop<-Pop[[l]][c("long","lat","clusterid", "x", "y")]
Size<-table(pop$clusterid)
popaggData<-aggregate(pop[,-c(3,5)], by=list(pop$clusterid), FUN=mean, na.rm=TRUE)
names(popaggData)[names(popaggData)=="Group.1"] <- "clusterid"
popaggData$Size<-Size
save(sampleData, file=paste("Data\\sampeData.RData",sep=""))
save(popoutData, file=paste("Data\\popoutData.RData",sep=""))
save(popaggData, file=paste("Data\\popaggData.RData",sep=""))


# Tests

# for gwlmmFit

rm(list=ls())



formula <- y~1+x|clusterid |long + lat

# test für centroid = FALSE
test<-gwlmm(formula, data = sampleData, maxit = 3)
pred<-predict(test)
predagg<-predict(test, popdata = popaggData, size = "Size")
preddisagg<-predict(test, popdata = popoutData, popAgg = FALSE)

# test für centroid = TRUE
testc<-gwlmm(formula, data = sampleData, centroid = TRUE, maxit = 3)
predc<-predict(testc)
predaggc<-predict(testc, popdata = popaggData, size = "Size")



#ROBUST ESTIMATION

# test für centroid = TRUE
testrc<-rgwlmm(formula, data = sampleData, centroid = TRUE, method = "Newton")
rpred_c<-predict(testrc)
rpredagg_c<-predict(testrc, popdata = popaggData, size = "Size")

# test für centroid = FALSE
testr<- rgwlmm(formula, data = sampleData, maxit = 3)

rpredagg<-predict(testr, popdata = popaggData, size = "Size")
rpreddisagg<-predict(testr, popdata = popoutData, popAgg = FALSE)









