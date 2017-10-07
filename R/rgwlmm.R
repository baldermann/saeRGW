# Outlier-robust estimation for the Geographically weighted linear mixed model (GWLMM)
#
# User interface for fitting a GWLMM model when the estimation process is distorted
# by extreme values (outliers).
#' @inheritParams gwlmm
#' @param k (numeric) defines the tuning constant for influience function (default = 1.346). Ses datails.
#' @param method (character) defines the iterative algorithm for approximating the variance
#' parameters. Possible values: \code{"Fix", "Newton"}, (default = \code{"Fix"}). See datails.
#' @param Start (list) optioanl list containing three obejcts for the starting
#' values for the robust approximation (default = NULL).
#' The three objects are: \code{betas} (Matrix with local coefficients);
#' \code{sigma.v} (numeric value for the variance of the random effects);
#' \code{sigma.e} (numeric value for the error term variance).
#'
#'
#'
#'@details
#'\itemize{
#'\item  \code{rgwlmm} uses robust ML estimation equations (following Sinha and Rao, 2009)
#' of the GWLMM that restrict the infuence of outliers on the parameter estimation using
#' Hubers's influence function (See Baldermann et al. (2016)). The tuning constant \code{k} defines the restriction. The smaller \code{k}, the stronger
#' the restriction and vice versa.Sinha and Rao (2009) recommended to set \code{k = 1.346}
#' following Huber (1964).
#' \item  \code{method}: using a fixed-point algorithm is recommended for estimating the variance parameters as it
#' is fast and stable compared to the Newton-Raphsom algorithm.
#'
#'}
#'
#'
#'@examples
#'##################################################################
#'# Outlier-robust estimation
#'\dontrun{
#'
#'# Model fit
#'rgwmodel<- rgwlmm(formula, data = sampleData)
#'
#'# In-sample prediction
#'rpred<-predict(rgwmodel)
#'#Small area preditions (mean) for aggregated population data
#'rpredagg<-predict(rgwmodel, popdata = popaggData, size = "Size")
#'#Small area preditions (mean) for unit-level population data
#'rpreddisagg<-predict(rgwmodel, popdata = popoutData, popAgg = FALSE)
#'
#'###########
#'
#'# Robust model fit when sample only contains centroid information
#'rgwmodel<- rgwlmm(formula, data = sampleData, centroid = TRUE)
#'
#'# In-sample prediction
#'rpred<-predict(rgwmodel)
#'#Small area means for aggregated population data
#'rpredagg<-predict(rgwmodel, popdata = popaggData, size = "Size")
#'
#'}
#'
#@rdname rfit
#@export
#rgwlmm <- function(formula, ...) UseMethod("rgwlmm")
#' @rdname fit
#' @export
rgwlmm <- function(formula,data, band=NULL, centroid=FALSE, maxit=100, tol=0.0001, k=1.345, Start =NULL, method = "Fix", ...){
    ff<- Formula(formula)
    clusterid1<- model.part(ff,data=data,rhs=2)[,1]
    data<-data[order(clusterid1),]
    mf.s<-model.frame(formula=ff, data = data)
    Y<- model.response(mf.s)
    X<- model.matrix(ff,data=mf.s,rhs=1)
    clusterid<- model.part(ff,data=mf.s,rhs=2)[,1]
    geoInfo <- model.part(ff,data=mf.s,rhs=3)
    #Centroid data set
      if(centroid==TRUE){
          temp <- cbind(geoInfo)
          centroids<-aggregate(temp, list(clusterid=clusterid), mean, na.rm=TRUE)
          rm(list=c("temp"))
        }
    Z  <-dummies::dummy(clusterid)
    colnames(Z)<-sort(unique(clusterid))
    ni_1              <-table(clusterid)
    n0                <-which(ni_1!=0)
    ni                <-ni_1[n0]
    #Number of areas in the sample. Does not have to be the same as in the Population (in case of non-sampled areas)
    m<-       length(unique(clusterid))
    n<-       sum(ni)

      # Starting values for variance parameters
       if(is.null(Start)){
            ex.lme     <-lmer(Y~0+X + (1|clusterid))
            sigmasq0.v           <- as.numeric(attr(VarCorr(ex.lme)$clusterid, "stddev"))^2
            if(sigmasq0.v<=0){sigmasq0.v<-1}
             sigmasq0             <- as.numeric(attr(VarCorr(ex.lme), "sc"))^2
            if(sigmasq0<=0){sigmasq0<-1}
             }else{
                 sigmasq0.v <-Start$sigma.v
                 sigmasq0   <-Start$sigma.e
             }

        # Compute Distances between
        if(centroid==TRUE){
            centroid.s      <-centroids[which(centroids$clusterid %in%clusterid  ==TRUE),c(-1)]
            distance        <-as.matrix(dist(centroid.s))
            uniquePos       <-lapply(sort(unique(clusterid)), where, x=clusterid)
              if(is.null(band)){
                      band    <-gwr.sel(Y ~ 0+X,coords=do.call(cbind,(lapply(centroid.s, rep, times = ni))))
              }else{band<-band}
        }else{
            distance        <-as.matrix(dist(geoInfo))
            uniquePos       <-as.list(1:n)
              if(is.null(band)){
                    band    <-gwr.sel(Y ~ 0+X ,coords=as.matrix(geoInfo))
                }else{band<-band}
        }

      #Weighting matrix
      W       <-as.data.frame(lapply(as.data.frame(gwr.Gauss(distance^2,band)),
                              function(x){ifelse(x<=1e-06,min(x[x>1e-06]),x)}))
      if(centroid ==TRUE){W<-as.data.frame(lapply(W,rep,times = ni))}

      #Starting values for local Betas:
        if(is.null(Start)){
          Var                 <- vi.inv(n=ni, e=sigmasq0, u=sigmasq0.v, Z=Z)
          estGwr<-do.call(rbind,mapply(gwrEst, W,Pos=uniquePos,
                      MoreArgs=list(X=X,y=Y, Inv=Var$V.Inv, V=Var$V), SIMPLIFY=FALSE))
          beta.gwr<-(estGwr[,1:dim(X)[2]])
          beta.gwr0         <-beta.gwr
        }else{
         beta.gwr <- Start$betas
         beta.gwr0<-beta.gwr
        }

        #Estimation of robust parameters
        Parameters<-ConvGWRobust(method = method, beta.gwr0=beta.gwr0,
                                       sigmasq0.v=sigmasq0.v,sigmasq0=sigmasq0,
                                       n =n,ni=ni,m=m,tol=tol,maxit=maxit,k=k,X=X,W=W,Y=Y,
                                       Z=Z, Cpp=FALSE, Pos = uniquePos)

        colnames(Parameters$betas)<-attr(X,"dimnames")[[2]]
        colnames(Parameters$vbeta2)<-attr(X,"dimnames")[[2]]

        result<- list(
                   call = match.call(),
                   Coefficients = Parameters$betas,
                   VarCoefficients = Parameters$vbeta2,
                   Variance = list("residual"= Parameters$sigmasq,
                                   "raneff"  =Parameters$sigmasq.v),
                   bandwidth = band,
                   nIter = Parameters$n.iter,
                   Model  = list(mf.s = mf.s, formula = formula,
                            Projection = Parameters$proj,
                            vbeta = Parameters$vbeta,
                            Weights = W,
                            centroid = centroid,
                            tunConst = k)
                    )
        class(result)<-"rgwlmm"
        result
}
#######################################################################################################################################
