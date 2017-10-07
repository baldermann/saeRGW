#' (Robust) Geographically weighted linear mixed model (GWLMM)
#'
#' User interface for fitting a GWLMM model. The function \code{gwlmm} fits a GWLMM via REML
#' whereas \code{rgwlmm} fits an outlier-robust GWLMM to data. The GWLMM is a random
#' intercespt model that take into account spatial non-stationarity in the model coefficients.
#'
#' @param formula (formula) a two-sided lineer formula object decribing
#' the fixed effects, the random intercept and the geographical information of the model.
#' The response is on the left of the ~ operator, and the x-variables on the right side
#' are seperated by a + opeartor. The identifier for the random intercept is seperated
#' with a vertical bar (|).After another (|), the two coordinates - longitude and latitude -
#' are seperated by a + operator.
#' Categorial x-variable need to be defined as factor-variables.
#' @param data (data.frame) A dataframe containing the variables named in \code{formula}.
#' @param band (numeric) A numeric value defining the bandwidth for the geographical weigths (default = NULL).
#' For a predifined bandwdth, insert value here.
#' @param centroid (logical) If coordinates in \code{formula} are constant within the area ID define cetroid = TRUE (default = FALSE)
#' @param maxit (integer) Defines the maximum number of iterations for the fitting process (default = 100).
#' @param tol (numeric) Defines the tolerance for the convergence of the fitting process (default = 1e-04).
#' @param ... not used
#'
#'
#' @return
#' The function \code{gwlmm} returns an object of class \code{gwlmm}. The function \code{rgwlmm} returns
#' an object of class \code{rgwlmm}. Both objects are lists containing the following elements
#' \itemize{
#' \item \code{call} (language) the call generating the value
#' \item \code{coefficients} (Matrix) the matrix of local regression coefficients
#' \item \code{varcoefficients} (Matrix) the matrix containing the
#' model-based variances of the local coefficients
#' \item \code{Variance} (numeric) the vector of fitted variance parameter(s)
#' \item \code{nIter} (numeric) number of iterations needed
#' for estimating the model parameters
#' \item \code{LRtest} (numeric) a named vector of test results for spatial
#' non-stationarity conating the likelihood, the likelihood ratio,
#'  the effective degrees of freedom and the p-value. (only in \code{gwlmm}-object)
#' \item \code{Model} (list) contains all model information for forther calculations
#' }
#'
#'@details
#'\itemize{
#'\item \code{gwlmm} additionally conducts a likelihood ratio  using the likelihood from the GWLMM
#' and a global LMM. The LMM is fitted using the \code{lmer} function from the \code{lme4}-package.
#' The test is saved in the \code{LRtest} value and is based on Fotheringham et al. (2009, p. 92).
#'
#' \item The bandwidth is optimized using the cross-validation procedure implemented
#'  in \code{gwr.sel} from the \code{spgwr}-package.
#'
#' \item The REML estimation implemented in the \code{gwlmm} function is conducted using the Nelder-Mead algirithm (\code{optim} from \code{stats}-package)
#' as suggested by Chandra et al. (2012).
#'}
#'
#'
#'@references
#'Baldermann, C., N. Salvati, T. Schmid (2016). Robust small area estimation under spatial
#'non-stationarity. Working Paper.
#'
#'Chandra, H., N. Salvati, R. Chambers, and N. Tzavidis (2012). Small area
#'estimation under spatial nonstationarity. Computational Statistics and Data Analysis 56 (10),
#'pp. 2875-2888.
#'
#'Fotheringham, A. S., C. Brunsdon, and M. Charlton (2002). Geographically Weighted Regression. West
#'Sussex: Wiley.
#'
#'Huber, P. J. (1964). Robust estimation of a location parameter. The Annals of Mathematical Statistical 35,
#'73 - 101.
#'
#'Bivand, R and D. Yu (2015). spgwr: Geographically Weighted Regression. R package
#'version 0.6-28. https://CRAN.R-project.org/package=spgw
#'
#'Sinha, S. K. and J. N. K. Rao (2009). Robust small area estimation. The Canadian Journal of Statistics 37
#'(3), 381 - 399.
#'
#'
#'
#'
#'@examples
#'# Data sets ?sampleData, ?popaggData and ?popoutData are
#'# implemented in the rsarGWR-package. See help files.
#'
#'\dontrun{
#'
#'formula <- y~1+x|clusterid |long + lat
#'
#'#Model fit
#'gwmodel<-gwlmm(formula, data = data)
#'
#'#In-sample predictions
#'pred<-predict(gwmodel)
#'
#'#Small area mean prediction for aggregated population data
#'predagg<-predict(gwmodel, popdata = popaggData, size = "Size")
#'
#'#Small area mean prediction for unit-level population data
#'preddisagg<-predict(gwmodel, popdata = popoutData, popAgg = FALSE)
#'}
#'
#@rdname fit
#@export
#gwlmm <- function(formula, ...) UseMethod("gwlmm")

#' @rdname fit
#' @export
#' @seealso \code{\link{predict.gwlmm}}, and \code{\link{predict.rgwlmm}}
gwlmm<-function(formula, data, band=NULL, centroid=FALSE, maxit=100, tol=0.0001, ...){
    ###################################################################################################################
    ###################################################################################################################
    ###################################################################################################################

    ff<- Formula(formula)
    #sort data
    clusterid1<- model.part(ff,data=data,rhs=2)[,1]
    data<-data[order(clusterid1),]

    mf.s<-model.frame(formula=ff, data = data)
    Y<- model.response(mf.s)
    X<- model.matrix(ff,data=mf.s,rhs=1)
    clusterid<- model.part(ff,data=mf.s,rhs=2)[,1]
    geoInfo <- model.part(ff,data=mf.s,rhs=3)

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

    #Number of areas in the sample.
    #Does not have to be the same as in the Population (in case of non-sampled areas)
    m<-       length(unique(clusterid))
    n<-       sum(ni)
    ex.lme               <-lmer(Y~0+X + (1|clusterid), REML=FALSE)



    #Starting values for variance parameters

    sigmasq0.v           <- as.numeric(attr(VarCorr(ex.lme)$clusterid, "stddev"))^2

    if(sigmasq0.v<=0){sigmasq0.v<-1}
    sigmasq0             <- as.numeric(attr(VarCorr(ex.lme), "sc"))^2

    if(sigmasq0<=0){sigmasq0<-1}



    # Compute Distances

    if(centroid==TRUE){
        distance        <-as.matrix(dist(centroids[,c(-1)]))
        uniquePos       <-lapply(sort(unique(clusterid)), where, x=clusterid)
        if(length(band)==0){
        band            <-gwr.sel(Y ~ 0+X,coords=as.matrix(geoInfo))
        }else{
        band=band
        }
    }else{
        distance        <-as.matrix(dist(geoInfo))
        uniquePos       <-as.list(1:n)
        if(length(band)==0){
        band            <-gwr.sel(Y ~ 0+X,coords=as.matrix(geoInfo))
        }else{
        band=band
        }
    }
    W <-as.data.frame(lapply(as.data.frame(gwr.Gauss(distance^2,band)),function(x){ifelse(x<=1e-06,min(x[x>1e-06]),x)}))

    if(centroid==TRUE){W<-as.data.frame(lapply(W,rep,times = ni))}

    sigmasq         <- NULL
    sigmasq.v       <- NULL
    beta0           <- NULL
    beta1           <- NULL
    loglik            <-NULL
    it              <-  1

    repeat{
        cat("Iteration sigma&beta",it,"\n",fill=T)
        sigmasq1            <- sigmasq0
        sigmasq1.v          <- sigmasq0.v
        Var                 <- vi.inv(n=ni, e=sigmasq0, u=sigmasq0.v, Z=Z)
        estGwr<-do.call(rbind,mapply(gwrEst, W,Pos=uniquePos,
                      MoreArgs=list(X=X,y=Y, Inv=Var$V.Inv, V=Var$V), SIMPLIFY=FALSE))
        beta.gwr<-estGwr[,c(1:dim(X)[2])]

        #with the Nelder-Mead procedure (ottimo) we can obtain the REML estimator of sigma2u and sigma 2 e with beta fixed.
        ottimo<-suppressWarnings(stats::optim(c(sigmasq0.v,sigmasq0),fn=logl,gr=grr, m=m, n=n, beta=beta.gwr, X=X, Y=Y, Z=Z,
                                  method="Nelder-Mead",
                                  control=list(fnscale=-1, trace = FALSE, maxit = 200) ))

        sigmasq0.v        <- ottimo$par[1]
        sigmasq0          <- ottimo$par[2]
        sigmasq.v         <- c(sigmasq.v, sigmasq0.v)
        sigmasq           <- c(sigmasq, sigmasq0)
        loglik            <- c(loglik, ottimo$value)

        erreur            <-abs(sigmasq0-sigmasq1)/abs(sigmasq1)+  abs(sigmasq0.v-sigmasq1.v)/abs(sigmasq1.v) #relative differenz
        #         erreur         <- abs(sigmasq0-sigmasq1) +  abs(sigmasq0.v-sigmasq1.v)
        #         erreur         <- (sigmasq0-sigmasq1)^2 +  (sigmasq0.v-sigmasq1.v)^2
        ifelse (it < maxit, ifelse(erreur > tol, {it=it+1}, {break}), {break})
    }

    vbeta   <-estGwr[,seq((dim(X)[2]+1),(dim(X)[2]*2))]
    proj    <-estGwr[,-seq(1,(dim(X)[2])*2)]
    colnames(beta.gwr)<-attr(X,"dimnames")[[2]]
    colnames(vbeta)<-attr(X,"dimnames")[[2]]
    LR<-as.numeric(-2*logLik(ex.lme, REML=F)+2*ottimo$value)
    df<-(2*sum(diag(proj))-sum(diag(t(proj)%*%proj))-dim(X)[2])
    p.val<-pchisq(LR,df=df, lower.tail = F)
    LR_test<-c("Likelihood" = ottimo$value ,"LR"=LR, "DF"= df, "P-Value" = p.val)



    conv           <- ifelse(it < maxit, 1, 0)
    n.iter         <- it

    result<-list(
            call = match.call(),
            Coefficients = beta.gwr,
            VarCoefficients = vbeta,
            Variance = list("residual"= sigmasq0, "raneff" =sigmasq0.v),
            LRtest = LR_test,
            bandwidth=band,
            nIter = it,
            Model = list(mf.s = mf.s,
                      formula=formula,
                      Projection = proj,
                      centroid = centroid)
            )

    class(result)<-"gwlmm"
    result
}
