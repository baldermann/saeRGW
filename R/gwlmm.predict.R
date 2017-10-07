#' Predict Method for (Robust) GWLMM Fits
#'
#' Interface for predicting small area means based on a
#' \code{\link{gwlmm}} ( or \code{\link{rgwlmm}}) model fit.
#'
#'
#'
#' @param object (formula) object of class \code{'gwlmm'} or \code{'rgwlmm'}
#' @param popdata (data.frame) an optional data frame (Default = NULL). If ommited,
#'  in-sample prdictions are estimated.
#' @param size (character) name for the column in \code{popdata}
#' containing the area-specific
#' population sizes (Default = NULL). Obligatory when \code{popAgg = TRUE}.
#' @param popAgg (logical) \code{TRUE} if popdata is aggregated and only contains
#' area level information
#' @param ... not used
#'
#'
#' @return
#' The function \code{predict.gwlmm} returns a list containing the predictions
#'
#' If \code{popdata = NULL}, a list with the following elements is returned.
#'
#' \itemize{
#' \item \code{call} (language) the call generating the value
#' \item \code{nIter} (numeric) number of iterations needed for the estimation of the
#' random effects (only for \code{rgwlmm}-objects)
#' \item \code{raneff} (numeric) named vector of random effects
#' \item \code{residuals} (numeric) vector of restimated residuals
#' \item \code{prediction} (numeric) vector of estimated in-sample predicted values
#' \item \code{xbeta} (numeric) vector of estimated fixed part of the in-sample predictions
#' }
#'
#' If \code{popdata} is a data.frame a list with the following elements is returned.
#'
#' \itemize{
#' \item \code{call} (language) the call generating the value
#' \item \code{nIter} (numeric) number of iterations needed for the estimation of the
#' random effects (only for \code{rgwlmm}-objects)
#' \item \code{areaEst} (data.frame) a data frame containing area-specific mean predictions and MSE estimates.
#' \item \code{ranEff} (data.farme) a data frame containing the random effects.
#' \item \code{sampleEst} (data.frame) a data frame containing the in-sample
#' predictions (prediction, residuals, xbeta)
#'}
#'
#'
#'@details
#'\itemize{
#'\item The argument \code{popdata} can have three different definitions:
#' (1) \code{popdata = NULL},
#' (2) \code{popdata} is a data frame for aggregated population information,
#' (3) \code{popdata} is a data frame for unit-level population information.
#' In case (1) only in-sample predictions are estimated. In case (2) the \code{size} argument is obligatory.
#' In case (3) \code{popAgg} must be \code{TRUE}. Population data can contain non-sampled areas.
#'
#'\item The method \code{predict} is implemented for three combinations of
#'geographic information in the sampled (S) and population (P)
#'data: (1) only centroid information in S and P; (2) unit-level geographic information in S and P;
#'(3)  unit-level geographic information in S and centriod informaton in P.
#'
#'\item For objects of class \code{gwlmm}: the MSE estimates for the small area means returned in \code{areaEst} are named \code{CCT} and \code{CCST}. The former is based
#'on pseudo-linearization following Chambers et al. (2011), and latter is linearization-based following
#'Chambers et al. (2014).
#'
#'\item For objects of class \code{rgwlmm}: the MSE estimates for the small area means returned in \code{areaEst}
#'are named \code{CCT} (and \code{CCT_bc})  and \code{CCST} ( and \code{CCST_bc}). The subscript \code{_bc} indicates the respective
#'MSE estimate for bias corrected robust area mean. The MSE estimates are based on the pseudo-linearizetion (CCT)
#'and the linearization-based approach in Chambers et al. (2014).
#'
#'\item For objects of class \code{rgwlmm}: the bias correction is based on the area-specific prediction error where the influence of extreme
#'values is restricted using Huber's influence function. \code{cbConst} defines this restriction but
#'should be less restrictive the \code{k} in \code{\link{rgwlmm}}. By default it is set \code{cbConst = 3}.
#'
#'}
#'@references
#'Chambers, R., J. Chandra, and N. Tzavidis (2011). On bias-robust mean
#'squared error estimation for pseudo-linear small area estimators. Survey
#'Methodology 37 (2), 153 - 170.
#'
#'Chambers, R., H. Chandra, N. Salvati, and N. Tzavidis (2014). Outlier
#'robust small area estimation. Journal of the Royal Statistical Society:
#'Series B 76 (1), 47 - 69.
#'
#'@examples
##'# Data sets ?sampleData, ?popaggData and ?popoutData are
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
#' @seealso \code{\link{rgwlmm}}, and \code{\link{gwlmm}}
#' @rdname predictFit
#' @export
predict.gwlmm<-function(object, popdata=NULL, size=NULL, popAgg  = TRUE, ...){

    X               <- model.matrix(Formula(object$Model$formula),object$Model$mf.s)
    Y               <- model.response(object$Model$mf.s)
    clusterid       <- droplevels(model.part(Formula(object$Model$formula),object$Model$mf.s,rhs=2)[,1])
    X.mean          <- aggregate(X, by=list(clusterid), FUN=mean, na.rm=TRUE)[,-1]
    geoInfo         <- model.part(Formula(object$Model$formula),object$Model$mf.s,rhs=3)
    ni              <- table(droplevels(clusterid))
    n               <- sum(ni)
    m               <- length(ni)
    Z               <- dummies::dummy(clusterid)
    colnames(Z)     <- sort(unique(clusterid))
    sigmasq         <- object$Variance$residual
    sigmasq.v       <- object$Variance$raneff
    Var             <- vi.inv(n=ni, e=sigmasq, u=sigmasq.v,Z=Z)

    # 1-projection -> produce residuals
    residuals       <- diag(1,n,n)-object$Model$Projection

    #Calculation of matrices needed
    B.1             <- diag(1/diag((t(Z)%*%Var$R.Inv%*%Z + Var$G.Inv)))
    B.2             <- t(Z)%*%Var$R.Inv
    B               <- B.1%*%B.2

    #Used later forsmall area prediction
    temp            <- B%*%residuals
    B.unshr         <- solve(t(Z)%*%Z)%*%t(Z)
    temp.un         <- B.unshr%*%residuals
    row.names(temp) <- row.names(t(Z))

    #Predictions:
    u.hat         <- temp%*%Y
    u.hat.un      <- temp.un%*%Y
    row.names(u.hat)<- row.names(t(Z))
    Zu            <- Z%*%u.hat
    Zu.un         <- Z%*%u.hat.un
    xbeta         <- c(object$Model$Projection %*%Y)
    pred.s        <- c(xbeta+Zu)
    pred.s.un     <- c(xbeta+Zu.un)
    res.s         <- Y-pred.s
    res.s.un      <- Y-pred.s.un

    if(is.null(popdata)){
     list(call = match.call(),
          raneff=u.hat,
          residuals=res.s,
          prediction = pred.s,
          xbeta = xbeta)
    }else{
################################################################################################
   #CCST: h1, Varaincane comming from the random Effects
   #derivate of H:
    p           <- dim(X)[2]
    dHdv        <- t(Z)%*%Var$R.Inv%*%Z  +Var$G.Inv
    Vt          <- sum(res.s^2)/(n-1)
    Vr          <- sum((res.s.un)^2)/(n-1)
    VarHv       <- Vt*t(Z)%*%Var$R.Inv%*%Var$R.Inv%*%Z
    Var.Cov.V   <- solve(dHdv)%*%VarHv%*%solve(dHdv)

    #CCST: h4 , Variance due to the estimtation of the variance components
    ZZt         <- Z%*%t(Z)
    #Estimate Variance of the variance parameters
    I1          <- 0.5*sum(diag(Var$V.Inv%*%ZZt%*%Var$V.Inv%*%ZZt))
    I2          <- 0.5*sum(diag(Var$V.Inv%*%Var$V.Inv%*%ZZt))
    I3          <- 0.5*sum(diag(Var$V.Inv%*%ZZt%*%Var$V.Inv))
    I4          <- 0.5*sum(diag(Var$V.Inv%*%Var$V.Inv))
    VC.s        <- solve(matrix(c(I1,I2,I3,I4),2,2))
    VC.u        <- Zu%*%t(Zu) + diag(sigmasq,n)
    dBds.v      <- B.1%*%Var$G.Inv%*%Var$G.Inv%*%B
    dBds.e      <- B.1%*%(t(Z)%*%Var$R.Inv%*%Var$R.Inv%*%Z)%*%B-
                   B.1%*%(t(Z)%*%Var$R.Inv%*%Var$R.Inv)
    V.sigma     <- dBds.v%*%VC.u%*%t(dBds.v)*VC.s[1,1] +dBds.e%*%VC.u%*%t(dBds.e)*VC.s[2,2] +
                    (dBds.v%*%VC.u%*%t(dBds.e) + dBds.e%*%VC.u%*%t(dBds.v))*VC.s[1,2]
################################################################################################
   #Small area prediction
    clusterid.r <- model.part(Formula(object$Model$formula),popdata,rhs=2)[,1]
    popdata     <- popdata[order(clusterid.r),]
    X.r         <- model.matrix(Formula(object$Model$formula),popdata)
    clusterid.r <- model.part(Formula(object$Model$formula),popdata,rhs=2)[,1]
    geoInfo.r   <- model.part(Formula(object$Model$formula),popdata,rhs=3)

    area_est    <- data.frame(area=sort(unique(clusterid.r)), est_sample=0, est_syn=0)
    MSE         <- list(area=sort(unique(clusterid.r)), v1=NULL, v1.un=NULL, bias =NULL,
                    bias.un=NULL,h1=NULL, h2=NULL,h3 = NULL,h4=NULL, h5 =NULL)
    if(object$Model$centroid == TRUE){centroid_s = TRUE}else{centroid_s = FALSE}
    if(centroid_s == TRUE){
          geoInfo<-aggregate(geoInfo, list(clusterid=clusterid), mean, na.rm=TRUE)
    }

    for (i in sort(unique(clusterid.r))){ #i=sort(unique(clusterid.r))[1]
            cat("Out of sample Estimation, area ",i,"\n",fill=T)
            i    <-as.character(i)
            thisr<-c(which(clusterid.r==i))
             if(popAgg !=TRUE){
                if(length(thisr)<2){next}
              }
            thiss<- c(which(clusterid==i))
              if(is.null(thiss)){thiss<-0}

            pos_u           <- which(row.names(u.hat)==i)
            thisu           <-ie(pos_u,temp[pos_u,])
            thisu.un        <-ie(pos_u,temp.un[pos_u,])
            pos             <-which(sort(unique(clusterid.r))==i)
            ir              <-rep(0,n)
            ir[thiss]       <-1

              if(popAgg ==TRUE){ #if(centroid_s==FALSE & centroid_r == TRUE)
                  Ri            <- popdata[,size][i] - length(Y[thiss])
                  thisCentroid  <- which(clusterid.r==i)
                  dist.cluster  <- spatstat::crossdist(geoInfo[,1],geoInfo[,2],geoInfo.r[thisCentroid,1],geoInfo.r[thisCentroid,2])
                  w.new         <- as.data.frame(lapply(as.data.frame(gwr.Gauss((dist.cluster)^2,object$bandwidth)),
                                   function(x){ifelse(x<=1e-06,min(x[x>1e-06]),x)}))
                    if(centroid_s==TRUE){
                            w.new   <- as.data.frame(rep(unlist(w.new),times = ni))
                            }
                  P<-mapply(pred.out,w=as.data.frame(w.new),Pos=thisr,SIMPLIFY=FALSE,
                               MoreArgs=list(Xout=X.r,X=X,Inv=Var$V.Inv, W1=diag(n),
                                 mXout=X.r[thisr,], V=Var$V))

                  P0        <-do.call(rbind, P)
                  PP        <-Ri*(P0[,c(1:n)] + thisu)
                  PP.un     <-Ri*(P0[,c(1:n)] + thisu.un)
                  PP1       <-PP
                  PP1.un    <-PP.un
                  welch     <-ir + PP
                  Ni        <-Ri + length(thiss)
              }else{
                  mXout<-apply(matrix(X.r[thisr,],nrow=length(thisr), ncol=dim(X)[2]),2,mean)
                  Ri            <- length(thisr)
                  dist.cluster  <- spatstat::crossdist(geoInfo[,1],geoInfo[,2],geoInfo.r[thisr,1],geoInfo.r[thisr,2])
                  thisr         < -unlist(thisr)
                  w.new         <- as.data.frame(lapply(as.data.frame(gwr.Gauss((dist.cluster)^2,object$bandwidth)),
                                   function(x){ifelse(x<=1e-06,min(x[x>1e-06]),x)}))

                  P<-mapply(pred.out,w=as.data.frame(w.new),Pos=thisr,SIMPLIFY=FALSE,
                                     MoreArgs=list(Xout=as.matrix(X.r,ncol=dim(X)[2]),
                                       X=X,Inv=Var$V.Inv, W1=diag(n), mXout=mXout, V=Var$V))

                  P0<-do.call(rbind, P)
                  PP<-t(apply(P0[,c(1:n)],1,'+',thisu))
                  PP.un <-t(apply(P0[,c(1:n)],1,'+',thisu.un))
                  PP1       <-apply(PP,2,sum)
                  PP1.un    <-apply(PP.un,2,sum)
                  welch     <-ir+as.vector(PP1)
                  Ni        <-Ri + length(thiss)
                }
         #Prediction
            area_est$est_sample[pos]    <- sum(welch*Y)/Ni
            area_est$est_syn[pos]       <- (sum(pred.s[thiss]) +sum(PP1*Y))/Ni
            est_syn.un      <-(ifelse(length(pred.s[thiss])==0,0,sum(pred.s.un[thiss]))+sum(PP1.un*Y))/Ni
            cost            <-(Ri/Ni)^2

        #MSE Estimation
            ####CCT :
            # CCT paper Gl. (7)
            MSE$v1[pos]         <-(1/Ni^2)*(sum((PP1^2 +(Ri/n))*res.s^2))
             # Variance with unshrunken mu_hat
            MSE$v1.un[pos]      <-(1/Ni^2)*(sum((PP1^2 +(Ri/n))*res.s.un^2))

            MSE$bias[pos]       <-(1/Ni)*sum(welch*pred.s)-area_est$est_syn[pos]
            MSE$bias.un[pos]    <-(1/Ni)*sum(welch*pred.s.un) - est_syn.un

            ####CCST:
            MSE$h1[pos]       <- cost*ie(pos_u,Var.Cov.V[pos_u,pos_u])
               #correct calculation
            if(popAgg ==TRUE){
                  MSE$h2[pos]         <- (Ri/Ni^2)*sum(P0[,-c(1:n)]*t(X.r[thisr,]))
            }else{
                  MSE$h2[pos]         <- (Ri/Ni^2)*sum(P0[,-c(1:n)]%*%mXout)
                 }
            MSE$h3[pos]       <- cost*sum(res.s.un^2)/((n-1)*(Ri))
            # h3 in the CCST-Paper
            MSE$h5[pos]       <- cost*ie(pos_u,V.sigma[pos_u,pos_u])

        } # Ends the prediction loop


     area_est$CCST   <-MSE$h1 + MSE$h2 +  MSE$h3 + MSE$bias.un^2 + MSE$h5
     area_est$CCT    <-MSE$v1.un + MSE$bias.un^2


      return(list(call = match.call(),
                  areaEst = area_est,
                  ranEff =  data.frame("area"= rownames(u.hat),
                            "u.hat"= c(u.hat)),
                  sampleEstt = data.frame("area"=clusterid,"prediction" = pred.s,
                             "residuals"= res.s,"xbeta" = xbeta)
                ))
  }
}
