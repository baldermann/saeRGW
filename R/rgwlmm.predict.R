## Outlier-robust estimation for the Geographically weighted linear mixed model (GWLMM)
#
# User interface for fitting a GWLMM model when the estimation process is distorted
# by extreme values (outliers).
#' @inheritParams predict.gwlmm
#' @param bcConst (numeric) needed for objects of class \code{rgwlmm}. Defines the tuning constant for influience function
#'  in the bias correction term (default = 3). Ses datails.
#' @param maxit (integer) needed for objects of class \code{rgwlmm}. Defines the maximum number of iterations for
#' the estiamtion of the random effects for \code{rgwlmm} objects (default = 100).
#' @param tol (numeric) needed for objects of class \code{rgwlmm}. Defines the tolerance for the convergence of the fitting process (default = 1e-04).
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
#'
#' @rdname predictFit
#' @export
predict.rgwlmm<-function(object, popdata = NULL, size=NULL, popAgg  = TRUE, bcConst = 3, maxit=100, tol=0.0001, ...){

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
    k               <- object$Model$tunConst
    const           <- 2 * pnorm(k) - 1 - 2 * k * dnorm(k) + 2 * (k^2) * (1-pnorm(k))
    sigmasq         <- object$Variance$residual
    sigmasq.v       <- object$Variance$raneff
    Var             <- vi.inv(n=ni, e=sigmasq, u=sigmasq.v,Z=Z)
    xbeta           <- c(object$Model$Projection %*%Y)
    eblup.v         <- c(Var$G %*% t(Z) %*% Var$V.Inv %*% (Y - xbeta))
    r               <- c(sqrt(Var$U.Inv) %*% (Y - xbeta))
    psi.r           <- hub.psi(r, k)$psi
    der.psi.r       <- hub.psi(r, k)$der.psi
    sum.v           <- NULL
    vv              <- NULL
    v0              <- eblup.v
    it              <- 1

    repeat
    {
        cat("Iteration area-effect",it,"\n",fill=T)
        J            <- c(rep(0,it-1), -1, 1)
        J            <- J[2:(it+1)]
        zv0            <- Z%*%v0
        r1             <- c(sqrt(Var$R.Inv) %*% (Y - xbeta - zv0))
        psi.rr          <- hub.psi(c(r1), k)
        psi.r1         <- psi.rr$psi
        der.psi.r1     <- psi.rr$der.psi
        v1             <- c(sqrt(Var$G.Inv) %*% v0)
        psi.v          <- hub.psi(v1, k)
        psi.v1         <- psi.v$psi
        der.psi.v1     <- psi.v$der.psi
        E.psi.der      <- pnorm(k) - pnorm(-k)

        H              <- t(Z) %*% sqrt(Var$R.Inv) %*%  psi.r1 - sqrt(Var$G.Inv) %*% psi.v1
        HH             <- t(Z) %*% sqrt(Var$R.Inv) %*% diag(rep(E.psi.der, n)) %*% sqrt(Var$R.Inv) %*% Z +
                        sqrt(Var$G.Inv) %*% diag(rep(E.psi.der, m)) %*% sqrt(Var$G.Inv)

        v0             <- c(solve(HH) %*% H)+v0
        vv             <- cbind(vv, v0)

        sum.v0         <- sum(abs(v0))
        sum.v          <- c(sum.v, sum.v0)

        erreur         <- abs(sum(J * sum.v))
        ifelse (it < maxit, ifelse(erreur > tol, {it=it+1}, {break}), {break})

    }

    conv.v  <-  ifelse(it < maxit, 1, 0)
    niter.v <-  it
#CCT: Weights for the CCT-MSE-Estimation:
    t1<-c(sqrt(Var$R.Inv) %*% (Y - xbeta - Z%*%v0))
    q1<-c(sqrt(Var$G.Inv) %*% v0)
    psi.t1 <-hub.psi(c(t1), k)
    psi.q1 <-hub.psi(c(q1), k)

    W2 <- diag(c(psi.t1$psi/t1))
    W3 <- diag(c(psi.q1$psi/q1))
    B.1    <- diag(1/diag((t(Z)%*%Var$R.Inv%*%W2%*%Z + Var$G.Inv%*%W3))) #is needed for h4!
    B.2    <-(t(Z)%*%Var$R.Inv%*%W2)  #is needed for h4!
    B      <-B.1%*%B.2
    B.unshr<-solve(t(Z)%*%Z)%*%t(Z)  # unshrunken weight
    residuals=diag(1,n,n)-object$Model$Projection
    temp<-B%*%residuals              #normal weights
    temp.un<-B.unshr%*%residuals     #unshrunken weights
#Predictions:
    u.hat       <-temp%*%Y
    u.hat.un    <-temp.un%*%Y

    row.names(u.hat)<-attr(Z,"dimnames")[[2]]
    Zu            <-Z%*%u.hat
    Zu.un         <-Z%*%u.hat.un
    pred.s        <-c(object$Model$Projection%*%Y+Zu)
    pred.s.un     <-c(object$Model$Projection%*%Y+Zu.un)
    res.s         <-Y-pred.s
    res.s.un      <-Y-pred.s.un

    if(is.null(popdata)){
        list( call = match.call(),
              nIter = niter.v,
              raneff=u.hat,
              residuals=res.s,
              prediction = pred.s,
              xbeta = xbeta)
    }else{
#scale parameter for the bias correction:
    phi_bc        <- median(abs(res.s))
################################################################################################
    #####CCST: h1 , Varaincane comming from the random Effects
    #derivative of H (robust ML erstimation equation) with respect to v:
    p           <- dim(X)[2]
    D1          <- diag(c(psi.t1$der.psi))
    D2          <- diag(c(psi.q1$der.psi))
    dHdv        <- t(Z)%*%sqrt(Var$R.Inv)%*%D1%*%sqrt(Var$R.Inv)%*%Z  +sqrt(Var$G.Inv)%*%sqrt(Var$G.Inv)
#     lambda1     <- 1+p/n*(var(der.psi.r)/(mean(der.psi.r)^2))
    lambda2     <- 1+p/n*(var(psi.t1$der.psi)/(mean(psi.t1$der.psi)^2))
#     vr          <- lambda1*sum(c(psi.r^2))/(n-p)
    Vt          <- lambda2*sum(c(psi.t1$psi)^2)/(n-p)
    VarHv       <- Vt*t(Z)%*%Var$R.Inv%*%Z
    Var.Cov.V   <- solve(dHdv)%*%VarHv%*%solve(dHdv)
    ####CCST: h4 , Variance due to the estimtation of the variance components
    ZZt         <- Z%*%t(Z)
    Q           <- t(psi.r)%*%Var$U.sqrt%*%Var$V.Inv
    Q2          <- Var$U.sqrt%*%Var$V.Inv
    Hsigma.v   <- (1/2)*(t(Z)%*%(psi.r*Q2%*%ZZt%*%t(Q2)%*%psi.r) - t(Z)%*%(const*diag(Var$V.Inv%*%ZZt)))#
    Hsigma.e   <- (1/2)*(t(Z)%*%(psi.r*Q2%*%t(Q2)%*%psi.r) - t(Z)%*%(const*diag(Var$V.Inv)))
#
    Hsigma      <- list(Hsigma.v,Hsigma.e)
    BB          <- matrix(mapply(function(a,b){(1/m)*sum(a*b)}, a=list(Hsigma[[1]],Hsigma[[1]],Hsigma[[2]],
                          Hsigma[[2]]),b=list(Hsigma[[1]],Hsigma[[2]],Hsigma[[1]],Hsigma[[2]])),2,2)
    D           <-hub.psi(r, k)$der.psi
    Q3          <-ZZt%*%Var$V.Inv
    Q4          <-Var$V.Inv%*%Var$V.Inv
    QQ1         <-t(t(Var$U.Inv%*%r*D)%*%Q2/(-2) + (1/2)*t(psi.r)%*%Var$U.Inv.sqrt%*%Var$V.Inv - t(psi.r)%*%Q2%*%Q3)*(Q3%*%Var$U.sqrt%*%psi.r)
    QQ2         <-t(t(Var$U.Inv%*%r*D)%*%Q2/(-2) + (1/2)*t(psi.r)%*%Var$U.Inv.sqrt%*%Var$V.Inv - t(psi.r)%*%Q2%*%Var$V.Inv)*t(Q2)%*%psi.r
    QQ3         <-t(t(Var$U.Inv%*%r*D)%*%Q2/(-2) + (1/2)*t(psi.r)%*%Var$U.Inv.sqrt%*%Var$V.Inv - t(psi.r)%*%Q2%*%Var$V.Inv)*(Q3%*%Var$U.sqrt%*%psi.r)

    C1          <-sum(t(Z)%*%(QQ1 + const*diag(Var$V.Inv%*%Q3%*%ZZt)/2)/m)
    C2          <-sum(t(Z)%*%(QQ2+const*diag(Q4)/2)/m)
    C3          <-sum(t(Z)%*%(QQ3 + const*diag(Var$V.Inv%*%Q3)/2)/m)
    C           <-matrix(c(C1,C3,C3,C2),2,2)
    #Estimated Variance of the variance parameters
    VC.s        <-solve(C)%*%BB%*%solve(C)/m
    VC.u        <-Zu%*%t(Zu) + diag(sigmasq,n)
    dBds.v      <-B.1%*%(Var$G.Inv%*%Var$G.Inv%*%W3 + Var$G.Inv%*%Var$G.Inv%*%diag(psi.q1$der.psi))%*%B/2
    dBds.e      <-B.1%*%(t(Z)%*%Var$R.Inv%*%Var$R.Inv%*%W2%*%Z + t(Z)%*%Var$R.Inv%*%Var$R.Inv%*%diag(psi.t1$der.psi)%*%Z)%*%B/2-
                  B.1%*%(t(Z)%*%Var$R.Inv%*%Var$R.Inv%*%W2 + t(Z)%*%Var$R.Inv%*%Var$R.Inv%*%diag(psi.t1$der.psi))/2
    V.sigma     <-dBds.v%*%VC.u%*%t(dBds.v)*VC.s[1,1] +dBds.e%*%VC.u%*%t(dBds.e)*VC.s[2,2] +
                    (dBds.v%*%VC.u%*%t(dBds.e) + dBds.e%*%VC.u%*%t(dBds.v))*VC.s[1,2]
#################################################################################
    #Small area prediction

    ##robustified dependent Variable
    psi.res.s         <-hub.psi(Y-xbeta,k*sqrt(sigmasq + sigmasq.v))
    Y.rob             <-xbeta + psi.res.s$psi
    res.rob           <-Y.rob-xbeta
    v.res.rob         <-(1 + p/n*var(psi.res.s$der.psi)/mean(psi.res.s$der.psi))*sum(res.rob^2)/(n-p)
    W1                <-diag(c(Y.rob/Y))
    clusterid.r       <- model.part(Formula(object$Model$formula),popdata,rhs=2)[,1]
    popdata           <- popdata[order(clusterid.r),]
    X.r               <- model.matrix(Formula(object$Model$formula),popdata)
    clusterid.r       <- model.part(Formula(object$Model$formula),popdata,rhs=2)[,1]
    geoInfo.r         <- model.part(Formula(object$Model$formula),popdata,rhs=3)
    area_est          <-data.frame(area=sort(unique(clusterid.r)), est_sample=0, est_syn=0,bc2=0)
    MSE               <-list(area=sort(unique(clusterid.r)), v1=NULL,v1.un=NULL,
                            bias =NULL,bias.un=NULL,h1=NULL, h2=NULL,h3 = NULL,h4=NULL, h5 =NULL)
    if(object$Model$centroid == TRUE){centroid_s = TRUE}else{centroid_s = FALSE}
    if(centroid_s==TRUE){
          geoInfo<-aggregate(geoInfo, list(clusterid=clusterid), mean, na.rm=TRUE)
    }

    for (i in sort(unique(clusterid.r))){ # i = sort(unique(clusterid.r))[1]
            cat("Out of sample Estimation, area ",i,"\n",fill=T)
            i     <- as.character(i)
            thisr <- c(which(clusterid.r==i))
            if(popAgg != TRUE & length(thisr) < 2 ){next}

            thiss <- c(which(clusterid==i))
            if(length(thiss)[1] == 0){thiss <- 0}


            pos_u           <- which(row.names(u.hat)==i)
            thisu           <- ie(pos_u,temp[pos_u,])
            thisu.un        <- ie(pos_u,temp.un[pos_u,])
            pos             <- which(sort(unique(clusterid.r))==i)
            ir              <- rep(0,n)
            ir[thiss]       <- 1

            if(popAgg ==TRUE){
                  Ri            <- popdata[,size][i] - length(Y[thiss])
                  thisCentroid  <- which(clusterid.r==i)
                  dist.cluster  <- spatstat::crossdist(geoInfo[,1],geoInfo[,2],
                                                      geoInfo.r[thisCentroid,1],geoInfo.r[thisCentroid,2])
                  w.new         <- as.data.frame(lapply(as.data.frame(gwr.Gauss((dist.cluster)^2,object$bandwidth)),
                                                        function(x){ifelse(x<=1e-06,min(x[x>1e-06]),x)}))
                  if(centroid_s!=TRUE & length(thiss)==1){gi <- 1}
                  if(centroid_s!=TRUE & length(thiss)>1){ gi <- ie(pos_u,apply((object$Model$Weights[thiss,thiss]),1,mean))}
                  if(centroid_s==TRUE){
                        w.new   <- as.data.frame(rep(unlist(w.new),times = ni))
                        gi      <- ie(pos_u,rep(1,length(thiss)))
                  }

                  P         <- mapply(pred.out,w=as.data.frame(w.new),Pos=thisr,SIMPLIFY=FALSE,
                                      MoreArgs=list(Xout=X.r,X=X,Inv=Var$V.Inv, W1=W1,
                                        mXout=X.r[thisr,], V=Var$V))



                  P0        <-do.call(rbind, P)
                  PP        <-Ri*(P0[,c(1:n)] + thisu)
                  PP.un     <-Ri*(P0[,c(1:n)] + thisu.un)
                  PP1       <-PP
                  PP1.un    <-PP.un
                  welch     <-ir + PP
                  Ni        <-Ri + length(thiss)

            }else{
                  mXout         <- apply(matrix(X.r[thisr,],nrow=length(thisr), ncol=dim(X)[2]),2,mean)
                  Ri            <- length(thisr)
                  dist.cluster  <- spatstat::crossdist(geoInfo[,1],geoInfo[,2],geoInfo.r[thisr,1],geoInfo.r[thisr,2])
                  thisr         <- unlist(thisr)
                  w.new         <- as.data.frame(lapply(as.data.frame(gwr.Gauss((dist.cluster)^2,
                                    object$bandwidth)),function(x){ifelse(x<=1e-06,min(x[x>1e-06]),x)}))
                  if(length(thiss)>1){
                      gi        <- ie(pos_u,apply((object$Model$Weights[thiss,thiss]),1,mean))
                      }else{gi  <- ie(pos_u,object$Model$Weights[thiss,thiss])}

                  P         <- mapply(pred.out,w=as.data.frame(w.new),Pos=thisr,SIMPLIFY=FALSE,
                                 MoreArgs=list(Xout=as.matrix(X.r,ncol=dim(X)[2]),X=X,Inv=Var$V.Inv, W1=W1, mXout=mXout, V=Var$V))
                  P0        <- do.call(rbind, P)
                  PP        <- t(apply(P0[,c(1:n)],1,'+',thisu))
                  PP.un     <- t(apply(P0[,c(1:n)],1,'+',thisu.un))
                  PP1       <- apply(PP,2,sum)
                  PP1.un    <- apply(PP.un,2,sum)
                  welch     <- ir+as.vector(PP1)
                  Ni        <- Ri + length(thiss)
            }

            area_est$est_sample[pos]    <- sum(welch*Y)/Ni
            area_est$est_syn[pos]       <- (sum(pred.s[thiss]) +sum(PP1*Y))/Ni
            area_est$bc2[pos]           <- ie(pos_u,sum(gi*phi_bc*hub.psi(res.s[thiss]/phi_bc,bcConst)$psi)/sum(gi))*Ri/Ni
            est_syn.un                  <-(ifelse(length(pred.s[thiss])==0,0,sum(pred.s.un[thiss]))+sum(PP1.un*Y))/Ni
            cost                        <-(Ri/Ni)^2

            ####CCT :
            MSE$v1[pos]         <-(1/Ni^2)*(sum((PP1^2 +(Ri/n))*res.s^2))
            MSE$v1.un[pos]      <-(1/Ni^2)*(sum((PP1^2 +(Ri/n))*res.s.un^2))

            MSE$bias[pos]       <-(1/Ni)*sum(welch*pred.s)-area_est$est_syn[pos]
            MSE$bias.un[pos]    <-(1/Ni)*sum(welch*pred.s.un) - est_syn.un

            ### CCT: adjustment for weights with bias correction:
            W4      <- hub.psi(res.s/phi_bc,bcConst)$psi/(res.s/phi_bc)
            ir_bc2  <- ie(pos_u,ir*(W4/sum(gi)))
            BB_bc2  <- ie(pos_u,t(gi*W4[thiss])%*%object$Model$Projection[thiss,]/sum(gi))
            thisu_bc2<-ie(pos_u,apply(gi*W4[thiss]%*%t(thisu),2,sum)/sum(gi))
            if(popAgg ==TRUE){
                PP1_bc2  <- PP + Ri*(ir_bc2-BB_bc2-thisu_bc2)
            }else{
                PP1_bc2  <- apply(t(apply(PP,1,'+',(ir_bc2-BB_bc2-thisu_bc2))),2,sum)
            }
            MSE$CCT_bc2[pos]<-(1/Ni^2)*(sum((PP1_bc2^2 +(Ri/n))*res.s.un^2))

            ####CCST:
            MSE$h1[pos]       <- cost*ie(pos_u,Var.Cov.V[pos_u,pos_u])
            if(popAgg ==TRUE){
                  MSE$h22bc[pos]      <-cost*ie(pos_u,sum((P0[-c(1:n)]-apply(matrix(object$Model$vbeta[thiss,],
                                                nrow=length(thiss)),2,mean))*(t(X.r[thisr,]) - X.mean[pos_u,])))
                  MSE$h2[pos]         <- (Ri/Ni^2)*sum(P0[,-c(1:n)]*t(X.r[thisr,]))
            }else{
                  MSE$h2[pos]         <- (Ri/Ni^2)*sum(P0[,-c(1:n)]%*%mXout)
                  MSE$h22bc[pos]      <- cost*ie(pos_u,sum((apply(P0[,-c(1:n)],2,mean)
                                              -apply(matrix(object$Model$vbeta[thiss,],nrow=length(thiss)),2,mean))*(t(mXout) - X.mean[pos_u,])))
               }
            MSE$h3[pos]       <- cost*sum(res.s.un^2)/((n-1)*(Ri))
            MSE$h4[pos]       <- MSE$bias.un[pos]
            MSE$h5[pos]       <- cost*ie(pos_u,V.sigma[pos_u,pos_u])  # h3 in the CCST-Paper
            if(length(thiss)>1){
            MSE$h7_bc2[pos]       <-ie(pos_u,sum(gi*(phi_bc*hub.psi(res.s.un[thiss]/phi_bc, bcConst)$psi)^2)/(sum(gi)*(sum(gi)-(sum(gi^2)/sum(gi)))))*cost
            }else{MSE$h7_bc2[pos]=0}
        }



    area_est$est_bc     <-area_est$est_sample + area_est$bc2
    area_est$CCST       <-MSE$h1 + MSE$h2 +  MSE$h3 + MSE$h4^2 + MSE$h5
    area_est$CCST_bc    <-MSE$h22bc +  MSE$h3 +  MSE$h7_bc2
    area_est$CCT        <-MSE$v1.un + MSE$bias.un^2
    area_est$CCT_bc     <-MSE$CCT_bc2


    return(list(call = match.call(),
                nIter    = niter.v,
                area_est = area_est,
                rand.eff =  data.frame("area"= rownames(u.hat), "u.hat"= c(u.hat)),
                sample_est=data.frame("area"=clusterid,"prediction" = pred.s,
                                      "residuals"= res.s, "xbeta" = xbeta)
                ))

    }
}


