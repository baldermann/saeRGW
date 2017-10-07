

###########################
#Functions needed
###########################
# library(abind)
#Clean simulation results
clean_res<-function(Results){

    f1<-function(x,y){
    sapply(y,function(x,y){(x[y])}, y=x)
    }
    f2 <-function(x,y){
    z<-lapply(x,f1,y=y)
    names(z)<-x
    z
    }
    f3<-function(x,y,z){
    l<-mapply(abind, x[[z]],along=y[[z]], MoreArgs = list(use.first.dimnames=T))
   }
    ldim<-function(x){
        sapply(x,function(x){length(dim(x))+1})
    }


#     names(Results)<-c(1:simruns)
    zz<-names(Results[[1]])
    zz.names<-lapply(Results[[1]], names)
    zz.dims<-sapply(Results[[1]], ldim)
    Results<-lapply(zz,function(x,y){lapply(y, function(x,y){x[[y]]}, y=x)}, y=Results)
    names(Results)=zz
    Results<-mapply(f2,x=zz.names,y=Results)
    Results<-lapply(zz,f3,x=Results,y=zz.dims)
    names(Results)=zz
    Results
}

clean_res2<-function(Results,problem=NULL, newNames=NULL){

    f1<-function(x,y){
    sapply(y,function(x,y){(x[y])}, y=x)
    }
    f2 <-function(x,y){
    z<-lapply(x,f1,y=y)
    names(z)<-x
    z
    }
    f3<-function(x,y,z){
    l<-mapply(abind, x[[z]],along=y[[z]], MoreArgs = list(use.first.dimnames=T))
   }
    ldim<-function(x){
        sapply(x,function(x){length(dim(x))+1})
    }

    f4<-function(x,names){
    if(is.list(x)==T){
            names(x)<-names
            as.data.frame(x)
            }else{as.data.frame(x)}
    }

    f5<-function(x, names){
    lapply(x,f4,names=names)
    }





#     names(Results)<-c(1:simruns)
    zz<-names(Results[[1]])
    zz.names<-lapply(Results[[1]], names)
    zz.dims<-sapply(Results[[1]], ldim)
    Results1<-lapply(zz,function(x,y){lapply(y, function(x,y){x[[y]]}, y=x)}, y=Results)
    names(Results1)=zz
    if(length(problem)!=0){
        Results1[[problem]]<-lapply(Results1[[problem]],f5,names=newNames)
        zz.dims[[problem]]<-ldim(Results1[[problem]][[1]])
    }
    Results2<-mapply(f2,x=zz.names,y=Results1)
    Results3<-lapply(zz,f3,x=Results2,y=zz.dims)
    names(Results3)=zz
    Results3
}





RMSEfun<-function(Estimates, areas, valid, fun, RMSE_True){

    enames<-names(Estimates)

    if(fun=="ERMSE"){
        ff1<-function(x,areas,val,y,RMSE_True=NULL, est){
            kk<-apply(sqrt(x[[est]][areas,y,val[[est]]]),1,mean,na.rm = T)
            names(kk)<-areas
            return(kk)
    }
    }
    if(fun=="EMSE"){
        ff1<-function(x,areas,val,y,RMSE_True=NULL, est){
            kk<-apply(x[[est]][areas,y,val[[est]]],1,mean,na.rm = T)
            names(kk)<-areas
            return(kk)
        }
    }
    if(fun=="RB"){
        ff1<-function(x,areas,val,y,RMSE_True, est){
            rmse_true<-matrix(rep(RMSE_True[,est],length(val[[est]])),ncol = length(val[[est]]))
            kk<-100*apply((sqrt(x[[est]][areas,y,val[[est]]]) - rmse_true)/rmse_true,1,mean,na.rm = T)
            names(kk)<-areas
            return(kk)
        }
    }
    if(fun=="RRMSE"){

        ff1<-function(x,areas,val,y,RMSE_True, est){
            rmse_true<-matrix(rep(RMSE_True[,est],length(val[[est]])),ncol = length(val[[est]]))
            kk<-100*sqrt(apply(((sqrt(x[[est]][areas,y,val[[est]]]) - rmse_true)/rmse_true)^2,1,mean,na.rm = T))
            names(kk)<-areas
            return(kk)
        }

    }
    ff2<-function(x,areas,val,est,RMSE_True){
        ll<-lapply(X=unlist(attr(x[[est]],"dimnames")[2]),FUN= ff1, x=x, areas=areas,val=val, est=est, RMSE_True=RMSE_True)
        names(ll)<-unlist(attr(x[[est]],"dimnames")[2])
        return(ll)
    }
    ll<-lapply(X=enames,FUN=ff2, x=Estimates, val = valid, areas=areas, RMSE_True=RMSE_True)
    names(ll)<-enames
    return(ll)

}

median.ordered <- function(x,...)
{
    levs <- levels(x)
    m <- median(as.integer(x))
    if(floor(m) != m)
    {
      # warning("Median is between two values; using the first one")
      m <- floor(m)
    }
    ordered(m, labels = levs, levels = seq_along(levs))
}




logl=function(est,m, n, beta, Y, Z,X){
#      browser()
    sigma.v<-est[1]
    sigma.e<-est[2]
    I=diag(1,m)
    V<-sigma.e*diag(1,n)+sigma.v*Z%*%I%*%t(Z) #Varianzmatrix
    #         Vi<-vi.inv.NoSp(n=ni, e=sigma.e, u=sigma.v)
    Vi<-solve(V)                              #Inverse der Varianzmatrix von Y
    b.s<-beta
    xbeta<-(X*b.s)%*%rep(1,dim(X)[2])
    y<-Y
    ee=eigen(V) #Eigenwerte der Varianzmatrix - zur einfachen Berechnung der Determinate
    -(n/2)*log(2*pi)-0.5*sum(log(ee$value))-(0.5)*t(y-xbeta)%*%Vi%*%(y-xbeta) #Log-Likelihood
}


Rlogl=function(est,m, n, beta, Y, Z,X){
          browser()
    sigma.v<-est[1]
    sigma.e<-est[2]
    I=diag(1,m)
    V<-sigma.e*diag(1,n)+sigma.v*Z%*%I%*%t(Z) #Varianzmatrix
    #         Vi<-vi.inv.NoSp(n=ni, e=sigma.e, u=sigma.v)
    Vi<-solve(V)
    #Inverse der Varianzmatrix von Y
    b.s<-beta
    xbeta<-(X*b.s)%*%rep(1,dim(X)[2])
    y<-Y
    ee=eigen(V) #Eigenwerte der Varianzmatrix - zur einfachen Berechnung der Determinate
    -((n-dim(X)[2])/2)*log(2*pi)-0.5*sum(log(ee$value))-(0.5)*t(y-xbeta)%*%Vi%*%(y-xbeta) - (0.5)*sum(log(eigen(t(X)%*%Vi%*%X)$value))#Log-Likelihood
}



#this function contains the first derivatives of the log-likelihood,s (that is gradient)

grr=function(est,m, n, beta, Y, Z,X){
    #         browser()
    sigma.v<-est[1]
    sigma.e<-est[2]

    s<-matrix(0,2,1)              # vector of the partial deriatives of log-likelihood with respcet to  rhospat,sigma.v,and sigma.e
    I<-diag(1,m)
    ZItZ<-Z%*%t(Z)
    V<-sigma.e*diag(1,n)+sigma.v*ZItZ #Varianzmatrix, jedoch ohne Räumlich Struktur - W-Matrix aus paper fehlt
    Vi<-solve(V)   #Inverse der Varianz
    Q <-Vi%*%ZItZ
    b.s<-beta
    xbeta<-(X*b.s)%*%rep(1,dim(X)[2])
    y<-Y
    #Ableitung nach Varianzparameter der Area-spezifischen random Effects
    s[1,1]<--sum(diag(Q))+((t(y-xbeta)%*%(Q%*%Vi)%*%(y-xbeta)))
    ##Ableitung nach Varianzparameter der individuellen Residuen
    s[2,1]<--sum(diag(Vi))+((t(y-xbeta)%*%(Vi%*%Vi)%*%(y-xbeta)))
    c(s[1,1],s[2,1])
}


Rgrr=function(est,m, n, beta, Y, Z,X){
#             browser()
    sigma.v<-est[1]
    sigma.e<-est[2]

    s<-matrix(0,2,1)              # vector of the partial deriatives of log-likelihood with respcet to  rhospat,sigma.v,and sigma.e
    I<-diag(1,m)
    ZItZ<-Z%*%t(Z)
    V<-sigma.e*diag(1,n)+sigma.v*ZItZ #Varianzmatrix, jedoch ohne Räumlich Struktur - W-Matrix aus paper fehlt
    Vi<-solve(V)   #Inverse der Varianz
    Q <-Vi%*%ZItZ
    b.s<-beta
    xbeta<-(X*b.s)%*%rep(1,dim(X)[2])
    y<-Y
    XvXi<-solve(t(X)%*%Vi%*%X)
    #Ableitung nach Varianzparameter der Area-spezifischen random Effects
    s[1,1]<--sum(diag(Q))+((t(y-xbeta)%*%(Q%*%Vi)%*%(y-xbeta))) + sum(diag(XvXi%*%t(X)%*%(Q%*%Vi)%*%X))
    ##Ableitung nach Varianzparameter der individuellen Residuen
    s[2,1]<--sum(diag(Vi))+((t(y-xbeta)%*%(Vi%*%Vi)%*%(y-xbeta)))  + sum(diag(XvXi%*%t(X)%*%(Vi%*%Vi)%*%X))
    c(s[1,1],s[2,1])
}


#Allgemeine Inverse einer Blockdiagonalen mit variierenden stichprobenumfängen in den Areas,
#sigma_u, sigma_e konstant
vi.inv<-function(n,e,u,Z){

    V.Inv<-as.matrix(bdiag(lapply(n, function(n) solve(diag(n)*e+matrix(u,n,n)))))
    m<-length(n)
    G              <- diag(rep(u, m))
    G.Inv          <- diag(1/rep(u, m))
    R.Inv          <- diag(1/rep(e, sum(n)))
    U.Inv          <- diag(1/rep(e+u,sum(n)))
    U.sqrt         <- diag(sqrt(rep(e+u,sum(n))))
    U.Inv.sqrt     <- diag(sqrt(diag(U.Inv)))
    V              <- Z%*%G%*%t(Z) + diag(e,sum(n))
    U              <-sqrt(V)
    list(V=V, V.Inv=V.Inv, G=G, G.Inv=G.Inv, R.Inv=R.Inv, U.Inv=U.Inv,
         U.sqrt= U.sqrt, U.Inv.sqrt=U.Inv.sqrt, U=U)
}


vi.inv.NoSp<-function(n,e,u){

    as.matrix(bdiag(lapply(n, function(n) solve(diag(n)*e+matrix(u,n,n)))))

}
#Allgemeine Inverse einer Blockdiagonalen mit variierenden stichprobenumfängen in den Areas und gewichtung,
#sigma_u, sigma_e konstant


vi.inv.NoSp_w <-function(n, e,u,w){
    e<-e
    u<-u
    w<-partition.vector(w,n)
    as.matrix(bdiag(mapply(function(n,w)solve(diag(1/w)*e+matrix(u,n,n)), n,w, SIMPLIFY=F)))
}


#Schätzung der lokalen Betas und der lokalen Projektionsvektoren fär die insample-Units. (nicht robust)
#Startwerte fär die approximation der robusten Schätzung
gwrEst<-function(w,Pos,X,y,Inv,V){
    w.i<-diag(w)
    Q<-t(X)%*%w.i%*%Inv
    t<-solve(Q%*%X)%*%Q
    v.beta<-diag(t%*%V%*%t(t))
    cbind(b=matrix(rep(t%*%y,each=length(Pos) ),nrow=length(Pos)),v.beta =matrix(rep(v.beta,each=length(Pos) ),nrow=length(Pos)), r=X[Pos,]%*%t)
}

# w = W[[1]]
# Pos = uniquePos[[1]]
# y=Y
# VInv=Var$V.Inv
# V=Var$V
# UInv=Var$U.Inv
# U=diag(diag(V))
# beta = beta.gwr


gwrRobustEst<-function(w,beta,X,y,VInv,V,UInv,k,U,Pos){
#           browser()
    nonzero         <-which(1/w!=Inf)
    w.inv           <-diag(1/w[nonzero])
    w               <-diag(w[nonzero])
    sqrtWinvU       <-diag(sqrt(diag(w.inv))*sqrt(diag(U[nonzero,nonzero])))
    Q               <- t(X[nonzero,])%*%w%*%VInv[nonzero,nonzero]%*%sqrtWinvU
    r               <- diag(1/diag(sqrtWinvU))%*%(y[nonzero]-X[nonzero,]%*%beta[Pos[1],])
    psi.fn          <- hub.psi(r, k)
    W1              <-diag(c(psi.fn$psi/r))
    A               <-solve(Q%*%W1%*%diag(1/diag(sqrtWinvU))%*%X[nonzero,])%*%Q%*%W1%*%diag(1/diag(sqrtWinvU))
    b               <-matrix(matrix(A%*%y[nonzero],1,dim(beta)[2])[rep(1,length(Pos)),],nrow=length(Pos))
    proj            <-X[Pos,]%*%A
    v.beta          <- A%*%V%*%t(A)
    cbind(b,
        X[Pos,]%*%v.beta,
        matrix(matrix(diag(v.beta),1,dim(beta)[2])[rep(1,length(Pos)),],nrow=length(Pos)),
        proj)
   }


#Huber psi

hub.psi <- function(x, k)
{
    psi <- ifelse(abs(x) <= k, x, sign(x) * k)
    der.psi <- ifelse(abs(x) <= k, 1, 0)
    list(psi = psi, der.psi = der.psi)
}


my.psi<-function(u,c){

    sm<-median(abs(u))/0.6745
    w <- psi.huber(u/sm,c)
    w*u
}





#Distanzmatritzen und Predictions für Out-Of-sample Units

pred.out.r<- function(w,Pos,X,Xout,y,Inv,u){
#     browser()
    w.i<-diag(w)
    Q<-t(X)%*%w.i%*%Inv
    (cbind(1,Xout[Pos])%*%solve(Q%*%X)%*%Q)%*%y +u
}


pred.out<-function(w,Pos,X,Xout,Inv,W1, mXout,V){
             # browser()
    x.i <-Xout[Pos,]

    w.i<-diag(w)
    Q<-t(X)%*%w.i%*%Inv #pxn
    QQ<-solve(Q%*%X)%*%Q
    cbind(x.i%*%QQ%*%W1, #1xn vector
        x.i%*%QQ%*%V%*%t(QQ)) #1Xp vector
    }



where<-function(x,y){which(x==y)}

ie<-function(x,y){
    if(length(x)!=0){
        y
    }else{0}
}


####################################################################################################
#################Convergenzalgorithmen


#Netwonverfahren (Sinha und Rao, 2008) und fixpunkt

ConvGWRobust<-function(method,beta.gwr0,sigmasq0.v,sigmasq0, ni, n,m, tol,maxit,k,X,W,Y, Z, Cpp,Pos){
    #if(method==NULL){stop('no method specified')}
    mehtodsPossible<-c("Newton","Fix")

#       browser()

    if(length(which(mehtodsPossible == method))==0) stop("method:", method ,"\n not defined or possible")
    sigmasq        <- NULL
    sigmasq.v      <- NULL
    beta0            <- NULL
    beta1            <- NULL

    const           <- 2 * pnorm(k) - 1 - 2 * k * dnorm(k) + 2 * (k^2) * (1-pnorm(k))
    K               <- diag(rep(const,n))
    it              <-  1
    repeat{
#           if(it==4){browser()}
#   browser()
        cat("Iteration sigma&beta",it,"\n",fill=T)
        #Taylor-Approxomation der Varianzparameter
        beta.gwr1           <- beta.gwr0
        sigmasq1            <- sigmasq0
        sigmasq1.v          <- sigmasq0.v
        #V                  <- diag(sigmasq0,n) + Z%*%diag(sigmasq0.v,m)%*%t(Z)
        Var                 <-vi.inv(n=ni, e=sigmasq0, u=sigmasq0.v, Z=Z)


        #Schätzung von lokalen Betas_ij

# # browser()


# browser()
# Pos=uniquePos
# w=W
# X=X
# y=Y
# VInv=Var$V.Inv
# UInv=Var$U.Inv
# k=k
# U=Var$U
# beta = beta.gwr0
# V=Var$V



if(Cpp == TRUE){
        k_val        <- pnorm(k) - pnorm(-k)
time1<-system.time(beta.gwr0    <-do.call(rbind,mapply(gwrRobustEstCpp,w=W,Pos=Pos,SIMPLIFY=FALSE,
                                        MoreArgs=list(X=X, y=Y, VInv=Var$V.Inv,
                                                      UInv=Var$U.Inv, k=k, U=Var$U, k_val = k_val, beta = beta.gwr0, V=Var$V))))
}else{
time2<-system.time(beta.gwr0            <-do.call(rbind,mapply(gwrRobustEst,w=W,Pos=Pos,SIMPLIFY=FALSE,
                                        MoreArgs=list(X=X, y=Y, VInv=Var$V.Inv,
                                                     UInv=Var$U.Inv, k=k, U=Var$U,beta = beta.gwr0, V=Var$V))))

}




        vbeta     <-beta.gwr0[,seq((dim(X)[2]+1),(dim(X)[2]*2))]
        vbeta2    <-beta.gwr0[,seq((dim(X)[2]*2+1),(dim(X)[2]*3))]
        proj      <-beta.gwr0[,-seq(1,(dim(X)[2])*3)]
        beta.gwr0 <-beta.gwr0[,c(1:dim(X)[2])]
        xbeta     <- (X*beta.gwr0)%*%rep(1,dim(X)[2])
        yhat      <- proj%*%Y



# Schätzung der Varianzparameter
        r               <- c(sqrt(Var$U.Inv) %*% (Y - xbeta))
        psi.fn          <- hub.psi(r, k)
        psi.r           <- psi.fn$psi
        der.psi.r       <- psi.fn$der.psi
        E.psi.der       <- pnorm(k) - pnorm(-k)


        qq              <- Var$V.Inv %*% Var$U.sqrt%*% psi.r
        ZZ              <- Z %*% t(Z) # Ableitung von V nach sigma_v

        if(method=="Newton"){
            q2             <- t(qq) %*% qq - sum(diag(Var$V.Inv %*% K))  # Vergleiche (3.74) mit Ableitung nach sigma_e
            dpsi.dsig      <- (-1/2) * der.psi.r * diag(Var$U.Inv) * r
            d1             <- t(dpsi.dsig) %*% Var$U.sqrt%*% Var$V.Inv
            d2             <- (1/2) * t(psi.r) %*% Var$U.Inv.sqrt%*% Var$V.Inv
            d3             <-  t(psi.r) %*% Var$U.sqrt%*% (Var$V.Inv%*%Var$V.Inv)
            der.qq         <- d1 + d2 - d3 # Komponenten aus der Mitte der Seite 53 nach "where ...."

            M2             <- c(2 * der.qq %*% qq +  sum(diag((Var$V.Inv%*%Var$V.Inv) %*% K))) # Seite 53 oben mit \partial V / \partial sigma_e=Id

            sigmasq0       <- sigmasq0 - c(solve(M2) %*% q2) # Vergleiche (3.76)
            sigmasq0       <- ifelse(sigmasq0 < 0, -sigmasq0/100, sigmasq0)
            sigmasq0       <- ifelse(sigmasq0 > 1000, 1, sigmasq0)
            sigmasq        <- c(sigmasq, sigmasq0)  # Vektor mit iterativen Sigma_e

            # Schätzung von Sigma_v

            dpsi.dsig.v    <- (-1/2) * der.psi.r * diag(Var$U.Inv) * r

            ZZ             <- Z %*% t(Z) # Ableitung von V nach sigma_v
            VdelV          <- Var$V.Inv %*% ZZ
            q3             <- t(qq) %*% ZZ %*% qq - sum(diag(VdelV %*% K)) # Vergleiche (3.74) mit Ableitung nach sigma_v

            d1.v           <- t(dpsi.dsig.v) %*% Var$U.sqrt%*% Var$V.Inv
            d2.v           <- (1/2) * t(psi.r) %*% sqrt(Var$U.Inv) %*% Var$V.Inv
            d3.v           <-  t(psi.r) %*% Var$U.sqrt%*% (VdelV %*% Var$V.Inv)
            der.qq.v       <- d1.v + d2.v - d3.v # Komponenten aus der Mitte der Seite 53 nach "where ...."

            M3             <- c(2 * der.qq.v %*% ZZ %*% qq + sum(diag(VdelV %*% VdelV %*% K))) # Seite 53 oben mit Ableitung nach sigma_v

            sigmasq0.v     <- sigmasq0.v - c(solve(M3) %*% q3) # Vergleiche (3.76)
            sigmasq0.v     <- ifelse(sigmasq0.v < 0, -sigmasq0.v/100, sigmasq0.v)
            sigmasq0.v     <- ifelse(sigmasq0.v > 1000, sigmasq0.v/100, sigmasq0.v)
            sigmasq.v      <- c(sigmasq.v, sigmasq0.v) # Vektor mit iterativen Sigma_v
        }

        if(method=="Fix"){

#             KK5<-sqrt(U) %*% V.inv
#             KK0<-t(psi.r)%*% KK5
            a1             <- t(qq)%*%qq
            a2             <- t(qq)%*%ZZ%*%qq
            a              <- matrix(c(a1[1],a2[1]), nrow = 2, ncol=1)

            KK2     <-K%*%Var$V.Inv%*%Var$V.Inv
            KK3     <-K%*%Var$V.Inv%*%ZZ%*%Var$V.Inv

            A1      <-sum(diag(KK2))
            A2      <-sum(diag(KK2%*%ZZ))
            A3      <-sum(diag(KK3))
            A4      <-sum(diag(KK3%*%ZZ))
            A       <-matrix(c(A1,A3,A2,A4), nrow = 2, ncol=2)

            sigma           <-solve(A)%*%a
            sigmasq0        <-sigma[1,]
            sigmasq0       <- ifelse(sigmasq0 < 0, -sigmasq0/100, sigmasq0)
#             sigmasq0       <- ifelse(sigmasq0 > 1000, 1, sigmasq0)
            sigmasq         <- c(sigmasq, sigmasq0)
            sigmasq0.v      <-sigma[2,]
            sigmasq0.v     <- ifelse(sigmasq0.v < 0, -sigmasq0.v/100, sigmasq0.v)
#             sigmasq0.v     <- ifelse(sigmasq0.v > 1000, sigmasq0.v/100, sigmasq0.v)
            sigmasq.v       <- c(sigmasq.v, sigmasq0.v)
        }

          erreur            <-abs(sigmasq0-sigmasq1)/sigmasq1+  abs(sigmasq0.v-sigmasq1.v)/sigmasq1.v #relative differenz
#         erreur         <- abs(sigmasq0-sigmasq1) +  abs(sigmasq0.v-sigmasq1.v)
#         erreur         <- (sigmasq0-sigmasq1)^2 +  (sigmasq0.v-sigmasq1.v)^2
#         erreur        <- (sigmasq0-sigmasq1)^2 +  (sigmasq0.v-sigmasq1.v)^2 # Stoppkriterium #abs(mean(apply(beta.gwr0-beta.gwr1,1,sum))) +
        ifelse (it < maxit, ifelse(erreur > tol, {it=it+1}, {break}), {break})
    }
    return(list(sigmasq.v= sigmasq0.v ,
                sigmasq  = sigmasq0 ,
                betas    = beta.gwr0,
                conv     = ifelse(it < maxit, 1, 0),
                n.iter   = it,
                proj     = proj,
                vbeta=vbeta,
                vbeta2 = vbeta2))
}



##########################################################################################################
#########################Startwerte#######################################################################


ConvREBLEUP<-function(beta0,sigmasq0.v,sigmasq0, ni, n,m, tol,maxit,k,X,Y,Z){

const             <- 2 * pnorm(k) - 1 - 2 * k * dnorm(k) + 2 * (k^2) * (1-pnorm(k))
K                 <- diag(rep(const,n))

#STEP 5-6:  Convergence
sigmasq        <- NULL
sigmasq.v      <- NULL
beta           <- NULL
it             <- 1
repeat{

    cat("Iteration sigma&beta",it,"\n",fill=T)

    beta1               <- beta0
    sigmasq1            <- sigmasq0
    sigmasq1.v          <- sigmasq0.v
    U.Inv               <- diag(1/rep(sigmasq0+sigmasq0.v,n))
    U                   <- diag(rep(sigmasq0+sigmasq0.v,n))
    sqrtU               <- sqrt(U)
    G                   <- diag(rep(sigmasq0.v, m))
    V                   <- Z%*%G%*%t(Z) + diag(sigmasq0,n)
    V.Inv               <- vi.inv.NoSp(n=ni, e=sigmasq0, u=sigmasq0.v)
    V.Inv.sq            <- V.Inv %*% V.Inv

    #Update for betas:

    sqrtU           <-diag(sqrt(diag(U)))
    Q               <-t(X)%*%V.Inv%*%sqrtU
    r               <-diag(1/diag(sqrtU))%*%(Y-X%*%beta1)
    psi.fn          <-hub.psi(r, k)
    psi.r           <- psi.fn$psi
    W1              <-diag(c(psi.fn$psi/r))
    A               <-solve(Q%*%W1%*%diag(1/diag(sqrtU))%*%X)%*%Q%*%W1%*%diag(1/diag(sqrtU))
    beta0           <-A%*%Y
#     proj            <-X%*%A
    beta            <-cbind(beta,beta0)
    xbeta           <- X%*%beta0

    #Update for Sigma
    der.psi.r       <- psi.fn$der.psi
    E.psi.der       <- pnorm(k) - pnorm(-k)


    qq              <- V.Inv %*% sqrtU %*% psi.r
    ZZ              <- Z %*% t(Z) # Ableitung von V nach sigma_v


    a1             <- t(qq)%*%qq
    a2             <- t(qq)%*%ZZ%*%qq
    a              <- matrix(c(a1[1],a2[1]), nrow = 2, ncol=1)

    KK2     <-K%*%V.Inv.sq
    KK3     <-K%*%V.Inv%*%ZZ%*%V.Inv

    A1      <-sum(diag(KK2))
    A2      <-sum(diag(KK2%*%ZZ))
    A3      <-sum(diag(KK3))
    A4      <-sum(diag(KK3%*%ZZ))
    AA       <-matrix(c(A1,A3,A2,A4), nrow = 2, ncol=2)

    sigma           <-solve(AA)%*%a
    sigmasq0        <-sigma[1,]
    sigmasq0       <- ifelse(sigmasq0 < 0, -sigmasq0/100, sigmasq0)
    #             sigmasq0       <- ifelse(sigmasq0 > 1000, 1, sigmasq0)
    sigmasq         <- c(sigmasq, sigmasq0)
    sigmasq0.v      <-sigma[2,]
    sigmasq0.v     <- ifelse(sigmasq0.v < 0, -sigmasq0.v/100, sigmasq0.v)
    #             sigmasq0.v     <- ifelse(sigmasq0.v > 1000, sigmasq0.v/100, sigmasq0.v)
    sigmasq.v       <- c(sigmasq.v, sigmasq0.v)#
    erreur          <- sum(abs((beta1-beta0)/beta1))+sum(abs((sigma - c(sigmasq1,sigmasq1.v))/c(sigmasq1,sigmasq1.v)))
#     erreur            <-abs(sigmasq0-sigmasq1)/sigmasq1+  abs(sigmasq0.v-sigmasq1.v)/sigmasq1.v
#     erreur          <- sum((beta1-beta0)^2)+sum((sigma - c(sigmasq1,sigmasq1.v))^2)
    ifelse (it < maxit, ifelse(erreur > tol, {it=it+1}, {break}), {break})
}
return(list(sigmasq.v= sigmasq0.v ,
            sigmasq  = sigmasq0 ,
            beta     = beta0,
            conv     = ifelse(it < maxit, 1, 0),
            n.iter   = it,
            A        = A))


}





