## File:  internals.R
## Description:  functions - spfod

Moran <- function(modelo,mtz,p)
{
 moranT <- moran.test(residuals(modelo),mat2listw(mtz))
 mT <-data.frame(paso=p,I=as.numeric(moranT$statistic),p.value=moranT$p.value)
 return(mT)
}

VPcalc <- function(mtz,OD,Mat)
{
  LaW <- mtz 
  MWoM  <- Mat %*% (1/2*(LaW+t(LaW))) %*% Mat
  ## vectores propios
    vpo  <- eigen(MWoM)$vectors
  ## saco el asociado al valor propio cero
  vpo<-vpo[,-which.min(apply(vpo,2,var))]
  colnames(vpo)<-paste(OD,1:ncol(vpo),sep='')
  return(vpo)
}

forcito <- function(modelo,vp,datos,vpmod){
    sig <- NULL
    for (i in 1:ncol(vp)){
        f.i<-as.formula(paste('.~.',colnames(vp)[i],sep='+'))
        d.i<-data.frame(datos,vp[,i,drop=FALSE])
        m.i<-update(modelo,f.i,data=d.i)
        sig <- c(sig,1-pchisq(-2*(logLik(modelo)[1]-logLik(m.i)[1]),1))
    }
    vpmax<-which.min(sig)
    datos<-data.frame(datos,vp[,vpmax,drop=FALSE])
    f<-as.formula(paste('.~.',colnames(vp)[vpmax],sep='+'))
    modelo<-update(modelo,f,data=datos)
    
    vpmod <- cbind(vpmod,vp[,vpmax])
    colnames(vpmod)[k] <- colnames(vp)[vpmax]
    vp <- vp[,-vpmax]
    
    return(list(modelo=modelo,vp=vp,datos=datos,vpmod=vpmod))
}

##
##  ---- End Of File ---
##
