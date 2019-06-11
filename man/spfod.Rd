\name{spfod}
\alias{spfod}
\title{spatial filtering origin destination}
\description{Una breve descripcion}
\usage{spfod(mod, W, alfa = 0.05, metodo = TRUE, residuos = TRUE)}
\arguments{
  \item{mod}{Objeto resultado de un lm o glm ajustado a los flujos}
  \item{W}{Matriz W (de los poligonos, NO de los flujos)}
  \item{alfa}{nivel de significacion}
  \item{metodo}{TRUE es chun}
  \item{residuos}{ TRUE se centra con respecto a los regresores}}
\details{Aca debemos explicar alo mas}
\value{El modelo con los vectores propios seleccionados y el historial de valores (y p-valores) del indice de moran}
\author{Eugenia Riano, Martin Hansz, Fernando Massa, Antonio Rey, Andres Castrillejo}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (mod, W, alfa = 0.05, metodo = TRUE, residuos = TRUE) 
{
    stopifnot(any(class(mod) \%in\% c("lm", "glm", "glmerMod", 
        "lmerMod")))
    if (any(c("glmerMod", "lmerMod") \%in\% class(mod))) 
        datos <- mod@frame
    else datos <- mod$model
    if (!"spdep" \%in\% loadedNamespaces()) 
        library(spdep)
    N <- nrow(datos)
    n <- nrow(W)
    X <- model.matrix(mod)
    Wd <- diag(n) \%x\% W
    Wo <- W \%x\% diag(n)
    if (residuos) {
        M <- diag(N) - X \%*\% solve(t(X) \%*\% X) \%*\% t(X)
    }
    else {
        M <- diag(N) - matrix(1, N, 1) \%*\% matrix(1, 1, N)/N
    }
    if (metodo) {
        Wod = 1/2 * (Wd + Wo)
        vpod <- VPcalc(Wod, OD = "od", Mat = M)
        m_od <- Moran(modelo = mod, mtz = Wod, p = 0)
        pvod <- m_od$p.value
        if (pvod > alfa) {
            cat("\n", "No se detecto autocorrelacion espacial", 
                "\n")
            return(list(mod = mod, m = NULL))
        }
        cat("\n")
        print(m_od)
        k <- 0
        minpv <- pvod
        columnas = ncol(vpod)
        while (minpv <= alfa) {
            k <<- k + 1
            vpmod <- NULL
            if (m_od$p.value[nrow(m_od)] <= alfa) {
                REMOD <- forcito(modelo = mod, vp = vpod, datos = datos, 
                  vpmod = vpmod)
                mod <- REMOD$modelo
                vpod <- REMOD$vp
                datos <- REMOD$datos
                vpmod <- REMOD$vpmod
                m_od_remod <- Moran(modelo = mod, mtz = Wod, 
                  p = k)
                m_od <- rbind(m_od, m_od_remod)
                cat("\n", "Origen-Destino", "\n")
                print(m_od)
                pvod <- m_od_remod$p.value
            }
            columnas <- columnas - 1
            sig <- rep(0, columnas)
            minpv <- pvod
        }
        return(list(mod = mod, m = list(OD = m_od), v = list(vpmod)))
    }
    else {
        Wod <- W \%x\% W
        vpo <- VPcalc(Wo, "Eo", Mat = M)
        vpd <- VPcalc(Wd, "Ed", Mat = M)
        vpod <- VPcalc(Wod, "Eod", Mat = M)
        m_o <- Moran(modelo = mod, mtz = Wo, p = 0)
        m_d <- Moran(modelo = mod, mtz = Wd, p = 0)
        m_od <- Moran(modelo = mod, mtz = Wod, p = 0)
        pvo <- m_o$p.value
        pvd <- m_d$p.value
        pvod <- m_od$p.value
        if (min(pvo, pvd, pvod) > alfa) {
            cat("\n", "No se detecto autocorrelacion espacial", 
                "\n")
            return(list(mod = mod, m = NULL))
        }
        cat("\n")
        print(m_o)
        print(m_d)
        print(m_od)
        k <- 0
        minpv <- min(pvo, pvd, pvod)
        columnas <- ncol(vpo)
        sig <- rep(0, columnas)
        vpmod <- NULL
        vpmo <- NULL
        vpmd <- NULL
        while (minpv <= alfa) {
            k <<- k + 1
            if (m_o$p.value[nrow(m_o)] <= alfa) {
                REMOD <- forcito(modelo = mod, vp = vpo, datos = datos, 
                  vpmod = vpmod)
                mod <- REMOD$modelo
                vpo <- REMOD$vp
                datos <- REMOD$datos
                vpmod <- REMOD$vpmod
                mk <- Moran(modelo = mod, mtz = Wo, p = k)
                m_o <- rbind(m_o, mk)
                cat("\n", "Origen", "\n")
                print(m_o)
                pvo <- mk$p.value
            }
            if (m_d$p.value[nrow(m_d)] <= alfa) {
                REMOD <- forcito(modelo = mod, vp = vpd, datos = datos, 
                  vpmod = vpmod)
                mod <- REMOD$modelo
                vpd <- REMOD$vp
                datos <- REMOD$datos
                vpmod <- REMOD$vpmod
                mk <- Moran(modelo = mod, mtz = Wd, p = k)
                m_d <- rbind(m_d, mk)
                cat("\n", "Destino", "\n")
                print(m_d)
                pvd <- mk$p.value
            }
            if (m_od$p.value[nrow(m_od)] <= alfa) {
                REMOD <- forcito(modelo = mod, vp = vpod, datos = datos, 
                  vpmod = vpmod)
                mod <- REMOD$modelo
                vpod <- REMOD$vp
                datos <- REMOD$datos
                vpmod <- REMOD$vpmod
                mk <- Moran(modelo = mod, mtz = Wod, p = k)
                m_od <- rbind(m_od, mk)
                cat("\n", "Origen-Destino", "\n")
                print(m_od)
                pvod <- mk$p.value
            }
            columnas <- columnas - 1
            minpv <- min(pvo, pvd, pvod)
        }
        return(list(mod = mod, m = list(origen = m_o, destino = m_d, 
            OD = m_od), v = list(vpmod)))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
