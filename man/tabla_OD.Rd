\name{tabla_OD}
\alias{tabla_OD}
\docType{data}
\title{Tabla de Origen Destino - Municipios}
\description{data.frame con las frecuencias de origen-destino. En este
  caso Municipios}
\usage{data("tabla_OD")}
\format{
  A data frame with 64 observations on the following 3 variables.
  \describe{
    \item{\code{Orig}}{a factor with levels \code{A} \code{B} \code{C} \code{CH} \code{D} \code{E} \code{F} \code{G}}
    \item{\code{Dest}}{a factor with levels \code{A} \code{B} \code{C} \code{CH} \code{D} \code{E} \code{F} \code{G}}
    \item{\code{counts}}{a numeric vector}
  }
}
\details{Podemos poner algo de la construccion}
\examples{
data(tabla_OD)
## maybe str(tabla_OD) ; plot(tabla_OD) ...
}
\keyword{datasets}
%> tabla_OD
%   Orig Dest counts
%1     A    A     60
%2     B    A      1
%3     C    A      2
%4    CH    A      1
%5     D    A      3
%6     E    A      0
%7     F    A     11
%8     G    A     10
