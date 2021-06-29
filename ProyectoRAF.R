# Proyecto Riesgo Actuarial y Financiero

# Librerías necesarias
library(ROI)
library(quantmod)
library(xts)
library(readr)
library(MASS)
library(MASSExtra)
library(fitdistrplus)
library(Matrix)
library(modopt.matlab)
library(xtable)
library(matlib)
library(gapminder)
library(flexsurv)
library(tidyverse)
library(ggplot2)
library(arsenal)
library(ordinal)
library(lubridate)
library(actuaryr)
library(tictoc)

# Se cargan las bases de datos

# Base de datos para el grupo 8
G8 <- as.data.frame(read_delim("C:/Users/victo/OneDrive/Escritorio/Proyecto Riesgo Actuarial y Financiero/Grupo_P8.txt", 
                 "|", escape_double = FALSE, col_types = cols(FECINICIO = col_datetime(format = "%Y-%m-%d"), 
                                                              FECFIN = col_datetime(format = "%Y-%m-%d"), 
                                                              PTD = col_logical(), PPD = col_logical(), 
                                                              PH = col_logical(), PPH = col_logical(), 
                                                              RC = col_logical()), trim_ws = TRUE))
# Base de datos de los Siniestros Históricos
SINH <- as.data.frame(read_delim("C:/Users/victo/OneDrive/Escritorio/Proyecto Riesgo Actuarial y Financiero/Siniestros_Hist.txt", 
                   "|", escape_double = FALSE, col_types = cols(FECHASIN = col_date(format = "%Y-%m-%d"), 
                                                                FECPAGOAMP = col_date(format = "%Y%m%d")), 
                              trim_ws = TRUE))

# Resumen de las dos bases de datos
summary(G8)
summary(SINH)


# Depuracion de los datos


# G8T
G8T <- G8
# Eliminar los registros que tenian prima suscrita = 0
G8T <- subset(G8, G8$VLRPRISUSCR > 0)
# Cambiar RC = T si efectivamente VLRPRSRC > 0
for (i in 1:nrow(G8T)) {
  if (G8T$VLRPRISUSCR[i] > 0) {
    G8T$RC = T
  }
}
# Escalar los datos
G8T$VLRASEGU <- G8T$VLRASEGU/10^6


# Se crea un data frame temporal en el caso de ser necesario volver al original
SINHT <- SINH
# Escoger unicamente los valores que tenian como estado final "pagado" 
SINHT <- subset(SINHT, SINHT$ESTADO_FINAL == "Pagado")
# Escoger unicamente los valores asegurados mayores que 0
SINHT <- subset(SINHT, SINHT$VLRASEGU > 0)
# Escoger unicamente los datos correspondientes al 2018
SINHT <- subset(SINHT, year(SINHT$FECHASIN) > 2017)
# Eliminar datos cuyo valor deducible sea mayor que valor siniestro
SINHT <- subset(SINHT, SINHT$VLRSININCUR >= SINHT$VLRDEDUCIBLE)
# Eliminar los datos cuyo valor pagado es 0
SINHT <- subset(SINHT, SINHT$VLRPAGADO. > 0)
# Eliminar los datos cuyo valor de siniestro es 0
SINHT <- subset(SINHT, SINHT$VLRSININCUR > 0)

# Escalar los datos
SINHT$VLRSININCUR <- SINHT$VLRSININCUR/10^6


#-------------------------------------------------------------------------------------
# Solucion de problemas

# 1. Cual es el valor esperado de S?
# Es necesario modelar X y N:

# Modelando X

# Encontrar la mejor distribucion para el posible monto de reclamacion
descdist(as.numeric(as.vector(unlist(SINHT$VLRSININCUR))), boot = 300,  discrete = F)
x <- fitdist(as.numeric(as.vector(unlist(SINHT$VLRSININCUR))), "gamma", method = "mme")
plot(x)

xmean <- 0.4273345/0.07793957
xvar <- 0.4273345/(0.07793957)^2


# Modelando N

# 500 simulaciones de 2000 muestras aleatorias para promediar la media y la varianza
x <- c()
y <- c()
for (i in 1:500) {
  MA <- SINHT[sample(nrow(SINHT), 2000), ]
  diasMA <- count(MA, MA$FECHASIN)
  a <- descdist(diasMA$n,boot = 50, discrete = T)
  x[i] <- a$mean 
  y[i] <- a$sd
}
mean(x)
mean(y)
mean(y)^(2)/mean(x)
rm(x)
rm(y)


# Como la varianza(n)/media(n) \approx 1, entonces n distribuye Poisson 
# Asumimos que el n de 300 polizas tambien distribuye poisson

# Determinar el parametro: media de las reclamaciones por dia, por lo cual
lambda <- sum(((G8T$VLRPRISUSCR/10^6)/1.1)/G8T$VLRASEGU)


#-------------------------------------------------------------------------------------
# 2. Probabilidad de que la reserva no alcance a cubrir el total de reclamaciones
# Necesitamos modelar la S del año.

# Reserva

resv <- 0.9*sum((G8T$VLRPRISUSCR/10^6)/(1.1))

# Algoritmos para el año

# TLC APLICADO
TLC = function(S_inf, S_sup, E_S, var_S) {
  sapply(S_inf:S_sup,
         function(s) {
           pnorm((s-E_S)/sqrt(var_S), 0, 1)
         })
}
TLC_v_an = TLC(0, sum(G8T$VLRASEGU), (lambda)*(xmean), (lambda)*mean((as.numeric(as.vector(unlist(SINHT$VLRSININCUR))))^2))

# Tabla
s = 0:(length(TLC_v_an)-1)


# Se generan los valores de las probabilidades
TLC_P_an = cbind(s,TLC_v_an)
for (i in 1:10) {
  print((1-TLC_P_an[floor((i/10)*(resv))+2,][2])*100)
}


# Panjer para Poisson APLICADO
Panjer.Poisson <- function(p, lambda) {
  if(sum(p) > 1 || any(p < 0)) stop("p parameter not a density")
  if(lambda*sum(p) > 727) stop("Underflow")
  cumul <- f <- exp(-lambda*sum(p))
  r <- length(p)
  s <- 0
  repeat {
    s <- s + 1
    m <- min(s, r)
    last <- (lambda/s)*sum(1:m*head(p, m)*rev(tail(f, m)))
    #sort ya hace lo de rev
    f <- c(f, last)
    cumul <- cumul + last
    if (cumul > 0.99999999)
      break
  }
  return(f)
}
P = integer(160)
P[1] = 0
for (i in 2:160) {
  P[i-1] <- dgamma(i-1, shape = 0.4273345, rate = 0.07793957)
}
PPan = Panjer.Poisson(P, lambda)

# Tabla
s = 0:(length(PPan)-1)
PPan = cbind(s, PPan)

# Se genera la columna para la distribucion
PPan <- cbind(PPan, integer(nrow(PPan)))
PPan[1, 3] <- PPan[1, 2]
PPan
for (i in 2:nrow(PPan)) {
  PPan[i, 3] <- PPan[i-1, 3] + PPan[i, 2]
}

# Se generan los valores de las probabilidades
for (i in 1:10) {
  print((1-PPan[floor((i/10)*(resv))+2,][3]))
}

#------------------------------------------------------------------------------------
# 3. Valor minimo de reserva para responder con el total de las reclamaciones con un c %
# de probabilidad

# Interpolacion lineal con TLC
TLC_P_an[85, ]
TLC_P_an[86, ]

# Interpolacion lineal con poisson
PPan[91, ]
PPan[92, ]

#------------------------------------------------------------------------------
# 4. Qué tanto tiempo se puede extender la cobertura de las polizas para que 
# la reserva siga cubriendo el monto total de reclamaciones con una
# probabilidad del 95% ? 






























#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------

# Sparse Vector 
SparseVec <- function(freq) {
  if (any(freq < 0)) stop("negative frequency")
  M <- length(freq)
  mu <- sum((1:M)*freq);
  sigma2 <- sum((1:M)^2*freq)
  
  ##mean and variance of the compound r.v.; see
  (3.4)
  MM <- ceiling(mu + 10* sqrt(sigma2)) + 6
  fs <- dpois(0:(MM-1), freq[1]) 
  
  ##density of S_1 = 1*N_1
  for (j in 2 :M) {
    MMM <- trunc((MM-1)/j)
    fj <- rep(0, MM)
    
    ##construct the density of j*N_j
    fj[(0:MMM)*j+1] <- dpois(0:MMM, freq[j])
    fs <- convolve(fs, rev(fj), type="o")
  }
  
  ##fs is the density of S_j = 1*N_1 + ... + j*N_j, j=2..M
  return(fs) 
}
#***********************************************

f <- SparseVec(c(1,2,1));
#***********************************************

#Anexo
# Tabla
s = 0:(length(f)-1)
f = cbind(s,f)
colnames(f) = c("s","f(s)")
#Mostrar 10 valores
head(f,10)

#-------------------------------------------------------------------------------------

# Panjer para Poisson
Panjer.Poisson <- function(p, lambda) {
  if(sum(p) > 1||any(p < 0)) stop("p parameter not a density")
  if(lambda*sum(p) > 727) stop("Underflow")
  cumul <- f <- exp(-lambda*sum(p))
  r <- length(p)
  s <- 0
  repeat {
    s <- s + 1
    m <- min(s, r)
    last <- (lambda/s)*sum(1:m*head(p,m)*rev(tail(f,m)))
    #sort ya hace lo de rev
    f <- c(f, last)
    cumul <- cumul + last
    if (cumul > 0.99999999)
      break
  }
  return(f)
}
PP = Panjer.Poisson(c(1/4,1/2,1/4),4)
#Anexo
library(knitr)
# Tabla
s = 0:(length(PP)-1)
PP = cbind(s,PP)
kable( PP , caption = "Recursividad de Panjer",
       align = c('c', 'c'),
       col.names = c(" s "," f(s) "),
       row.names = FALSE,
       digits = 2000)

#-------------------------------------------------------------------------------------

# Transformada Rapida de Fourier
n <- 40;
p <- rep(0,n);
p
p[2:3] <- 0.5; 
p
lab <- 1
f_fft <- Re(fft(exp(lab*(fft(p)-1)), inverse = TRUE))/n
#Anexo
library(knitr)
# Tabla
s = 0:(length(f_fft)-1)
f_fft = cbind(s,f_fft)
kable( f_fft , caption = "Recursividad de Panjer con bluce For",
       align = c('c', 'c'),
       col.names = c(" s "," f(s) "),
       row.names = FALSE,
       digits = 5000)












