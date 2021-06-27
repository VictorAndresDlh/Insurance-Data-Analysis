# Proyecto Riesgo Actuarial y Financiero

# Librerías necesarias
library(ROI)
library(quantmod)
library(xts)
library(readr)
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
library(MASS)
library(MASSExtra)
library(actuaryr)

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
G8T$VLRASEGU <- ceiling(G8T$VLRASEGU/10^3)


# SINHT
SINHT <- SINH
# Escoger unicamente los valores que tenian como estado final "pagado" 
SINHT <- subset(SINHT, SINHT$ESTADO_FINAL == "Pagado")
# Escoger unicamente los valores asegurados mayores que 0
SINHT <- subset(SINHT, SINHT$VLRASEGU > 0)
# Escoger unicamente los datos correspondientes al 2018
SINHT <- subset(SINHT, year(SINHT$FECHASIN) > 2017)
# Eliminar los primeros 3000 y ultimos 3000
SINHT <- SINHT[order(SINHT$VLRSININCUR),]
SINHT <- head(SINHT, -3000)
SINHT <- tail(SINHT, -3000)
# Escalar los datos
SINHT$VLRSININCUR <- ceiling(SINHT$VLRSININCUR/10^3)


# 500 simulaciones de 1600 muestras aleatorias para promediar la media y la varianza
x <- c()
y <- c()
for (i in 1:500) {
  MA <- SINHT[sample(nrow(SINHT), 1600), ]
  diasMA <- count(MA, MA$FECHASIN)
  a <- descdist(diasMA$n,boot = 50, discrete = T)
  x[i] <- a$mean 
  y[i] <- a$sd
}
mean(x)
mean(y)
mean(y)^(2)/mean(x)


# Como la varianza(n)/media(n) \approx 1, entonces n distribuye Poisson 
# Asumimos que el n de 300 polizas tambien distribuye poisson
# Determinar el parametro: media de las reclamaciones por dia, por lo cual
lambda <- 300/365


# Encontrar la mejor distribucion para el posible monto de reclamacion
descdist(as.numeric(as.vector(unlist(SINHT$VLRSININCUR))), boot = 300,  discrete = F)
  x <- fitdist(as.numeric(as.vector(unlist(SINHT$VLRSININCUR))), "gamma", method = "mme")
plot(x)

xmean <- 5.328861/(0.001809636) 
xvar <- 5.328861/(0.001809636)^2


# Algoritmos

#-------------------------------------------------------------------------------------

# TLC
TLC = function(S_inf,S_sup,E_S,var_S) {
  sapply(S_inf:S_sup,
         function(s){
           pnorm((s-E_S)/sqrt(var_S),0,1)
         })
}
TLC_v = TLC(0, sum(tail(sort(G8T$VLRASEGU), 30)), (lambda*30)*(xmean), (lambda*30)*(xvar)^(2)+((xmean)^(2))*(lambda*30))
#Anexo
library(knitr)
# Tabla
s = 0:(length(TLC_v)-1)
TLC_P = cbind(s,TLC_v)
TLC_PJ = kable(TLC_P,
               caption = "Aproximación con el Teorema Central del Límite", 
               align = c('c', 'c'), 
               col.names = c(" s ", " F(s) -> TLC "), 
               row.names = FALSE,
               digits = 5000)
TLC_PJ

#-------------------------------------------------------------------------------------

# Sparse Vector 
SparseVec <- function(freq) {
  if (any(freq<0)) stop("negative frequency")
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
  if(sum(p)>1||any(p<0)) stop("p parameter not a density")
  if(lambda*sum(p) > 727) stop("Underflow")
  cumul <- f <- exp(-lambda*sum(p))
  r <- length(p)
  s <- 0
  repeat {
    s <- s+1
    m <- min(s, r)
    last <- (lambda/s)*sum(1:m*head(p,m)*rev(tail(f,m)))
    #sort ya hace lo de rev
    f <- c(f,last)
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
p[2:3] <- 0.5; 
lab <- 1
f_fft <- Re(fft(exp(lab*(fft(p)-1)), inverse=TRUE))/n
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

#-------------------------------------------------------------------------------------

#SUMA DE POISSON NO INDEPENDIENTES
SPNI = function(s, j1, j2, j3, n1, n2, n3) {
  SUMT2 = SS1 <- rep(0, s+1)
  SUMT1 = matrix(SS1, nrow = length(SS1), ncol = length(SS1))
  for(j in 0:s) {
    for(i in 0:j) {
      SUMT1[i+1, j+1] = ifelse(i%% n1 == 0 & (j-i) %% n2 == 0, dpois(i/n1, j1)*dpois((j-i)/n2, j2), 0)
    }
    SUMT2[j+1] = ifelse((s-j)%%n3 == 0 , dpois((s-j)/n3, j3),0 )
    S1 = as.numeric(colSums(SUMT1)%*%SUMT2)};
    S1
  #exp(j1+j2+j3)
}

# s=3 , Pr[X=1,2,3]=1/4,2/4,1/4, 1*N_1 + 2*N_2 + 3*N_3
SPNI(3, 1, 2, 1, 1, 2, 3)




TLC = function(S_inf,S_sup,E_S,var_S) {
  sapply(S_inf:S_sup,
         function(s){
           pnorm((s-E_S)/sqrt(var_S),0,1)
         })
}
TLC_v = TLC(0, 20, 10, 3*3+9*3)
s = 0:(length(TLC_v)-1)
TLC_P = cbind(s,TLC_v)
TLC_PJ = kable(TLC_P,
               caption = "Aproximación con el Teorema Central del Límite", 
               align = c('c', 'c'), 
               col.names = c(" s ", " F(s) -> TLC "), 
               row.names = FALSE,
               digits = 5000)
TLC_PJ


























vecprob = c()
for (i in 1:max(SINHT$VLRSININCUR)) {
  vecprob[i] <- dgamma(i, shape = 3.573161e-01, scale = 1/7.801438e-08)
}

panjernb <- function(r, p, vecp) {
  a <- 1-p
  b <- (1-p)*(r-1)
  c <- c()
  if (vecp[1] == 0) {
    c[1] <- 0
  }
  else {
    c[1] <- (p/(1-(1-p)*exp(log(vecp[1]))))^(r)
  }
  d <- c()
  while (sum(c) == 0.99) {
    for (s in 2:243) {
      for (h in 1:i) {
        d[h] <- a+((b*h)/s)*vecp[h+1]*c[s-h]
        c[s] <- (1/(1-a*vecp[1]))*sum(d)
      }
    }
  }
  return(c)
}


sum(G8T$VLRPRISUSCR)

1 - pnbinom(40, size = 26.03031, mu = 40.34247)

# Se seleccionan aleatoriamente 3000 polizas 100 veces para determinar
# la media y la varianza

x <- c()
y <- c()
for (i in 1:100) {
  MA <- SINHT[sample(nrow(SINHT), 3000), ]
  diasMA <- count(MA, MA$FECHASIN)
  a <- descdist(diasMA$n,boot = 50, discrete = T)
  x[i] <- a$mean 
  y[i] <- a$sd
}
mean(x)
mean(y)

(mean(y)^2)/mean(x)

# Como var(n)/E(n) = 1.2 \approx = 1. Entonces n-poisson
# Para determinar lambda, se realizan 500 simulaciones para muestras
# aleatorias de 300 polizas

x <- c()
for (i in 1:500) {
  MA <- SINHT[sample(nrow(SINHT), 300), ]
  diasMA <- count(MA, MA$FECHASIN)
  a <- descdist(diasMA$n,boot = 50, discrete = T)
  x[i] <- a$mean 
}
mean(x)




