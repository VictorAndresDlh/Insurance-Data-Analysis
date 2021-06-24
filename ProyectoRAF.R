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
# Se generan nuevos dataframe G8T y SINHT
G8T <- G8
SINHT <- SINH

# G8T
# Eliminar los registros que tenian prima suscrita = 0
G8T <- subset(G8, G8[3] > 0)
# Cambiar RC = T si efectivamente VLRPRSRC > 0
for (row in 1:nrow(G8T)) {
  if (G8T$VLRPRISUSCR > 0) {
    G8T$RC = T
  }
}
rm(row)

# SINHT
# Escoger unicamente los valores que tenian como estado final "pagado" 
SINHT <- subset(SINHT, SINHT$ESTADO_FINAL == "Pagado")
# Escoger unicamente los valores asegurados mayores que 0
SINHT <- subset(SINHT, SINHT$VLRASEGU > 0)
# Escoger unicamente los datos correspondientes al 2018
SINHT <- subset(SINHT, year(SINHT$FECHASIN) > 2017)

# Encontrar la mejor distribucion para el posible monto de reclamacion
descdist(as.numeric(as.vector(unlist(SINHT$VLRSININCUR))), boot = 300,  discrete = F)
x <- fitdist(as.numeric(as.vector(unlist(SINHT$VLRSININCUR))), "gamma", method = "mme")
plot(x)

# Proceso de conteo
dias <- count(SINHT, SINHT$FECHASIN)
semana <- count(SINHT, week(SINHT$FECHASIN))
meses <- count(SINHT, month(SINHT$FECHASIN))

descdist(dias$n,boot = 50, discrete = T)
plot(fitdist(dias$n, "pois"))
plot(fitdist(dias$n, "nbinom"))
descdist(semana$n,boot = 50, discrete = T)
plot(fitdist(dias$n, "pois"))
plot(fitdist(semana$n, "nbinom", method = "mge"))
