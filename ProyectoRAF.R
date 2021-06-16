# Proyecto Riesgo Actuarial y Financiero

# Librerías necesarias
library(ROI)
library(quantmod)
library(xts)
library(readr)
library(fitdistr)
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

# Se cargan las bases de datos

# Base de datos para el grupo 8
G8 <- read_delim("C:/Users/victo/OneDrive/Escritorio/Proyecto Riesgo Actuarial y Financiero/Grupo_P8.txt", 
                 "|", escape_double = FALSE, col_types = cols(FECINICIO = col_datetime(format = "%Y-%m-%d"), 
                                                              FECFIN = col_datetime(format = "%Y-%m-%d"), 
                                                              PTD = col_logical(), PPD = col_logical(), 
                                                              PH = col_logical(), PPH = col_logical(), 
                                                              RC = col_logical()), trim_ws = TRUE)
# Base de datos de los Siniestros Históricos
SINH <- read_delim("C:/Users/victo/OneDrive/Escritorio/Proyecto Riesgo Actuarial y Financiero/Siniestros_Hist.txt", 
                   "|", escape_double = FALSE, col_types = cols(FECHASIN = col_date(format = "%Y-%m-%d"), 
                                                                FECPAGOAMP = col_date(format = "%Y%m%d")), 
                              trim_ws = TRUE)

str(G8)
summary(G8)

# ANALISIS DE G8
# Histograma y distribucion de Valor de la prima suscrita incluida RC
G8 %>%
  select(VLRPRISUSCR) %>%
  ggplot(aes(x = VLRPRISUSCR)) +
  geom_histogram(colour = "white")

descdist(as.numeric(unlist(G8[,3])), boot = 300, discrete  = FALSE)
plot(fitdist(as.numeric(unlist(G8[,3])), "norm"))

# Histograma y distribucion de Valor de la prima suscrita incluida RC
G8 %>%
  select(VLRASEGU) %>%
  ggplot(aes(x = VLRASEGU)) +
  geom_histogram(colour = "white")

descdist(as.numeric(unlist(G8[,4])), boot = 300, discrete  = FALSE)
plot(fitdist(as.numeric(unlist(G8[,4]))/10^9, "beta"))

# Histograma y distribucion de Valor de la prima suscrita incluida RC
G8 %>%
  select(VLRASEGURC) %>%
  ggplot(aes(x = VLRASEGURC)) + 
  geom_histogram(colour = "white")

descdist(as.numeric(unlist(G8[,5])), boot = 300, discrete  = FALSE)
plot(fitdist(as.numeric(unlist(G8[,5])), "unif"))

#ANALISIS DE SINH
# Histograma y distribucion de Valor de la prima suscrita incluida RC
SINH %>% 
  select(VLRPRIMAPAG) %>% 
  ggplot(aes(x = VLRPRIMAPAG)) +
  geom_histogram(colour = "white")

descdist(as.numeric(unlist(SINH[,2])), boot = 300, discrete  = FALSE)
plot(fitdist(as.numeric(unlist(SINH[,2])), "gamma", method = "mle", lower=c(0,0), start=list(scale=1,shape=1)))

# Histograma y distribucion de Valor de la prima suscrita incluida RC
SINH %>% 
  select(VLRPRISUSCR) %>% 
  ggplot(aes(x = VLRPRISUSCR)) +
  geom_histogram(colour = "white")

descdist(as.numeric(unlist(SINH[,3])), discrete  = FALSE)

# Histograma y distribucion de Valor de la prima suscrita incluida RC
SINH %>% 
  select(VLRASEGU) %>% 
  ggplot(aes(x = VLRASEGU)) +
  geom_histogram(colour = "white")

descdist(as.numeric(unlist(SINH[,4])), discrete  = FALSE)

# Histograma y distribucion de Valor de la prima suscrita incluida RC
SINH %>% 
  select(VLRPAGADO.) %>% 
  ggplot(aes(x = VLRPAGADO.)) +
  geom_histogram(colour = "white")

descdist(as.numeric(unlist(SINH[,6])), discrete  = FALSE)

# Histograma y distribucion de Valor de la prima suscrita incluida RC
SINH %>% 
  select(VLRSININCUR) %>% 
  ggplot(aes(x = VLRSININCUR)) +
  geom_histogram(colour = "white")

descdist(as.numeric(unlist(SINH[,7])), discrete  = FALSE)

# Histograma y distribucion de Valor de la prima suscrita incluida RC
SINH %>% 
  select(VLRDEDUCIBLE) %>% 
  ggplot(aes(x = VLRDEDUCIBLE)) +
  geom_histogram(colour = "white")

descdist(as.numeric(unlist(SINH[,8])), discrete  = FALSE)

# Histograma y distribucion de Valor de la prima suscrita incluida RC
SINH %>% 
  select(VLRRECOBRO) %>% 
  ggplot(aes(x = VLRRECOBRO)) +
  geom_histogram(colour = "white")

descdist(as.numeric(unlist(SINH[,9])), discrete  = FALSE)

# Histograma y distribucion de Valor de la prima suscrita incluida RC
SINH %>% 
  select(VLRRSVACONSTIAMP) %>% 
  ggplot(aes(x = VLRRSVACONSTIAMP)) +
  geom_histogram(colour = "white")

descdist(as.numeric(unlist(SINH[,10])), discrete  = FALSE)

# Histograma y distribucion de Valor de la prima suscrita incluida RC
SINH %>% 
  select(VLRRSVAPAGAMP) %>% 
  ggplot(aes(x = VLRRSVAPAGAMP)) +
  geom_histogram(colour = "white")

descdist(as.numeric(unlist(SINH[,11])), boot = 300, discrete  = F)
