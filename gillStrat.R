########## Para identificar a pasta da biblioteca (fernandosm olnly) ##########
# setting path
.libPaths("/home/fernandosm/Tools/R/x86_64-pc-linux-gnu-library/2.15")
.libPaths()
########## Para identificar a pasta da biblioteca (fernandosm olnly) ##########

# Biblioteca carrega pacote para plots e stochastic simulation algorithm
library(GillespieSSA)
library(ggplot2)

# Reaction channels:
# R1: 
# R2: 
# R3: 
# Rn:
# 
# Propensity functions
# a1(x) = 
# a2(x) = 
# a3(x) = 
# an(x) = 

# Parameters
a1  <-  
parms <- 
sCM <- # state change matrix
x0 <- # condição inicial
tFinal <- 2000;

############# direct method #############

############# optmized tau-leap method #############

############# fixed tau-leap method #############