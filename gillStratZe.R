########## Para identificar a pasta da biblioteca (fernandosm olnly) ##########
# setting path
.libPaths("/home/fernandosm/Tools/R/x86_64-pc-linux-gnu-library/2.15")
.libPaths()
########## Para identificar a pasta da biblioteca (fernandosm olnly) ##########

# Biblioteca carrega pacote para plots e stochastic simulation algorithm
library(GillespieSSA)
library(ggplot2)

#### Reaction channels ####
# Para Vocalizante e não-vocalizantes (1)
# R1: Não-vocalizante que se torna vocalizante
# R2: Vocalizante que se torna não-vocalizante
# Exclusivo para não-vocalizante (Zé)
# R3: Insucesso reprodutivo de y com a presença de x
# R4: Sucesso reprodutivo de y, com par k5
# R5: Sucesso reprodutivo de x com com par k3
# Não-vocalizante e parasita
# R6: Insucessso reprodutivo de y com presença de z
# R7: Insucessso reprodutivo de x com presença de z
# R8: Sucessso reprodutivo de z
# R9: Sucessso reprodutivo de z com presença de x
# R10: Insucessso reprodutivo de z com presença de y
#
#### Propensity functions ####
# a1(y) = k1*X
# a2(x) = k2*X
# a3(x,y) = kc*X*Y
# a4(y) = k5*Y
# a5(x) = k3*X
# a6(x) = kd*Y*Z
# a7(x) = ke*X*Z
# a8(x) = k6*Z
# a9(x) = kf*X*Z
# a10(x) = kg*Y*Z

# Parameters
propFunctionsA  <-  c('k1*Y','k2*(N-Y-Z)','kc*(N-Y-Z)*Y','k5*Y','k3*(N-Y-Z)',
                      'kd*Y*Z','ke*(N-Y-Z)*Z','k6*Z','kf*(N-Y-Z)*Z','kg*Y*Z') # N=x+y+z
parms <- c(k1=0.5,k2=0.5,k5=0.03,k3=0.06,kc=0.1,N=400,
           kd=0.015,ke=0.000001,k6=0.01,kf=0.03,kg=0.03);
sCM <- matrix(c(-1,1,-1,1,1,-1,-1,0,0,0,
                0,0,0,0,0,0,0,1,1,-1),nrow=2,byrow=TRUE);# state change matrix
x0 <- c(Y=100,Z=10)# condição inicial
tFinal <- 1;

############# direct method #############
out.sim <- ssa(x0=x0, a=propFunctionsA, nu=sCM, parms, tf = tFinal, method="D",
               simName = "Sapo model", verbose = F,
               consoleInterval = 10, maxWallTime = 300)
out.sim
# str(out.sim) # estrutura dos dados
results <- cbind(out.sim$data,X=(400-out.sim$data[,2]-out.sim$data[,3]))
colnames(results)[1] <- "timeSeries"
results[1:10,]

plot(results[,1],results$,ylim=c(0,410),col="green",pch=16,
     xlab='tempo',ylab='Número de Indivíduos',cex.lab=1.4,cex.main=1.6,main='Sapo - Estocástico\nMétodo Direto')
points(results[,1],results[,3],col="black",pch=16)
legend(x='topright',c("Não-Voc","Voc"), cex=1, 
       col=c("green","black"),pch=c(16,16),bty='n')


############# optmized tau-leap method (fast) #############
out.sim2 <- ssa(x0=x0, a=propFunctionsA, nu=sCM, parms, tf = tFinal, method="OTL",
                simName = "Sapo model", verbose = F,
                epsilon = 0.03, # only applicable for OTL
                nc = 10,   # only applicable for OTL
                hor = NaN,  # only applicable for OTL
                dtf = 10,   # only applicable for OTL
                nd = 10,  # only applicable for OTL
                consoleInterval = 10, maxWallTime = 300);

# str(out.sim) # estrutura dos dados
results <- cbind(out.sim2$data,x=400-out.sim2$data[,2])
colnames(results)[1] <- "timeSeries"
results[1:10,]

plot(results[,1],results[,2],ylim=c(0,410),col="green",pch=16,
     xlab='tempo',ylab='Número de Indivíduos',cex.lab=1.4,cex.main=1.6,main='Sapo - Estocástico\nTau-leap eficiente')
points(results[,1],results[,3],col="black",pch=16)
legend(x='topright',c("Não-Voc","Voc"), cex=1, 
       col=c("green","black"),pch=c(16,16),bty='n')

############# fixed tau-leap method (faster) #############
out.sim3 <- ssa(x0=x0, a=propFunctionsA, nu=sCM, parms, tf = tFinal, method="ETL",
                simName = "SIR model 2", verbose = F,tau=0.05,
                consoleInterval = 10, maxWallTime = 300);

# str(out.sim) # estrutura dos dados
results <- cbind(out.sim3$data,x=400-out.sim3$data[,2])
colnames(results)[1] <- "timeSeries"
results[1:10,]

plot(results[,1],results[,2],ylim=c(0,410),col="green",pch=16,
     xlab='tempo',ylab='Número de Indivíduos',cex.lab=1.4,cex.main=1.6,main='Sapo - Estocástico\nTau-leap fixo')
points(results[,1],results[,3],col="black",pch=16)
legend(x='topright',c("Não-Voc","Voc"), cex=1, 
       col=c("green","black"),pch=c(16,16),bty='n')
