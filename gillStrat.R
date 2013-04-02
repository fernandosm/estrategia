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
# Exclusivo para não-vocalizante (2a) e (2b)
# R3: Sucesso reprodutivo de x, com parâmetro lambda
# R4: Sucesso reprodutivo de y, com parâmetro gama# 
#
#### Propensity functions ####
# a1(x) = k1*y
# a2(x) = k2*x
# a3a(x) = (1-e^(-lambda*y))*x
# a4a(x) = (1-e^(-gama*x))*y
# a3b(x) = (1-e^(-lambda*y*x))
# a4b(x) = (1-e^(-gama*x*y))

# Parameters
propFunctionsA  <-  c('k1*y','k2*(N-y)','(1-exp((-lambda*y)))*(N-y)','(1-exp((-gama*(N-y))))*y') # N=x+y
parms <- c(k1=1,k2=1,lambda=0.03,gama=0.005,N=400);
sCM <- matrix(c(-1,1,1,1),nrow=1,byrow=TRUE);# state change matrix
x0 <- c(y=200)# condição inicial
tFinal <- 10;

############# direct method #############
out.sim <- ssa(x0=x0, a=propFunctionsA, nu=sCM, parms, tf = tFinal, method="D",
               simName = "Sapo model", verbose = F,
               consoleInterval = 10, maxWallTime = 300)

# str(out.sim) # estrutura dos dados
results <- cbind(out.sim$data,x=400-out.sim$data[,2])
colnames(results)[1] <- "timeSeries"
results[1:10,]

plot(results[,1],results[,2],ylim=c(0,410),col="green",pch=16,
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
