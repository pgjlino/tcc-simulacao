Sys.setenv(JAGS_HOME= "C:\\Program Files\\JAGS\\JAGS-4.3.1")
#library(rjags)
#library(R2jags)
library(runjags)
library(ggplot2)
library(patchwork)
library(coda)
library(mirt) #gerar dados simulados TRI

#===============================================
# Definir os diretórios
#===============================================
# Limpar o ambiente global
rm(list=ls())


n.amostra <- c(50,300,1000)
n.itens <- c(5,15,30)
# Caminho das pastas
main.dir<-setwd("C:/Users/pauli/OneDrive/Documentos/TCC 1")

# Determinando os diretórios
Parm.dir<-paste(main.dir, '/Parametros reais',sep='')
Dados.dir<-paste(main.dir,'/DadosGerados',sep='')
Dados.Tri.dir <- matrix(NA,3,3)
colnames(Dados.Tri.dir) <- c("i=5","i=15","i=30")
rownames(Dados.Tri.dir) <- c("n=50","n=300","i=1000")
for (n in 1:length(n.amostra)) {
  for (i in 1:length(n.itens)) {
    Dados.Tri.dir[n,i] <- paste(Dados.dir,'/TRI/n',n.amostra[n],
                                'i',n.itens[i],sep='')
    
  }
  
}
Results.dir<-paste(main.dir, '/Resultados',sep='')

getwd()



#===============================================
# Gerando parâmetros reais
#===============================================
#---------------------------
#Thetas reais (habilidades)
#---------------------------
for (n in n.amostra) {
  theta<-matrix(NA,n,1)
  colnames(theta)<-c("theta")
  set.seed(214856+n)
  theta[,1] <- rnorm(n)
  write.csv(theta,file=(paste(Parm.dir,"/theta",n,".csv",sep=""))
            , row.names=F)
}

#---------------------------
#Parâmetros dos itens reais
#---------------------------

for (j in n.itens) {
  param.itens <- matrix(NA,j,3)
  colnames(param.itens)<-c("a",paste('b',seq(2,1),sep='')) 
  row.names(param.itens)<-paste("I",seq(1,j),sep='')
  set.seed(45162+j)
  param.itens <- cbind(a = runif(j,.75,1.33),
                       b2 = runif(j, 0.75, 1.75), 
                       b1 = runif(j, -1.75, -0.75))
  write.csv(param.itens,file=(paste(Parm.dir,"/param.itens",j,".csv",sep="")))
}
#===============================================
# Carregando os parâmetros reais e as sementes
#===============================================
thetas.reais <- list()
for (r in 1:length(n.amostra)) {
  n<-n.amostra[r]
  theta<-read.csv(file=paste(Parm.dir,"/theta",n,".csv",sep=""))
  thetas.reais[[length(thetas.reais) + 1]] <-theta
  names(thetas.reais)[length(thetas.reais)] <- paste("theta",n,sep="")
}

param.itens.reais <- list()
for (r in 1:length(n.itens)) {
  j <-n.itens[r]
  itens<-read.csv(file=paste(Parm.dir,"/param.itens",j,".csv",sep=""))
  param.itens.reais[[length(param.itens.reais) + 1]] <-itens[,2:4]
  names(param.itens.reais)[length(param.itens.reais)] <- paste("param.itens",j,sep="")
}

#===============================================
# Gerando os dados da TRI e estimando TRI apenas
#===============================================

Gerar_dados_TRI <- "
data {
for(i in 1:N){
    for(j in 1:J){
      X[i,j] ~ dordered.logit(mur[i,j],br[j,1:(c[j]-1)])  # Verossimilhança
      mur[i,j] <- ar[j] * thetar[i]  # Preditor linear
    }
  }
}
model{
  for(i in 1:N){
    for(j in 1:J){
      X[i,j] ~ dordered.logit(mu[i,j],b[j,1:(c[j]-1)])  # Verossimilhança
      mu[i,j] <- a[j] * theta[i]  # Preditor linear
    }
    theta[i] ~ dnorm(0,1)  # priori variável latente
  }
    for(j in 1:J){
    a[j] ~ dlnorm(0,16) #log normal prior
}
  for(j in 1:J){
    for(k in 1:(c[j]-1)){
      b[j,k] ~ dnorm(0,0.01)  # priori dificuldades
    }
  }
}
"

simula.TRI <- function(nitem, namostra) {
  n <- match(namostra,n.amostra) 
  i <- match(nitem,n.itens) 
  parâmetros_estimados <- matrix(NA,500,(namostra+nitem*3))
  sd.par_estimados <- matrix(NA,500,(namostra+nitem*3))
  for (m in 1:500) {
    c <- rep(3,nitem)
    #Ajustando####
    bettas.inic <- matrix(rep(c(-1.5,1.5),nitem),
                          nrow=nitem, ncol=max(c)-1, byrow=T)
    
    inic.model =list(list(a=rep(1,nitem), b=bettas.inic))
    
    br=cbind(param.itens.reais[[i]]$b1,param.itens.reais[[i]]$b2)
    
    data.grm = list(N=namostra,c=c,J=nitem,ar=param.itens.reais[[i]]$a,
                    br=br,thetar=thetas.reais[[n]]$theta)
 
    saída.MCMC <- run.jags(Gerar_dados_TRI, monitor=c("X","theta","a","b"),
                           data=data.grm,n.chains=1, method="rjags",
                           inits=inic.model, adapt=500, burnin=500,
                           sample=1000,modules="glm")
    resul <- summary(saída.MCMC)
    
    X <- matrix(resul[1:(namostra*nitem),2],namostra,nitem,byrow = F)
    write.csv(X,file=(paste(Dados.Tri.dir[n,i],"/replica",m,".csv",sep="")))
    
    
    #Extraindo estimativas EAP dos parametros
    parâmetros_estimados[m,] <- resul[(namostra*nitem+1):(namostra*nitem+namostra+3*nitem),4]
    sd.par_estimados[m,] <- resul[(namostra*nitem+1):(namostra*nitem+namostra+3*nitem),5]
    
    print(m)
    write.csv(parâmetros_estimados,
              file=(paste(Results.dir,"/ApenasTRI/estim.par.n",namostra,
                          "i",nitem,".csv",sep="")))
    write.csv(sd.par_estimados,
              file=(paste(Results.dir,"/ApenasTRI/sd.estim.par.n",namostra,
                          "i",nitem,".csv",sep="")))
  }
}


for (n in c(50,300,1000)) {
  for (i in c(5,15,30)) {
    tempo <- system.time({
      simula.TRI(i,n)
    })[3]
    
    write.csv(tempo,
              file=(paste(Results.dir,"/ApenasTRI/tempo",n,
                          "i",i,".csv",sep="")))
    
  }
  
}


#---------------------------
#Modelo de Regressão logistica 
#Gerar e estimação ingênua:
#---------------------------
Gerar_dados_reg <- "
data {
 for (i in 1:N) {
    y[i] ~ dbin((ilogit(0.5 + theta_real[i])), 1)
}
}
model
{
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dbin((ilogit(beta_0 + beta_1*theta[i])), 1)
    
  }

  # Priors
beta_0 ~ dnorm(0,0.001)
beta_1 ~ dnorm(0,0.001)
}
"

#logit(p) = b_0 + b_1*theta
#yi ~ binom(1,p)


simula.regressao <- function(nitem, namostra) {
  n = c(50,300,1000)
  j = c(15,30)
  betas.estimados <- matrix(NA,500,2)
  sd_betas_estimados <- matrix(NA,500,2)
  for (m in 1:500) {
    resultados <-  read.csv(file=paste(Results.dir,"/ApenasTRI/estim.par.n",namostra,
                                       "i",nitem,".csv",sep=""),row.names=1, quote="'")
    x <- as.vector(resultados[m,1:namostra])
    
    #Estima os parâmetros
    
    inits1 <- list(beta_0=0, beta_1=0)
    
    dados <- list(N = namostra,theta = x,theta_real=thetas.reais[[n]]$theta)
    
    saída.MCMC <- run.jags(model=Gerar_dados_reg, 
                           monitor=c("y","beta_0", "beta_1"),
                           data=dados, n.chains=1, method="rjags", 
                           inits=list(inits1),adapt=500, burnin=500,
                           sample=1000,modules="glm")
    resul <- summary(saída.MCMC)
    y <- matrix(resul[1:namostra,2],namostra,1,byrow = T)
    
    write.csv(y,file=(paste(Dados.dir,"/Regressão logistica/Y.n=",namostra,".",
                            nitem,"/replica",m,".csv",sep="")))
    betas.estimados[m,] <- resul[(namostra+1):(namostra+2),4]
    sd_betas_estimados[m,] <- resul[(namostra+1):(namostra+2),5]
    print(m)
    write.csv(sd_betas_estimados,
              file=(paste(Results.dir,"/Regressão logística/theta estimado/sd.betas n"
                          ,namostra,"i",nitem,".csv",sep="")))
    write.csv(betas.estimados,
              file=(paste(Results.dir,"/Regressão logística/theta estimado/betas n"
                          ,namostra,"i",nitem,".csv",sep="")))
    
  }
  
}

for (n in c(50,300,1000)) {
  for (i in c(5,15,30)) {
    
    tempo <- system.time({
      simula.regressao(i,n)
    })[3]
    write.csv(tempo,
              file=(paste(Results.dir,"/Regressão logística/theta estimado/tempo"
                          ,n,"e",i,".csv",sep="")))
    
  } 
}


#---------------------------
#Estimação conjunta
#---------------------------

TRI.Regres <-"
model{
for(i in 1:N){
  for(j in 1:J){
    X[i,j]~dordered.logit(mu[i,j],b[j,1:(c[j]-1)])
    mu[i,j] <- a[j]*theta[i]

  }
  theta[i]~dnorm(0,1)
}
  for(j in 1:J){
    a[j] ~ dlnorm(0,16) #log normal prior
}
for(j in 1:J){
  for(k in 1:(c[j]-1)){
    b[j,k]~dnorm(0,0.01) #prior dos limites
  }
}
for (i in 1:N) {
    y[i] ~ dbin((ilogit(beta_0 + beta_1*theta[i])), 1)
    
  }
beta_0 ~ dnorm(0,0.001)
beta_1 ~ dnorm(0,0.001)

}
"



Conjunta <- function(namostra, nitem) {
  n.amostra <- c(50,300,1000)
  n.itens <- c(5,15,30)
  n <- match(namostra,n.amostra) 
  i <- match(nitem,n.itens)
  par_estimados <- matrix(NA,500,(2+namostra+nitem*3))
  sd.par_estimados <- matrix(NA,500,(2+namostra+nitem*3))

  for (m in 1:500) {
    #Ajustando####
    c <- rep(3,nitem)
    x <- read.csv(file=paste(Dados.dir,"/TRI/n",
                             namostra,"i",nitem,"/replica",m,".csv",sep=""), row.names=1)
    Y <- read.csv(file=paste(Dados.dir,"/Regressão logistica/Y.n=",
                             namostra,".",nitem,"/replica",m,".csv",sep=""), row.names=1)
    
    data = list(N=namostra,c=c,J=nitem,X=as.matrix(x),y=as.matrix(Y))
    
    bettas.inic <- matrix(rep(c(-1.5,1.5),nitem),
                          nrow=nitem, ncol=max(c)-1, byrow=T)
    
    inic.model =list(list(a=rep(1,nitem), b=bettas.inic,beta_0=0, beta_1=0))
      saída.MCMC <- run.jags(TRI.Regres, monitor=c("beta_0","beta_1","theta","a","b"),
                             data=data,n.chains=1, method="rjags",
                             inits=inic.model, adapt=500, burnin=500,
                             sample=1000,modules = "glm")
      resul <- summary(saída.MCMC)
      
    #Extraindo estimativas EAP dos parametros
    par_estimados[m,] <- resul[,4]
    sd.par_estimados[m,] <- resul[,5]
    print(m)
    write.csv(par_estimados,
              file=(paste(Results.dir,"/Tri e Regressão/param n"
                          ,namostra,"i",nitem,".csv",sep="")))
    write.csv(sd.par_estimados,
              file=(paste(Results.dir,"/Tri e Regressão/sd.param n"
                          ,namostra,"i",nitem,".csv",sep="")))
  }
}

for (n in c(50,300,1000)) {
  for (i in c(5,15,30)) {
    
    tempo <- system.time({
      Conjunta(n,i)
    })[3]
    write.csv(tempo,
              file=(paste(Results.dir,"/Tri e Regressão/betas n"
                          ,n,"i",i,"tempo.csv",sep="")))
    
  } 
}

#---------------------------
#Resultados
#---------------------------

calcular_eqm <- function(est_pontual, sd, valor_real) {
  # Cálculo do viés
  vies <- est_pontual - valor_real
  
  # Cálculo do EQM
  eqm <- vies^2 + sd^2
  
  # Retorna o EQM
  return(eqm)
}
vies<-function (x,y) {
  vies = sum((x - y),na.rm = TRUE)/(length(x))
  return (vies)
}


plot.resultados.total <- function(beta){
  V1 <- numeric()
  V2 <- numeric()
  V3 <- numeric()
  for (n in c(50,300,1000)) {
    for (i in c(5,15,30)) {
      V1 <-c(V1,read.csv(paste(main.dir,"/Resultados/Regressão logística/theta estimado/betas n",n,
                               "i",i,".csv",sep=""),row.names=1, quote="'")[,beta])
      V2 <- c(V2,rep(i,500))
      V3 <- c(V3,rep(n,500))
      
    }
  }
  
  dados1 <- cbind(V1,V2,V3)
  dados1 <- as.data.frame(dados1)
  dados1[,2] <- as.factor(dados1[,2])
  dados1[,3] <- as.factor(dados1[,3])
  dados1[,4] <- as.factor(rep("Estimação Ingênua",nrow(dados1)))
  
  V1 <- numeric()
  V2 <- numeric()
  V3 <- numeric()
  for (n in c(50,300,1000)) {
    for (i in c(5,15,30)) {
      V1 <-c(V1,read.csv(paste(main.dir,"/Resultados/Tri e Regressão/param n",n,
                               "i",i,".csv",sep=""),row.names=1, quote="'")[,beta])
      V2 <- c(V2,rep(i,500))
      V3 <- c(V3,rep(n,500))
      
    }
  }
  
  dados2 <- cbind(V1,V2,V3)
  dados2 <- as.data.frame(dados2)
  dados2[,2] <- as.factor(dados2[,2])
  dados2[,3] <- as.factor(dados2[,3])
  dados2[,4] <- as.factor(rep("Estimação Conjunta",nrow(dados2)))
  
  dados <- rbind(dados1,dados2)
  if(beta==1){
    nome= expression(hat(beta)[0])
    plot <- ggplot(dados,aes(x = V3, y = V1,fill = V2)) +
      geom_errorbar(stat = "boxplot", width = .2, 
                    position = position_dodge(width=0.5))+
      geom_boxplot(width=0.5,outlier.shape = 1,outlier.size = 2) +
      scale_fill_manual(values = c("#7DBDDE","#677ECB","#324486"),
                        labels = c("5 Itens","15 Itens","30 Itens"),
                        name = "Quantidade de Itens")+
      labs(y = nome,x ="Tamanho Amostral")+ 
      scale_x_discrete(labels=c("n=50","n=300","n=1000"))+
      geom_hline(yintercept = 0.5,color="red",linetype='dotted',size=.8)+
      theme_bw()
    
  }
  else{
    nome=expression(hat(beta)[1])
    
    plot <- ggplot(dados,aes(x = V3, y = V1,fill = V2)) +
      geom_errorbar(stat = "boxplot", width = .2, 
                    position = position_dodge(width=0.5))+
      geom_boxplot(width=0.5,outlier.shape = 1,outlier.size = 2) +
      scale_fill_manual(values = c("#7DBDDE","#677ECB","#324486"),
                        labels = c("5 Itens","15 Itens","30 Itens"),
                        name = "Quantidade de Itens")+
      labs(y = nome,x ="Tamanho Amostral")+ 
      scale_x_discrete(labels=c("n=50","n=300","n=1000"))+
      geom_hline(yintercept = 1,color="red",linetype='dotted',size=.8)+
      theme_bw()
  }
  
  return(plot+facet_grid(. ~ V4))
}
plot.resultados.total(2)

plot.resultados.total.eqm <- function(beta){
  V1 <- numeric()
  V2 <- numeric()
  V3 <- numeric()
  b <- ifelse(beta==1,0.5,1)
  for (n in c(50,300,1000)) {
    for (i in c(5,15,30)) {
      sd <-read.csv(paste(main.dir,"/Resultados/Regressão logística/theta estimado/sd.betas n",n,
                          "i",i,".csv",sep=""),row.names=1, quote="'")[,beta]
      est <-read.csv(paste(main.dir,"/Resultados/Regressão logística/theta estimado/betas n",n,
                           "i",i,".csv",sep=""),row.names=1, quote="'")[,beta]
      V1 <- c(V1,sqrt(calcular_eqm(est,sd,b)))
      V2 <- c(V2,rep(i,500))
      V3 <- c(V3,rep(n,500))
      
    }
  }
  
  dados1 <- cbind(V1,V2,V3)
  dados1 <- as.data.frame(dados1)
  dados1[,2] <- as.factor(dados1[,2])
  dados1[,3] <- as.factor(dados1[,3])
  dados1[,4] <- as.factor(rep("Estimação Ingênua",nrow(dados1)))
  
  V1 <- numeric()
  V2 <- numeric()
  V3 <- numeric()
  for (n in c(50,300,1000)) {
    for (i in c(5,15,30)) {
      sd <-read.csv(paste(main.dir,"/Resultados/Tri e Regressão/sd.param n",n,
                          "i",i,".csv",sep=""),row.names=1, quote="'")[,beta]
      est <-read.csv(paste(main.dir,"/Resultados/Tri e Regressão/param n",n,
                           "i",i,".csv",sep=""),row.names=1, quote="'")[,beta]
      V1 <- c(V1,sqrt(calcular_eqm(est,sd,b)))
      V2 <- c(V2,rep(i,500))
      V3 <- c(V3,rep(n,500))
    }
  }
  
  dados2 <- cbind(V1,V2,V3)
  dados2 <- as.data.frame(dados2)
  dados2[,2] <- as.factor(dados2[,2])
  dados2[,3] <- as.factor(dados2[,3])
  dados2[,4] <- as.factor(rep("Estimação Conjunta",nrow(dados2)))
  
  dados <- rbind(dados1,dados2)
  if(beta==1){
    nome= expression(RMSE ~ hat(beta)[0])
    plot <- ggplot(dados,aes(x = V3, y = V1,fill = V2)) +
      geom_errorbar(stat = "boxplot", width = .2, 
                    position = position_dodge(width=0.5))+
      geom_boxplot(width=0.5,outlier.shape = 1,outlier.size = 2) +
      scale_fill_manual(values = c("#7DBDDE","#677ECB","#324486"),
                        labels = c("5 Itens","15 Itens","30 Itens"),
                        name = "Quantidade de Itens")+
      labs(y = nome,x ="Tamanho Amostral")+ 
      scale_x_discrete(labels=c("n=50","n=300","n=1000"))+
      theme_bw()
    
  }
  else{
    nome=expression(RMSE ~ hat(beta)[1])
    
    plot <- ggplot(dados,aes(x = V3, y = V1,fill = V2)) +
      geom_errorbar(stat = "boxplot", width = .2, 
                    position = position_dodge(width=0.5))+
      geom_boxplot(width=0.5,outlier.shape = 1,outlier.size = 2) +
      scale_fill_manual(values = c("#7DBDDE","#677ECB","#324486"),
                        labels = c("5 Itens","15 Itens","30 Itens"),
                        name = "Quantidade de Itens")+
      labs(y = nome,x ="Tamanho Amostral")+ 
      scale_x_discrete(labels=c("n=50","n=300","n=1000"))+
      theme_bw()
  }
  
  return(plot+facet_grid(. ~ V4))
}
plot.resultados.total.eqm(2)

#TRI##
media.vies.estimativas.TRI <- function(nitem, namostra){
  n <- match(namostra,n.amostra) 
  i <- match(nitem,n.itens) 
  t_real <- thetas.reais[[n]]$theta
  param_itens_reais <- as.matrix(param.itens.reais[[i]])
  resultados <-  read.csv(file=paste(main.dir,"/Resultados/ApenasTRI/estim.par.n",namostra,
                                     "i",nitem,".csv",sep=""),row.names=1, quote="'")
  theta.resul <- resultados[,1:namostra]
  a.resul <- resultados[,(namostra+1):(namostra+nitem)]
  b.resul <- resultados[,(namostra+1+nitem):ncol(resultados)]
  vies.theta <- rep(NA,500)
  for (m in 1:500) {
    vies.rep <- rep(NA,namostra)
    for (r in 1:namostra) {
      vies.rep[r] <- vies(theta.resul[m,r],t_real[r])
    }
    vies.theta[m] <- mean(abs(vies.rep)) 
  }
  
  vies.a <- rep(NA,500)
  for (m in 1:500) {
    vies.rep <- rep(NA,nitem)
    for (r in 1:nitem) {
      vies.rep[r] <- vies(a.resul[m,r],param_itens_reais[r,1])
    }
    vies.a[m] <- mean(abs(vies.rep)) 
  }
  
  vies.b1 <- rep(NA,500)
  for (m in 1:500) {
    vies.rep <- rep(NA,nitem)
    for (r in 1:nitem) {
      vies.rep[r] <- vies(b.resul[m,r],param_itens_reais[r,3])
    }
    vies.b1[m] <- mean(abs(vies.rep)) 
  }
  
  vies.b2 <- rep(NA,500)
  for (m in 1:500) {
    vies.rep <- rep(NA,nitem)
    for (r in 1:nitem) {
      vies.rep[r] <- vies(b.resul[m,(r+nitem)],param_itens_reais[r,2])
    }
    vies.b2[m] <- mean(abs(vies.rep)) 
  }
  
  return(list(vies.theta,vies.a,vies.b1,vies.b2))
}
media.sd.estimativas.TRI <- function(nitem, namostra){
  n <- match(namostra,n.amostra) 
  i <- match(nitem,n.itens) 
  resultados <-  read.csv(file=paste(main.dir,"/Resultados/ApenasTRI/sd.estim.par.n",namostra,
                                     "i",nitem,".csv",sep=""),row.names=1, quote="'")
  theta.resul <- resultados[,1:namostra]
  a.resul <- resultados[,(namostra+1):(namostra+nitem)]
  b.resul <- resultados[,(namostra+1+nitem):ncol(resultados)]
  sd.theta <- rep(NA,500)
  for (m in 1:500) {
    sd.theta[m] <- mean(as.matrix(theta.resul[m,]))
  }
  
  sd.a <- rep(NA,500)
  for (m in 1:500) {
    sd.a[m] <- mean(as.matrix(a.resul[m,])) 
  }
  
  sd.b1 <- rep(NA,500)
  for (m in 1:500) {
    sd.b1[m] <- mean(as.matrix(b.resul[m,1:nitem])) 
  }
  
  sd.b2 <- rep(NA,500)
  for (m in 1:500) {
    sd.b2[m] <- mean(as.matrix(b.resul[m,(nitem+1):(2*nitem)]))
  }
  
  return(list(sd.theta,sd.a,sd.b1,sd.b2))
}
#media.eqm.estimativas.TRI <- function(nitem, namostra){
n <- match(namostra,n.amostra) 
i <- match(nitem,n.itens) 
t_real <- thetas.reais[[n]]$theta
param_itens_reais <- as.matrix(param.itens.reais[[i]])
repli <- read.csv(paste(main.dir,"/Replicas/n",namostra,"e",nitem,".csv",sep=""),row.names=1, quote="'")[,1]
sd <-  read.csv(file=paste(Results.dir,"/ApenasTRI/sd.estim.par.n",namostra,
                           "i",nitem,".csv",sep=""),row.names=1, quote="'")[repli,]
pont <-  read.csv(file=paste(Results.dir,"/ApenasTRI/estim.par.n",namostra,
                             "i",nitem,".csv",sep=""),row.names=1, quote="'")[repli,]
theta.resul <- pont[,1:namostra]
sd.theta.resul <- sd[,1:namostra]
a.resul <- pont[,(namostra+1):(namostra+nitem)]
sd.a.resul <- sd[,(namostra+1):(namostra+nitem)]
b.resul <- pont[,(namostra+1+nitem):ncol(pont)]
sd.b.resul <- sd[,(namostra+1+nitem):ncol(sd)]
eqm.theta <- rep(NA,length(repli))
for (m in 1:length(repli)) {
  eqm.theta[m] <- mean(as.matrix(sqrt(calcular_eqm(theta.resul[m,],sd.theta.resul[m,],t_real))))
}
eqm.a <- rep(NA,length(repli))
for (m in 1:length(repli)) {
  eqm.a[m] <- mean(as.matrix(sqrt(calcular_eqm(a.resul[m,],sd.a.resul[m,],param_itens_reais[,1]))))
}

eqm.b1 <- rep(NA,length(repli))
for (m in 1:length(repli)) {
  eqm.b1[m] <- mean(as.matrix(sqrt(calcular_eqm(b.resul[m,(1:nitem)],sd.b.resul[m,(1:nitem)],param_itens_reais[,3]))))
}

eqm.b2 <- rep(NA,length(repli))
for (m in 1:length(repli)) {
  eqm.b2[m] <- mean(as.matrix(sqrt(calcular_eqm(b.resul[m,((nitem+1):(2*nitem))],sd.b.resul[m,((nitem+1):(2*nitem))],param_itens_reais[,2]))))
}

return(list(eqm.theta,eqm.a,eqm.b1,eqm.b2))
}

media.vies.estimativas.TRI.conj <- function(nitem, namostra){
  n <- match(namostra,n.amostra) 
  i <- match(nitem,n.itens) 
  t_real <- thetas.reais[[n]]$theta
  param_itens_reais <- as.matrix(param.itens.reais[[i]])
  resultados <-  read.csv(file=paste(main.dir,"/Resultados/Tri e Regressão/param n",namostra,
                                     "i",nitem,".csv",sep=""),row.names=1, quote="'")[,-c(1,2)]
  theta.resul <- resultados[,1:namostra]
  a.resul <- resultados[,(namostra+1):(namostra+nitem)]
  b.resul <- resultados[,(namostra+1+nitem):ncol(resultados)]
  vies.theta <- rep(NA,500)
  for (m in 1:500) {
    vies.rep <- rep(NA,namostra)
    for (r in 1:namostra) {
      vies.rep[r] <- vies(theta.resul[m,r],t_real[r])
    }
    vies.theta[m] <- mean(abs(vies.rep)) 
  }
  
  vies.a <- rep(NA,500)
  for (m in 1:500) {
    vies.rep <- rep(NA,nitem)
    for (r in 1:nitem) {
      vies.rep[r] <- vies(a.resul[m,r],param_itens_reais[r,1])
    }
    vies.a[m] <- mean(abs(vies.rep)) 
  }
  
  vies.b1 <- rep(NA,500)
  for (m in 1:500) {
    vies.rep <- rep(NA,nitem)
    for (r in 1:nitem) {
      vies.rep[r] <- vies(b.resul[m,r],param_itens_reais[r,3])
    }
    vies.b1[m] <- mean(abs(vies.rep)) 
  }
  
  vies.b2 <- rep(NA,500)
  for (m in 1:500) {
    vies.rep <- rep(NA,nitem)
    for (r in 1:nitem) {
      vies.rep[r] <- vies(b.resul[m,(r+nitem)],param_itens_reais[r,2])
    }
    vies.b2[m] <- mean(abs(vies.rep)) 
  }
  
  return(list(vies.theta,vies.a,vies.b1,vies.b2))
}
media.sd.estimativas.TRI.conj <- function(nitem, namostra){
  n <- match(namostra,n.amostra) 
  i <- match(nitem,n.itens) 
  resultados <-  read.csv(file=paste(main.dir,"/Resultados/Tri e Regressão/sd.param n",namostra,
                                     "i",nitem,".csv",sep=""),row.names=1, quote="'")[,-c(1,2)]
  theta.resul <- resultados[,1:namostra]
  a.resul <- resultados[,(namostra+1):(namostra+nitem)]
  b.resul <- resultados[,(namostra+1+nitem):ncol(resultados)]
  sd.theta <- rep(NA,500)
  for (m in 1:500) {
    sd.theta[m] <- mean(as.matrix(theta.resul[m,]))
  }
  
  sd.a <- rep(NA,500)
  for (m in 1:500) {
    sd.a[m] <- mean(as.matrix(a.resul[m,])) 
  }
  
  sd.b1 <- rep(NA,500)
  for (m in 1:500) {
    sd.b1[m] <- mean(as.matrix(b.resul[m,1:nitem])) 
  }
  
  sd.b2 <- rep(NA,500)
  for (m in 1:500) {
    sd.b2[m] <- mean(as.matrix(b.resul[m,(nitem+1):(2*nitem)]))
  }
  
  return(list(sd.theta,sd.a,sd.b1,sd.b2))
}
#media.eqm.estimativas.TRI.conj <- function(nitem, namostra){
n <- match(namostra,n.amostra) 
i <- match(nitem,n.itens) 
t_real <- thetas.reais[[n]]$theta
param_itens_reais <- as.matrix(param.itens.reais[[i]])
repli <- read.csv(paste(main.dir,"/Replicas/n",namostra,"e",nitem,".csv",sep=""),row.names=1, quote="'")[,1]
sd <-  read.csv(file=paste(Results.dir,"/Tri e Regressão/sd.param n",namostra,
                           "i",nitem,".csv",sep=""),row.names=1, quote="'")[repli,-c(1,2)]
pont <-  read.csv(file=paste(Results.dir,"/Tri e Regressão/param n",namostra,
                             "i",nitem,".csv",sep=""),row.names=1, quote="'")[repli,-c(1,2)]
theta.resul <- pont[,1:namostra]
sd.theta.resul <- sd[,1:namostra]
a.resul <- pont[,(namostra+1):(namostra+nitem)]
sd.a.resul <- sd[,(namostra+1):(namostra+nitem)]
b.resul <- pont[,(namostra+1+nitem):ncol(pont)]
sd.b.resul <- sd[,(namostra+1+nitem):ncol(sd)]
eqm.theta <- rep(NA,length(repli))
for (m in 1:length(repli)) {
  eqm.theta[m] <- mean(as.matrix(sqrt(calcular_eqm(theta.resul[m,],sd.theta.resul[m,],t_real))))
}
eqm.a <- rep(NA,length(repli))
for (m in 1:length(repli)) {
  eqm.a[m] <- mean(as.matrix(sqrt(calcular_eqm(a.resul[m,],sd.a.resul[m,],param_itens_reais[,1]))))
}

eqm.b1 <- rep(NA,length(repli))
for (m in 1:length(repli)) {
  eqm.b1[m] <- mean(as.matrix(sqrt(calcular_eqm(b.resul[m,(1:nitem)],sd.b.resul[m,(1:nitem)],param_itens_reais[,3]))))
}

eqm.b2 <- rep(NA,length(repli))
for (m in 1:length(repli)) {
  eqm.b2[m] <- mean(as.matrix(sqrt(calcular_eqm(b.resul[m,((nitem+1):(2*nitem))],sd.b.resul[m,((nitem+1):(2*nitem))],param_itens_reais[,2]))))
}

return(list(eqm.theta,eqm.a,eqm.b1,eqm.b2))
}

boxplot.vies <- function(i,n,p){
  V1 <- numeric()
  V2 <- numeric()
  V3 <- numeric()
  for (k in 1:length(n)) {
    for (m in 1:length(i)) {
      v <- media.vies.estimativas.TRI(i[m],n[k])[[p]]
      V1 <- c(V1,v)
      V2 <- c(V2,rep(m,length(v)))
      V3 <- c(V3,rep(k,length(v)))
    }
    
  }
  
  dados1 <- cbind(V1,V2,V3)
  dados1 <- as.data.frame(dados1)
  dados1[,2] <- as.factor(dados1[,2])
  dados1[,3] <- as.factor(dados1[,3])
  dados1[,4] <- as.factor(rep("Estimação Ingênua",nrow(dados1)))
  
  V1 <- numeric()
  V2 <- numeric()
  V3 <- numeric()
  for (k in 1:length(n)) {
    for (m in 1:length(i)) {
      v <- media.vies.estimativas.TRI.conj(i[m],n[k])[[p]]
      V1 <- c(V1,v)
      V2 <- c(V2,rep(m,length(v)))
      V3 <- c(V3,rep(k,length(v)))
    }
    
  }
  
  dados2 <- cbind(V1,V2,V3)
  dados2 <- as.data.frame(dados2)
  dados2[,2] <- as.factor(dados2[,2])
  dados2[,3] <- as.factor(dados2[,3])
  dados2[,4] <- as.factor(rep("Estimação Conjunta",nrow(dados2)))
  dados <- rbind(dados1,dados2)
  
  boxplot <- ggplot(dados,aes(x = V3, y = V1,fill = V2)) +
    geom_errorbar(stat = "boxplot", width = .2, 
                  position = position_dodge(width=0.5))+
    geom_boxplot(width=0.5,outlier.shape = 1,outlier.size = 2) +
    scale_fill_manual(values = c("#7DBDDE","#677ECB","#324486"),
                      labels = c("5 Itens","15 Itens","30 Itens"),
                      name = "Quantidade de Itens")+
    labs(y = expression('Vieses dos Estimadores'),x ="Tamanho Amostral")+ 
    #geom_text(aes(label=ifelse(abs(V1)>0.75,V1,"")), hjust = -1)+
    scale_x_discrete(labels=c("n=50","n=300","n=1000"))+
    theme_bw()
  
  return(boxplot+facet_grid(. ~ V4))
}
boxplot.sd <- function(i,n,p){
  V1 <- numeric()
  V2 <- numeric()
  V3 <- numeric()
  for (k in 1:length(n)) {
    for (m in 1:length(i)) {
      v <- media.sd.estimativas.TRI(i[m],n[k])[[p]]
      V1 <- c(V1,v)
      V2 <- c(V2,rep(m,length(v)))
      V3 <- c(V3,rep(k,length(v)))
    }
    
  }
  
  dados1 <- cbind(V1,V2,V3)
  dados1 <- as.data.frame(dados1)
  dados1[,2] <- as.factor(dados1[,2])
  dados1[,3] <- as.factor(dados1[,3])
  dados1[,4] <- as.factor(rep("Estimação Ingênua",nrow(dados1)))
  
  V1 <- numeric()
  V2 <- numeric()
  V3 <- numeric()
  for (k in 1:length(n)) {
    for (m in 1:length(i)) {
      v <- media.sd.estimativas.TRI.conj(i[m],n[k])[[p]]
      V1 <- c(V1,v)
      V2 <- c(V2,rep(m,length(v)))
      V3 <- c(V3,rep(k,length(v)))
    }
    
  }
  
  dados2 <- cbind(V1,V2,V3)
  dados2 <- as.data.frame(dados2)
  dados2[,2] <- as.factor(dados2[,2])
  dados2[,3] <- as.factor(dados2[,3])
  dados2[,4] <- as.factor(rep("Estimação Conjunta",nrow(dados2)))
  dados <- rbind(dados1,dados2)
  
  boxplot <- ggplot(dados,aes(x = V3, y = V1,fill = V2)) +
    geom_errorbar(stat = "boxplot", width = .2, 
                  position = position_dodge(width=0.5))+
    geom_boxplot(width=0.5,outlier.shape = 1,outlier.size = 2) +
    scale_fill_manual(values = c("#7DBDDE","#677ECB","#324486"),
                      labels = c("5 Itens","15 Itens","30 Itens"),
                      name = "Quantidade de Itens")+
    labs(y = expression('SD dos Estimadores'),x ="Tamanho Amostral")+ 
    scale_x_discrete(labels=c("n=50","n=300","n=1000"))+
    theme_bw()
  
  return(boxplot+facet_grid(. ~ V4))
}

boxplot.vies(c(5,15,30),c(50,300,1000),2)+labs(y = expression(bar(B)[abs](hat(bold(a)))),x ="Tamanho Amostral")
boxplot.sd(c(5,15,30),c(50,300,1000),2)+labs(y = expression(bar(S)(hat(bold(a)))),x ="Tamanho Amostral")

boxplot.vies(c(5,15,30),c(50,300,1000),3)+labs(y = expression(bar(B)[abs](hat(bold(b))[. %.% 1])),x ="Tamanho Amostral")
boxplot.sd(c(5,15,30),c(50,300,1000),3)+labs(y = expression(bar(S)(hat(bold(b))[. %.% 1])),x ="Tamanho Amostral")

boxplot.vies(c(5,15,30),c(50,300,1000),4)+labs(y = expression(bar(B)[abs](hat(bold(b))[. %.% 2])),x ="Tamanho Amostral")
boxplot.sd(c(5,15,30),c(50,300,1000),4)+labs(y = expression(bar(S)(hat(bold(b))[. %.% 2])),x ="Tamanho Amostral")

boxplot.vies(c(5,15,30),c(50,300,1000),1)+labs(y = expression(bar(B)[abs](hat(bold(theta)))),x ="Tamanho Amostral")
boxplot.sd(c(5,15,30),c(50,300,1000),1)+labs(y = expression(bar(S)(hat(bold(theta)))),x ="Tamanho Amostral")

