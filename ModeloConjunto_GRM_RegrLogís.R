X <- read.csv("~/dados_ENADE_amostra1000.csv", row.names=1)

#Modelo JAGS Estimação Conjunta
Estimação.conjunta <-"
model{
  for(i in 1:N){
    for(j in 1:J){
      X[i,j] ~ dordered.logit(mu[i,j],b[j,1:(c[j]-1)])  # Verossimilhança
      mu[i,j] <- a[j] * theta[i]  # Preditor linear
    }
    y[i] ~ dbin(p[i], 1)
  logit(p[i]) <- beta_0 + beta_1*theta[i]
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
beta_0 ~ dnorm(0.0,1E-3)
beta_1 ~ dnorm(0.0,1E-3)
}
"


simula.conj <- function(dados.TRI, Y) {

  Y <- as.numeric(Y)-1
  
  d <- as.matrix(dados.TRI)
  for (i in 1:ncol(d)) {
    d[,i] <- as.numeric(d[,i])
  }
  c <- c(3,3,3,2,2)
  
  #Ajustando####
  
  dados <- list(N = length(Y),c=c,J=ncol(d),X=as.matrix(d),y = Y)
  
  bettas.inic <- matrix(c(0,1,
                          0,1,
                          0,1,
                          0,NA,
                          0,NA),
                        nrow=ncol(d), ncol=max(c)-1, byrow=T)
  
  inic.model =list(list(a=rep(1,ncol(d)), b=bettas.inic,beta_0=0, beta_1=0,
                        .RNG.seed=10, .RNG.name="base::Mersenne-Twister"),
                   list(a=rep(5,ncol(d)), b=bettas.inic+2,beta_0=5, beta_1=5,
                        .RNG.seed=20, .RNG.name="base::Mersenne-Twister"),
                   list(a=rep(10,ncol(d)), b=bettas.inic+4,beta_0=-5, beta_1=-5,
                        .RNG.seed=30, .RNG.name="base::Mersenne-Twister"))
  
  
  saída.MCMC <- run.jags(model=Estimação.conjunta, 
                         monitor=c("beta_0", "beta_1","a","b","theta"),
                         data=dados, n.chains=3, method="parallel", 
                         inits=inic.model,adapt=1000, burnin=30000,
                         sample=10000,modules="glm")
  return(saída.MCMC)
}

resultado <- simula.conj(X[,1:5],X[,6])
