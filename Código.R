# Punto 1 -- Estimar modelo de tasa spot
# Modelo Spot: CIR - Vasicek
# Tasa spot TES 10


#Leer datos TES
TES = read.table("TES.dat",header=TRUE,stringsAsFactors=FALSE)
TES <- TES[, c("fecha","ya10")]
r<- TES$ya10

ts.plot(r)

# Fechas 

fechas = rev(as.Date(TES$fecha,"%d/%m/%Y",by="days"))
r=rev(TES$ya10)/100


np = length(r)
n1= which(fechas=="2010-01-04")
y = r[n1:np]
cbind(fechas,y)
np1 = length(y)

ejex.mes = seq(fechas[1],fechas[np1],by="month")
ejex.ano = seq(fechas[1],fechas[np1],by="year")

plot(fechas,y, xaxt="n", panel.first = grid(),type='l',
     ylab='TES10')
axis.Date(1, at=ejex.mes, format="%m/%y")
axis.Date(1, at=ejex.ano, labels = FALSE, tcl = -0.2)


#--------------------------------------------------------Estimación CIR

library(yuima)

cir <- setModel(drift="alpha*(mu-x)", diffusion="sigma*sqrt(x)")

Delta = 1/252
mod.cir <- setYuima(data=setData(r, delta=Delta), model=cir)


source("OU_Calibrate_ML.r")
source("OU_Calibrate_LS.r")

delta = 1/252
(M2 = OU_Calibrate_ML(r,delta))
(M3 = OU_Calibrate_LS(r,delta))

#Parametros

(cbind(c(M2[3],M2[1],M2[2]),c(M3[3],M3[1],M3[2])))



param.init <- list(alpha = M2$lambda, 
                   mu = M2$mu, sigma = M2$sigma)

low.par <- list(alpha=0, mu=0, sigma=0)
upp.par <- list(alpha=500, mu=10,sigma=50)

#------estimacion QMLE

bancol.cir <- qmle(mod.cir, start = param.init,
                   lower = low.par, upper = upp.par,method="L-BFGS-B")

summary(bancol.cir)

pars=bancol.cir@coef


(4*pars[2]*pars[3]/pars[1]^2 > 2)


#-------Simulación

T = np1/252

sampling <- setSampling(Initial = 0.0, 
                        Terminal = T, n = np1, delta=Delta)

bancol.cir.s<-setYuima(model=cir, sampling=sampling)

X<-simulate(bancol.cir.s,xinit=r[1],true.parameter=pars)

Xt=X@data@original.data[-1]
t = seq(0,T,Delta)
plot(fechas,Xt,type='l',ylim=c(0,0.1))
abline(h=0)

for(j in 1:30){
  X<-simulate(bancol.cir.s,xinit=r[1],true.parameter=pars)
  Xt=X@data@original.data[-1]
  lines(fechas,Xt,col='gray')}
lines(fechas,r,col='red',lwd=2)



#-------------------------------PUNTO 2: Modelo MBG para COLCAP
#### Estimación

library(PerformanceAnalytics) # para: Return.calculate()
library(fAssets)  # para: assetsMeanCov()


COL = read.table("colcap.prn",header=TRUE,stringsAsFactors=FALSE)

fechascol <- seq(from = as.Date("2009-01-01"), by="day",length.out=243)

ejex.mes = seq(fechascol[1],fechascol[243],by="month")
ejex.ano = seq(fechascol[1],fechascol[243],by="year")

#-----------estimacion MBG con Yuima

Delta <- 1/243

mbg <- setModel(drift="mu*x", diffusion="sigma*x")

mod <- setYuima(model=mbg, data=setData(COL, delta=Delta))

set.seed(123)
fit <- qmle(mod, start=list(mu=1, sigma=1),
            lower=list(mu=0.1, sigma=0.1),
            upper=list(mu=100, sigma=10))
summary(fit)

y = COL[1:243,]

X <- diff(log(y))
mean(X)
#X <- as.numeric(na.omit(diff(log(COL))))

(alpha <- mean(X)/Delta)
(sigma <- sqrt(var(X)/Delta))
(mu <- alpha +0.5*sigma^2)


#------simulacion
#--------periodo de simulacion y frecuencia

n <- 242
samp <- setSampling(Terminal=T, n=n, delta = Delta)
mod.acciones <- setYuima(model=mbg, sampling=samp)
#-----------simular trayectorias
X = simulate(mod.acciones, 
             true.par=list(mu=mu,sigma=sigma), 
             xinit=COL[1])
Xt=X@data@original.data




plot(fechascol[1:(n+1)],y[1:(n+1)], xaxt="n", panel.first = grid(),
     type='l',col="red",xlab="Fecha", ylab="Valor de la acción")
axis.Date(1, at=ejex.mes, format="%m/%y")
axis.Date(1, at=ejex.ano, labels = FALSE, tcl = -0.2)

for(j in 1:30){
  
  X = simulate(mod.acciones, 
               true.par=list(mu=mu,sigma=sigma), 
               xinit=COL[1])
  Xt=X@data@original.data

lines(fechascol[1:(n+1)],Xt,col='grey')
}  


# PUNTO 3: Modelo de anualidad financiada con portafolio balanceado

require(yuima)

#--------consumo lineal
#ct = function(t)C*(1+rho*floor(t))

#---------consumo geometrico

ct = function(t)C*exp(deltaq*floor(t))
deltaq = mu.p/3

#---------------------valor del contrato por formula 
source("CIrZeroBondPrice.r")

#-------------composicion portafolio balanceado

w = 0.5
#--------------parametros de mbg
mu.p = 0.5115856
sigma.p = 0.1398896

rho = mu.p*0.3
#--------------parametros CIR
sigma = 0.04721173
alpha = 0.99959366
mu =  0.06899461
theta = -(1-w)

lambda = 0
r0 = 0.0903
x = r0
C = 12*10


fn = function(t){ct(t)*exp(w*(-mu.p + w*sigma.p^2)*t)*CIRZeroBondPrice(t, alpha, mu, sigma,lambda, x)}

library(rmutil)


(V = int(fn,0,T))



#-----------definir el modelo

mod2D <- setModel(drift = c("alpha*(mu-x1)",
                            "((1-w)*x1+w*mu.p)*x2-ct(t)"),  
                  diffusion = matrix( c( "sigma*sqrt(x1)","","","w*sigma.p*x2"), 2, 2),
                  solve.variable=c("x1","x2"), state.variable=c("x1","x2"), time.variable="t")
str(mod2D)


T <- 10
n <- T*252
Delta = 1/252
samp <- setSampling(Terminal=T, n=n, delta = Delta)

mod.cir.2D <- setYuima(model=mod2D, sampling=samp)

#-----------simular trayectorias
X = simulate(mod.cir.2D, 
             true.par=list(alpha=alpha,mu=mu,sigma=sigma,mu.p=mu.p,sigma.p=sigma.p), 
             xinit=c(x,1.25*V))
rt=X@data@original.data[,1]
St=X@data@original.data[,2]

par(mfrow=c(1,1))
t = seq(0,T,Delta)
plot(t,rt,type='l')
plot(t,St,type='l',ylim=c(-100,1500))
abline(h=0)
for(j in 1:30){
  X = simulate(mod.cir.2D, 
               true.par=list(alpha=alpha,mu=mu,sigma=sigma,mu.p=mu.p,sigma.p=sigma.p), 
               xinit=c(x,1.25*V))
 
  St=X@data@original.data[,2]
  lines(t,St,col='gray')
}

#-----------densidad tiempo de default
N = 1000
T0 = double(N)

for(j in 1:N){
  X = simulate(mod.cir.2D, 
               true.par=list(alpha=alpha,mu=mu,sigma=sigma,mu.p=mu.p,sigma.p=sigma.p), 
               xinit=c(x,1.25*V))
  St=X@data@original.data[,2]
  T0[j]= ifelse(sum(St < 0)>0,t[min(which(St<0))],NA)
}
par(mfrow=c(2,2))
plot(t,St,type='l')
abline(h=0)
hist(na.omit(T0),40)
plot(density(na.omit(T0)))

#-----------probabilidad de default


( p =1- sum(is.na(T0))/N)


#---------- Punto 4 

#------------funcion para la prob default 

prob.default = function(w){
  
  theta = -(1-w)
  V = int(fn,0,T)

  mod2D <- setModel(drift = c("alpha*(mu-x1)",
                              "((1-w)*x1+w*mu.p)*x2-ct(t)"),  
                    diffusion = matrix( c( "sigma*sqrt(x1)","","","w*sigma.p*x2"), 2, 2),
                    solve.variable=c("x1","x2"), state.variable=c("x1","x2"), time.variable="t")
  str(mod2D)
  
  
  T <- 10
  n <- T*252
  Delta = 1/252
  samp <- setSampling(Terminal=T, n=n, delta = Delta)
  
  mod.cir.2D <- setYuima(model=mod2D, sampling=samp)
  
  N = 1000
  T0 = double(N)
  
  for(j in 1:N){
    X = simulate(mod.cir.2D, 
                 true.par=list(alpha=alpha,mu=mu,sigma=sigma,mu.p=mu.p,sigma.p=sigma.p), 
                 xinit=c(x,1.25*V))
    St=X@data@original.data[,2]
    T0[j]= ifelse(sum(St < 0)>0,t[min(which(St<0))],NA)
  }
  
  ( p =1- sum(is.na(T0))/N)
  
  return(p)
}

c = seq(0.1,1,0.1)
b = double(1)
b[0]=0
j=double(1)
j
for (i in c) {
  
  j =c(j,i)
}

j
b=c(0,0,0,0.013,0.087,0.16,0.202,0.216,0.304, 0.313,0.33)
cbind(j,b)

