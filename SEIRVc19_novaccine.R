library(RColorBrewer)
library(deSolve)
library("tidyverse")
library("readxl")
library('socialmixr')

c19age <- function(time, state, pars) {
  with(
     as.list(c(state, pars)), {
  	n=pars$nages; gamma=pars$gamma; gam_a = pars$gam_a; gam_h = pars$gam_h
   h = pars$h; f = pars$f; sigma = pars$sigma; eps = pars$eps; mrate  = pars$mrate;
   gam_p = pars$gam_p
   rzero = pars$rzerovar
   
   beta=rzero[time]*pars$beta;
	
  	S=as.numeric(state[1:n])
  	i=1; E=as.numeric(state[(i*n+1):((i+1)*n)])
  	i=2; A=as.numeric(state[(i*n+1):((i+1)*n)])
  	i=3; I=as.numeric(state[(i*n+1):((i+1)*n)])
  	i=4; R=as.numeric(state[(i*n+1):((i+1)*n)])
  	i=5; H=as.numeric(state[(i*n+1):((i+1)*n)])
  	i=6; P=as.numeric(state[(i*n+1):((i+1)*n)])
  	i=7; D=as.numeric(state[(i*n+1):((i+1)*n)])
  	N=S+E+A+I+R+H+P

  	Imat=matrix(as.numeric(I)/as.numeric(N),nrow=n,ncol=n,byrow=T)
  	Emat=matrix(as.numeric(E)/as.numeric(N),nrow=n,ncol=n,byrow=T)
  	Amat=matrix(as.numeric(A)/as.numeric(N),nrow=n,ncol=n,byrow=T)
  	Smat=matrix(as.numeric(S),nrow=n,ncol=n)
  	
  	newinfections = as.numeric(rowSums(Imat*Smat*beta)+eps*rowSums(Amat*Smat*beta))  #this has dimension nages	
    dS <- -newinfections 
    dE <- newinfections -sigma*E                            #latent
    dA <- f*sigma*E - gam_a*A                               #asymptomatic
    dI <- (1-f)*sigma*E - gamma * I                 #symptomatic
    dH <- h*gamma*I - gam_h*H                               #hospitalized
    dR <- (1-h)*gamma * I + (1-mrate)*gam_h*H + gam_a*A   #recovered
    dP <- mrate*gam_h*H - gam_p*P                                       #died
    dD <- gam_p*P
    
    return(list(c(dS, dE, dA, dI, dR, dH, dP, dD)))
     })}

agestats=data.frame(agegroup=c('0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39','40-44','45-49',
                               '50-54','55-59','60-64','65-69','70-74','75-79','80+'),
                    pyramid=c(0.058,0.063,0.058,0.055,0.063,0.068,0.068,0.066,0.061,0.068,0.070,0.064,0.054,0.050,0.049,0.033,0.049),
                    IFR=c(0.0,0.0,0.0,0.0001,0.0002,0.0004,0.0006,0.001,0.001,0.002,0.005,0.008,0.017,0.027,0.043,0.062,0.096),
                    Hosp=c(0.001,0.001,0.001,0.002,0.005,0.01,0.016,0.023,0.029,0.039,0.058,0.072,0.102,0.117,0.146,0.177,0.18),
                    HFR=c(0.038,0.038,0.038,0.038,0.038,0.038,0.038,0.04,0.045,0.056,0.078,0.113,0.169,0.232,0.291,0.348,0.535)
)

nages=17;
SIMTIME=365*1.5

iperiod=2.7; 
gamma=rep(1/iperiod,nages); 
sigma=rep(1/(5.2),nages); 
gam_a=rep(1/(iperiod),nages); 
gam_h=rep(1/3,nages)
gam_p=rep(1/8,nages); 
eps=rep(1,nages)

h=agestats$Hosp #risk of hospitalisation
f=seq(0.75,0.25,-0.5/(nages-1)) #fraction asymptomatic
mrate=agestats$HFR #Hospital fatality rate
deathrate=agestats$IFR #infection fatality rate

C1=contact_matrix(polymod, countries = "United Kingdom", age.limits = c(0,seq(5,75,5), 79),quiet=T)
betamatrix=C1$matrix
mev=max(Re(eigen(iperiod*betamatrix)$value))

rzero=3.89
betamatrix = betamatrix/mev
beta1=betamatrix
rzerovar = rep(rzero,SIMTIME)

Nuk=55e6
S=Nuk*agestats$pyramid
E=rep(0,nages)
A=rep(0,nages)
I=rep(0,nages)
H=rep(0,nages)
R=rep(0,nages)
P=rep(0,nages)
D=rep(0,nages)
N=S+E+A+I+H+R+P

#start with a small amount of infection in the third age group:

S[3]=S[3]-1
I[3]=I[3]+1

init <- c(S=S, E=E, A=A, I=I, R=R, H=H, P=P, D=D)
params <- list(nages=nages,beta = betamatrix, gamma = gamma,gam_a=gam_a,gam_h=gam_h,gam_p=gam_p,f=f,h=h,mrate=mrate,sigma=sigma,eps=eps,rzerovar=rzerovar)
times <- seq(1, SIMTIME, by = 1)
out <- as.data.frame(ode(y = init, times = times, func = c19age, parms = params))
time=out[,1]
S=out[,2:(nages+1)]
i=1; E=out[(i*nages+2):((i+1)*nages+1)]
i=2; A=out[(i*nages+2):((i+1)*nages+1)]
i=3; I=out[(i*nages+2):((i+1)*nages+1)]
i=4; R=out[(i*nages+2):((i+1)*nages+1)]
i=5; H=out[(i*nages+2):((i+1)*nages+1)]
i=6; P=out[(i*nages+2):((i+1)*nages+1)]
i=7; D=out[(i*nages+2):((i+1)*nages+1)]
N=S+E+A+I+H+R+P

Infected=rowSums(I)
peaktime=which.max(Infected)
totdeaths=diff(rowSums(D))
inhosp=rowSums(H)+rowSums(P)
peaktime=which.max(rowSums(I))

march23=which.min(abs(totdeaths[1:peaktime]-54))#54 is the number of deaths on 23rd March
march30=which.min(abs(totdeaths[1:peaktime]-200))#54 is the number of deaths on 30th March
print(paste("lockdown starts on",march23))


breakpoints=seq(march30,SIMTIME,28)
breakpoints=seq(march30,SIMTIME,28)
rtstages=rep(0.6*rzerovar[1],length(breakpoints))
rtstages[1:6]=c(0.8,0.9,1.0,1.2,1.2,1.2)
rtstages[8]=0.5*rzerovar[1]
rtstages[9]=1.3

rtstages[10:14]=1.6
rtstages[15:17]=0.7*rzerovar[1]
   
rovertime=rzerovar
for(x in 1:length(rtstages))
{
   startdate=breakpoints[x]
   init2 <- cbind(S[startdate,], E[startdate,], A[startdate,], I[startdate,], R[startdate,], H[startdate,], P[startdate,],D[startdate,])

   rt=rtstages[x]
   rzerovar[(startdate+1):SIMTIME] = rt
   betamatrix=C1$matrix
   mev=max(Re(eigen(iperiod*betamatrix)$value))
   betamatrix2 = betamatrix/mev
   params <- list(nages=nages,beta = betamatrix2, gamma = gamma,gam_a=gam_a,gam_h=gam_h,gam_p=gam_p,f=f,h=h,mrate=mrate,sigma=sigma,eps=eps,rzerovar=rzerovar)
   times2 <- seq(startdate, SIMTIME, by = 1)
   out2 <- as.data.frame(ode(y = as.numeric(init2), times = times2, func = c19age, parms = params))

   time2=out2[,1]
   S2=out2[,2:(nages+1)]
   i=1; E2=out2[(i*nages+2):((i+1)*nages+1)]
   i=2; A2=out2[(i*nages+2):((i+1)*nages+1)]
   i=3; I2=out2[(i*nages+2):((i+1)*nages+1)]
   i=4; R2=out2[(i*nages+2):((i+1)*nages+1)]
   i=5; H2=out2[(i*nages+2):((i+1)*nages+1)]
   i=6; P2=out2[(i*nages+2):((i+1)*nages+1)]
   i=7; D2=out2[(i*nages+2):((i+1)*nages+1)]
      
   S[startdate:SIMTIME,] = S2
   E[startdate:SIMTIME,] = E2
   A[startdate:SIMTIME,] = A2
   I[startdate:SIMTIME,] = I2
   R[startdate:SIMTIME,] = R2
   H[startdate:SIMTIME,] = H2
   P[startdate:SIMTIME,] = P2
   D[startdate:SIMTIME,] = D2
   
   N2=S2+E2+A2+I2+R2+H2+P2

   for(tx in time2)
   {
      Smat=matrix(as.numeric(S[tx,]/N[tx,]),nrow=nages,ncol=nages)
      rovertime[tx+1] = rzerovar[tx]*max(Re(eigen((iperiod*Smat*betamatrix2))$values))
   }
}
N=S+E+A+I+R+H+P


sus=rowSums(S)
inf=rowSums(I)+rowSums(A)

propS=rowSums(S)/rowSums(N)
immune=rowSums(R)/rowSums(N)
popsize=rowSums(N)
cumdeaths4=rowSums(D)[1:(length(times)-1)]

totdeaths4=c(0.1,0.1,diff(rowSums(D)[1:(length(times)-1)]))
inhosp=rowSums(H) + rowSums(P)

mycols=rainbow(length(breakpoints)+1,alpha=0.3)
bp3=c(-100,-100,rep(breakpoints,each=4),SIMTIME+200,SIMTIME+200)+as.Date("23/03/2020",format="%d/%m/%Y")-march23
y1=rep(c(0,1000,1000,0),times=7)

#pdf("figs/model1.pdf",width=10,height=7)
par(mfrow=c(2,2),mar=c(4,5,1,1))
x1=times[1:(SIMTIME-1)]+as.Date("23/03/2020",format="%d/%m/%Y")-march30
y1=rovertime[1:(SIMTIME-1)]

plot(x1,y1,type="l",lwd=3,cex.lab=1.5,cex.axis=1.3,xlab="",ylab=expression(R[eff]))
X1=times[breakpoints]+as.Date("23/03/2020",format="%d/%m/%Y")-march23
#axis(1,at=X1+1,labels=format(X1+1,"%b"),cex.axis=0.9,las=2)

for(r in 0:length(breakpoints)){polygon(x=bp3[(r*4+1):(4*(r+1))],y=c(-200,max(y1)*1.5,max(y1)*1.5,-200),col=mycols[r+1],border=F)}

lines(x1[1:(SIMTIME-1)],y1,type="l",lwd=3,cex.lab=1.5,cex.axis=1.3)

abline(h=1,lty=2)
legend('top',c(expression(R[eff])),bty="n",cex=1.5)
#legend('topright',paste(c(rzero,rtstages[1:4])),pt.bg=mycols,pch=22,
#       title=expression(R[t]),bg="white",cex=1.3,pt.cex=2)


plot(x1,inhosp[1:(SIMTIME-1)],col=1,type="l",lwd=3,cex.lab=1.5,cex.axis=1.3,xlab="",ylab="# in hospital")
X1=times[breakpoints]+as.Date("23/03/2020",format="%d/%m/%Y")-march23
#axis(1,at=X1+1,labels=format(X1+1,"%d %b"),cex.axis=0.9,las=2)

for(r in 0:length(breakpoints)){polygon(x=bp3[(r*4+1):(4*(r+1))],y=c(-200000,max(inhosp)*1.5,max(inhosp)*1.5,-200000),col=mycols[r+1],border=F)}

lines(x1,inhosp[1:(SIMTIME-1)],col=1,type="l",lwd=3)

legend('topright',c('Hospital\noccupancy'),bty="n",cex=1.5)

plot(x1,totdeaths4[1:(SIMTIME-1)],col="black",type="l",lwd=3,
     ylim=range(totdeaths4),ylab="Deaths per day",xlab="",cex.lab=1.5,cex.axis=1.3)
X1=times[breakpoints]+as.Date("23/03/2020",format="%d/%m/%Y")-march23

axis(1,at=X1+1,labels=format(X1+1,"%b"),cex.axis=0.9,las=2)

for(r in 0:length(breakpoints)){
   polygon(x=bp3[(r*4+1):(4*(r+1))],y=c(-200,max(totdeaths4)*1.5,max(totdeaths4)*1.5,-200),col=mycols[r+1],border=F)
   }

lines(x1,totdeaths4[1:(SIMTIME-1)],col="black",type="l",lwd=3,
      ylim=range(totdeaths4,1100),ylab="Deaths per day")

legend('topright',c('Deaths'),bty="n",cex=1.5)

#number of people infected
plot(x1,inf[1:(SIMTIME-1)],type="l",lwd=3,ylab="# people infected",xlab="",cex.lab=1.5,cex.axis=1.3)
X1=times[breakpoints]+as.Date("23/03/2020",format="%d/%m/%Y")-march23

#axis(1,at=X1+1,labels=format(X1+1,"%d %b"),cex.axis=0.9,las=2)

for(r in 0:length(breakpoints)){polygon(x=bp3[(r*4+1):(4*(r+1))],y=c(-200000000,max(inf)*1.5,max(inf)*1.5,-200000000),col=mycols[r+1],border=F)}

lines(x1,inf[1:(SIMTIME-1)],type="l",lwd=3,ylab="# people infected",xlab="",cex.lab=1.5,cex.axis=1.3)
legend('topright',c('Number\ninfected'),bty="n",cex=1.5)


#dev.off()
