source("setup.R")   # library(seacarb) + the buffer functions from ../R
set.seed(42)
rich<-function(f,x,e){d1<-(f(x+e)-f(x-e))/(2*e);d2<-(f(x+2*e)-f(x-2*e))/(4*e);(4*d1-d2)/3}
N<-300
Tv<-runif(N,-1.8,32); Sv<-runif(N,28,38)
ATv<-runif(N,2000,2450)*1e-6; CTv<-ATv*runif(N,0.83,0.99)
Pv<-runif(N,0,2.5)*1e-6; Siv<-runif(N,0,120)*1e-6
bd<-function(S,T,CT,AT,P,Si) buffderiv(15,AT,CT,S=S,T=T,Pt=P,Sit=Si,
     k1k2="l",kf="dg",ks="d",pHscale="T",b="u74",warn="n")
worst<-c(dBe_dT=0,dBe_dS=0,dBe_dCT=0,dBe_dAT=0); wres<-0; nbad<-0
for(i in 1:N){
  a<-tryCatch(bd(Sv[i],Tv[i],CTv[i],ATv[i],Pv[i],Siv[i]),error=function(e)NULL)
  if(is.null(a)||!is.finite(a$Be)){nbad<-nbad+1; next}
  wres<-max(wres,abs(a$alk_residual)/ATv[i])
  f<-function(S,T,CT,AT) bd(S,T,CT,AT,Pv[i],Siv[i])$Be
  nu<-c(dBe_dT =rich(function(x) f(Sv[i],x,CTv[i],ATv[i]),Tv[i],1e-3),
        dBe_dS =rich(function(x) f(x,Tv[i],CTv[i],ATv[i]),Sv[i],1e-3),
        dBe_dCT=rich(function(x) f(Sv[i],Tv[i],x,ATv[i]),CTv[i],CTv[i]*1e-5),
        dBe_dAT=rich(function(x) f(Sv[i],Tv[i],CTv[i],x),ATv[i],ATv[i]*1e-5))
  for(k in names(worst)) worst[k]<-max(worst[k],abs(100*(a[[k]]-nu[[k]])/nu[[k]]))
}
cat(sprintf("N=%d  failed=%d\n",N,nbad))
cat(sprintf("worst |alk residual|/AT : %.2e\n",wres))
for(k in names(worst)) cat(sprintf("worst %-8s vs FD : %.2e %%\n",k,worst[k]))
cat(sprintf("\nOVERALL WORST: %.2e %%  -> %s\n",max(worst),
    ifelse(max(worst)<1e-4,"PASS","FAIL")))
