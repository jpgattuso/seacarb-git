source("setup.R")   # library(seacarb) + the buffer functions from ../R
source("ad_be.R")   # AD implementation
set.seed(1); N<-2000
S<-runif(N,30,37); T<-runif(N,-1.8,30); AT<-runif(N,2100,2400)*1e-6
CT<-AT*runif(N,.85,.97); Pt<-runif(N,0,2.5)*1e-6; Sit<-runif(N,0,120)*1e-6
t_carb<-system.time(cb<-carb(15,AT,CT,S=S,T=T,P=0,Pt=Pt,Sit=Sit,k1k2="l",kf="dg",ks="d",
                             pHscale="T",b="u74",warn="n"))[["elapsed"]]
h0<-10^(-cb$pH)
t_ana<-system.time(a<-buffderiv(15,AT,CT,S=S,T=T,Pt=Pt,Sit=Sit,k1k2="l",kf="dg",ks="d",
                                pHscale="T",b="u74",warn="n"))[["elapsed"]]
t_ad <-system.time(g<-be_ad(T,S,CT,AT,Pt,Sit,h0))[["elapsed"]]
w<-max(abs(100*(g$dBe_dT-a$dBe_dT)/a$dBe_dT), abs(100*(g$dBe_dS-a$dBe_dS)/a$dBe_dS),
       abs(100*(g$dBe_dCT-a$dBe_dCT)/a$dBe_dCT), abs(100*(g$dBe_dAT-a$dBe_dAT)/a$dBe_dAT))
cat(sprintf("vectorised AD vs analytic, %d random points: worst %.2e %%\n\n",N,w))
cat(sprintf("  carb() alone (the nonlinear solve)  : %6.3f s\n",t_carb))
cat(sprintf("  buffderiv  (carb + analytic derivs) : %6.3f s\n",t_ana))
cat(sprintf("  AD (vectorised duals), h given      : %6.3f s\n",t_ad))
cat(sprintf("  -> AD derivative cost is %.1fx the analytic derivative cost,\n",t_ad/(t_ana-t_carb+1e-9)))
cat(sprintf("     but only %.2fx the cost of the carb() solve you must pay anyway.\n",t_ad/t_carb))
