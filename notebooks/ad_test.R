source("setup.R")   # library(seacarb) + the buffer functions from ../R
source("ad_be.R")   # AD implementation
cases <- data.frame(
  lab=c("polar","temperate","notebook_base","subtropical","warm_salty_gyre","southern_ocean"),
  T=c(-1,10,20,25,29,5), S=c(34,34.8,35,35,36.5,33.5),
  CT=c(2150,2050,2000,2000,1950,2200)*1e-6, AT=c(2300,2290,2300,2300,2380,2320)*1e-6,
  Pt=c(1.5,.5,.2,.1,.05,1.8)*1e-6, Sit=c(60,10,3,2,1,70)*1e-6, stringsAsFactors=FALSE)

cat("AD (dual numbers + implicit function theorem)  vs  hand-derived analytic (buffderiv)\n\n")
cat(sprintf("%-16s %14s %14s %14s %14s\n","case","dBe/dT","dBe/dS","dBe/dCT","dBe/dAT"))
worst<-0
for(i in 1:nrow(cases)){
  a <- buffderiv(15, cases$AT[i], cases$CT[i], S=cases$S[i], T=cases$T[i],
                 Pt=cases$Pt[i], Sit=cases$Sit[i], k1k2="l", kf="dg", ks="d",
                 pHscale="T", b="u74", warn="n")
  h0 <- 10^(-carb(15,cases$AT[i],cases$CT[i],S=cases$S[i],T=cases$T[i],P=0,
                  Pt=cases$Pt[i],Sit=cases$Sit[i],k1k2="l",kf="dg",ks="d",
                  pHscale="T",b="u74",warn="n")$pH)
  g <- be_ad(cases$T[i],cases$S[i],cases$CT[i],cases$AT[i],cases$Pt[i],cases$Sit[i],h0)
  rd <- c(100*(g$dBe_dT-a$dBe_dT)/a$dBe_dT, 100*(g$dBe_dS-a$dBe_dS)/a$dBe_dS,
          100*(g$dBe_dCT-a$dBe_dCT)/a$dBe_dCT, 100*(g$dBe_dAT-a$dBe_dAT)/a$dBe_dAT)
  worst<-max(worst,max(abs(rd)))
  cat(sprintf("%-16s %14.2e %14.2e %14.2e %14.2e   (%% diff)\n",cases$lab[i],rd[1],rd[2],rd[3],rd[4]))
  cat(sprintf("%-16s Be: AD=%.10f  analytic=%.10f\n","",g$Be,a$Be))
}
cat(sprintf("\nWORST AD vs analytic: %.2e %%\n\n", worst))

## ---- cost -----------------------------------------------------------------
N<-2000
S<-runif(N,30,37); T<-runif(N,-1.8,30); AT<-runif(N,2100,2400)*1e-6
CT<-AT*runif(N,.85,.97); Pt<-runif(N,0,2.5)*1e-6; Sit<-runif(N,0,120)*1e-6
t_carb <- system.time(cb <- carb(15,AT,CT,S=S,T=T,P=0,Pt=Pt,Sit=Sit,k1k2="l",kf="dg",
                     ks="d",pHscale="T",b="u74",warn="n"))[["elapsed"]]
h0v <- 10^(-cb$pH)
t_ana <- system.time(buffderiv(15,AT,CT,S=S,T=T,Pt=Pt,Sit=Sit,k1k2="l",kf="dg",ks="d",
                     pHscale="T",b="u74",warn="n"))[["elapsed"]]
t_ad  <- system.time(for(i in 1:N) be_ad(T[i],S[i],CT[i],AT[i],Pt[i],Sit[i],h0v[i]))[["elapsed"]]
cat(sprintf("cost for N = %d points\n", N))
cat(sprintf("  carb() alone (the solve)            : %6.2f s\n", t_carb))
cat(sprintf("  buffderiv (carb + analytic derivs)  : %6.2f s\n", t_ana))
cat(sprintf("  AD, scalar dual numbers, given h    : %6.2f s   (%.0fx the analytic route)\n",
            t_ad, t_ad/t_ana))
