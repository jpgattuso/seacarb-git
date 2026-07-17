source("setup.R")   # library(seacarb) + the buffer functions from ../R
cases <- data.frame(
  T=c(-1,2,10,20,25,29,28,15,25,5), S=c(34,34.5,34.8,35,35,36.5,33,30,35,33.5),
  CT=c(2150,2100,2050,2000,2000,1950,1900,1950,2100,2200)*1e-6,
  AT=c(2300,2280,2290,2300,2300,2380,2150,2100,2300,2320)*1e-6,
  Pt=c(1.5,1.2,.5,.2,.1,.05,.3,.8,.1,1.8)*1e-6,
  Sit=c(60,40,10,3,2,1,5,15,2,70)*1e-6,
  case=c("polar","subpolar","temperate","notebook_base","subtropical",
         "warm_salty_gyre","warm_fresh","brackish_coastal","high_CO2_2100",
         "southern_ocean"), stringsAsFactors=FALSE)
r <- buffderiv(15, cases$AT, cases$CT, S=cases$S, T=cases$T, Pt=cases$Pt,
               Sit=cases$Sit, k1k2="l", kf="dg", ks="d", pHscale="T",
               b="u74", warn="n")
out <- cbind(case=cases$case, T=cases$T, S=cases$S,
             CT_umol=cases$CT*1e6, AT_umol=cases$AT*1e6,
             Pt_umol=cases$Pt*1e6, Sit_umol=cases$Sit*1e6,
             pH=-log10(-1/r$dh_dAT*0+10^0)*0,  # placeholder replaced below
             r)
out$pH <- NULL
out <- cbind(out, pH=as.numeric(carb(15,cases$AT,cases$CT,S=cases$S,T=cases$T,
             Pt=cases$Pt,Sit=cases$Sit,k1k2="l",kf="dg",ks="d",pHscale="T",
             b="u74",warn="n")$pH))
write.csv(out, "reference_values.csv", row.names=FALSE)
print(format(out[,c("case","pH","Be","dBe_dT","dBe_dS","dBe_dCT","dBe_dAT")],
             digits=10))
