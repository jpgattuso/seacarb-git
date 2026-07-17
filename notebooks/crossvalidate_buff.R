source("setup.R")   # library(seacarb) + the buffer functions from ../R
cases <- data.frame(
 lab=c("no nutrients","Pt,Sit only","polar high Si","NH4t only","HSt only",
       "everything","everything cold","everything warm"),
 T  =c(20, 20, -1, 20, 20, 18, -1.5, 30),
 S  =c(35, 35, 34, 35, 35, 35, 34, 36),
 AT =c(2300,2300,2300,2300,2300,2300,2300,2380)*1e-6,
 CT =c(2000,2000,2150,2000,2000,2000,2150,1950)*1e-6,
 Pt =c(0, 0.2, 1.5, 0, 0, 2.5, 2.0, 0.1)*1e-6,
 Sit=c(0, 3.0, 60,  0, 0, 120, 100, 2)*1e-6,
 NH4=c(0, 0,   0,   3, 0, 3.0, 2.0, 0.5)*1e-6,
 HSt=c(0, 0,   0,   0, 2, 2.0, 1.5, 0.3)*1e-6, stringsAsFactors=FALSE)

cat("Be and Rf from all five routines. Rf = DIC/(CO2*Be), so Be is recovered from\n")
cat("buffesm's R and buffer's BetaD as DIC/(CO2*R).  All use k1k2='l', kf='dg'.\n\n")
cat(sprintf("%-17s %13s %11s %11s %11s %11s %11s\n","case","Be(buffsun)",
            "zwg","esm","buffer","buffderiv","max|dev|"))
worst<-0
for(i in 1:nrow(cases)){
 a<-list(flag=15,var1=cases$AT[i],var2=cases$CT[i],S=cases$S[i],T=cases$T[i],P=0,
         Pt=cases$Pt[i],Sit=cases$Sit[i],NH4t=cases$NH4[i],HSt=cases$HSt[i],
         k1k2="l",kf="dg",ks="d",pHscale="T",b="u74",warn="n")
 sp <- if (cases$NH4[i]!=0 || cases$HSt[i]!=0)
        do.call(carbfull,a[c("flag","var1","var2","S","T","P","Pt","Sit","NH4t","HSt","k1k2","kf","ks","pHscale","b")])
       else do.call(carb,a[c("flag","var1","var2","S","T","P","Pt","Sit","k1k2","kf","ks","pHscale","b","warn")])
 s<-do.call(buffsun,a); z<-do.call(buffzwg,a); e<-do.call(buffesm,a); u<-do.call(buffer,a)
 d<-buffderiv(15,cases$AT[i],cases$CT[i],S=cases$S[i],T=cases$T[i],Pt=cases$Pt[i],
              Sit=cases$Sit[i],NH4t=cases$NH4[i],HSt=cases$HSt[i],
              k1k2="l",kf="dg",ks="d",pHscale="T",b="u74",warn="n")
 Be<-c(sun=s$Be, zwg=z$Be, esm=sp$DIC/(sp$CO2*e$R), buf=sp$DIC/(sp$CO2*u$BetaD), der=d$Be)
 dev<-max(abs(Be-Be[["sun"]]))/Be[["sun"]]; worst<-max(worst,dev)
 cat(sprintf("%-17s %13.8f %11.8f %11.8f %11.8f %11.8f %11.1e\n", cases$lab[i],
     Be[["sun"]],Be[["zwg"]],Be[["esm"]],Be[["buf"]],Be[["der"]],dev))
}
cat(sprintf("\nWORST relative deviation in Be across all five routines: %.2e\n",worst))
cat("\nNow Rf:\n")
cat(sprintf("%-17s %13s %11s %11s %11s %11s\n","case","Rf(buffsun)","zwg RF","esm R","buffer BetaD","max|dev|"))
worstR<-0
for(i in 1:nrow(cases)){
 a<-list(flag=15,var1=cases$AT[i],var2=cases$CT[i],S=cases$S[i],T=cases$T[i],P=0,
         Pt=cases$Pt[i],Sit=cases$Sit[i],NH4t=cases$NH4[i],HSt=cases$HSt[i],
         k1k2="l",kf="dg",ks="d",pHscale="T",b="u74",warn="n")
 s<-do.call(buffsun,a); z<-do.call(buffzwg,a); e<-do.call(buffesm,a); u<-do.call(buffer,a)
 R<-c(s$Rf,z$RF,e$R,u$BetaD); dev<-max(abs(R-R[1]))/R[1]; worstR<-max(worstR,dev)
 cat(sprintf("%-17s %13.8f %11.8f %11.8f %11.8f %11.1e\n",cases$lab[i],R[1],R[2],R[3],R[4],dev))
}
cat(sprintf("\nWORST relative deviation in Rf across the four: %.2e\n",worstR))
