library(data.table)

## yogyakarta data
data <- fread("./yogyakarta_seropositivity.csv")
load("foi_and_r0.RData")

foi.fn <- function(foi,ages=data$age_high) 1-exp(-foi*ages)
fit.fn <- function(par) sum(abs(data$seropos_mean/100-foi.fn(par)))


foi.bound <- c(low=0,upp=2)
optim.out <- optim(c(foi=foi.val),fit.fn,method="Brent",
                   lower=foi.bound["low"],upper=foi.bound["upp"])

pdf("foi_comparison.pdf")
curve(foi.fn(optim.out$par,x),from=0,to=11,xaxt="n",xlab="Age",ylab="Seropositivity",
      lwd=2,bty="n",ylim=c(0,1))
points(data$age_high,data$seropos_mean/100)
for (ii in 1:nrow(data)){
    lines(rep(data$age_high[ii],2),c(data$seropos_low[ii],data$seropos_high[ii])/100)
}
axis(1,at=data$age_high,tick=FALSE,
     labels=paste0(data$age_low,"-",data$age_high," years"))
axis(1,at=c(data$age_low,max(data$age_low)+ 2),
     labels=FALSE)
curve(foi.fn(4*foi.val,x),from=0,to=11,add=TRUE,col="red",lty=1,lwd=2)
legend("topleft",legend=paste0(c("Best fit (","Cattarino ("),
                             round(c(optim.out$par,foi.val*4),3),
                             ")"),
       lty=1,lwd=2,col=c("black","red"),bty="n")
dev.off()

## Save
foi.yogya <- optim.out$par / 4
save(foi.yogya,file="foi_yogya.RData")
