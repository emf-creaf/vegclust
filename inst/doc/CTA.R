## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load libraries, echo = T--------------------------------------------
library(vegclust)
library(RColorBrewer)

## ------------------------------------------------------------------------
#Description of sites and surveys
sites = c(1,1,1,2,2,2)
surveys=c(1,2,3,1,2,3)

## ------------------------------------------------------------------------
#Raw data table
xy<-matrix(0, nrow=6, ncol=2)
xy[2,2]<-1
xy[3,2]<-2
xy[4:6,1] <- 0.5
xy[4:6,2] <- xy[1:3,2]
xy[6,1]<-1
cbind(sites,surveys,xy)

## ------------------------------------------------------------------------
#Distance matrix
D = dist(xy)
D

## ----pcoa, fig = TRUE, fig.height=4, fig.width=5, fig.align = "center"----
trajectoryPCoA(D, sites, surveys, traj.colors = c("black","red"), lwd = 2)

## ------------------------------------------------------------------------
trajectoryLengths(D, sites, surveys)
trajectoryAngles(D, sites, surveys)

## ------------------------------------------------------------------------
trajectoryProjection(D, 1:3, 1:3)
trajectoryProjection(D, 4:6, 4:6)

## ------------------------------------------------------------------------
trajectoryProjection(D, 3, 4:6)

## ------------------------------------------------------------------------
segmentDistances(D, sites, surveys)$Dseg

## ------------------------------------------------------------------------
trajectoryDistances(D, sites, surveys, distance.type = "Hausdorff")
trajectoryDistances(D, sites, surveys, distance.type = "DSPD")

## ----load avoca, echo=T--------------------------------------------------
data("avoca")

## ----distance, echo=TRUE-------------------------------------------------
avoca_D_man = vegdiststruct(avoca_strat, method="manhattan", transform = function(x){log(x+1)})

## ----avoca pcoa, echo=T, fig=TRUE, fig.height=6, fig.width=6, fig.align = "center"----
trajectoryPCoA(avoca_D_man,  avoca_sites, avoca_surveys,
               traj.colors = brewer.pal(8,"Accent"), 
               axes=c(1,2), length=0.1, lwd=2)
legend("topright", bty="n", legend = 1:8, col = brewer.pal(8,"Accent"), lwd=2)


## ----trajectory lengths, echo=T------------------------------------------
trajectoryLengths(avoca_D_man, avoca_sites, avoca_surveys)

## ----int1, echo=FALSE----------------------------------------------------
plotTrajDiamDist<-function(cli = 7) {
l = colnames(avoca_strat[[1]])
ncl = 14
m197072= avoca_strat[avoca_surveys==1][[cli]]["NOTCLI",2:ncl]
m197072[m197072<1] = NA
m1974 = avoca_strat[avoca_surveys==2][[cli]]["NOTCLI",2:ncl]
m1974[m1974<1] = NA
m1978 = avoca_strat[avoca_surveys==3][[cli]]["NOTCLI",2:ncl]
m1978[m1978<1] = NA
m1983 = avoca_strat[avoca_surveys==4][[cli]]["NOTCLI",2:ncl]
m1983[m1983<1] = NA
m1987 = avoca_strat[avoca_surveys==5][[cli]]["NOTCLI",2:ncl]
m1987[m1987<1] = NA
m1993 = avoca_strat[avoca_surveys==6][[cli]]["NOTCLI",2:ncl]
m1993[m1993<1] = NA
m1999 = avoca_strat[avoca_surveys==7][[cli]]["NOTCLI",2:ncl]
m1999[m1999<1] = NA
m2004 = avoca_strat[avoca_surveys==8][[cli]]["NOTCLI",2:ncl]
m2004[m2004<1] = NA
m2009 = avoca_strat[avoca_surveys==9][[cli]]["NOTCLI",2:ncl]
m2009[m2009<1] = NA


plot(m197072, type="l", ylim=c(1,200), log="y",
       xlab="", ylab="Number of individuals (log)", main=paste0("Trajectory ",cli), 
       axes=FALSE, col=gray(0.8), lwd=2)
axis(2, las=2)
axis(1, at=1:(ncl-1), labels=l[2:ncl], las=2)
lines(m1974, col=gray(0.7), lwd=2)
lines(m1978, col=gray(0.6), lwd=2)
lines(m1983, col=gray(0.5), lwd=2)
lines(m1987, col=gray(0.4), lwd=2)
lines(m1993, col=gray(0.3), lwd=2)
lines(m1999, col=gray(0.2), lwd=2)
lines(m2004, col=gray(0.1), lwd=2)
lines(m2009, col=gray(0), lwd=2)
legend("topright", bty="n", lwd=2,col=gray(seq(0.8,0, by=-0.1)), legend=c("1970/72","1974","1978","1983", "1987", "1993","1999","2004","2009"))
}

## ----trajectory angles, echo=T-------------------------------------------
trajectoryAngles(avoca_D_man, avoca_sites, avoca_surveys)

## ----trajectory 3 DBH dist, echo=F, fig.height=4, fig.width=8, fig.align = "center"----
par(mfrow=c(1,2))
trajectoryPCoA(avoca_D_man,  avoca_sites, avoca_surveys,
               selection= 3,
               length=0.1, lwd=2)
plotTrajDiamDist(3)

## ----trajectory 4, echo=T, fig.height=4, fig.width=8, fig.align = "center"----
par(mfrow=c(1,2))
trajectoryPCoA(avoca_D_man,  avoca_sites, avoca_surveys,
               selection= 4,
               length=0.1, lwd=2)
plotTrajDiamDist(4)

## ----avoca DT, echo=FALSE------------------------------------------------
avoca_D_traj_man = trajectoryDistances(avoca_D_man, avoca_sites, distance.type="DSPD", verbose=FALSE)
print(round(avoca_D_traj_man,3))

## ----avoca DT PCoA, echo=FALSE, fig = TRUE, fig.height=5, fig.width=6, fig.align="center"----
cmd_D2<-cmdscale(avoca_D_traj_man, add=TRUE, eig=TRUE, k=7)
x<-cmd_D2$points[,1]
y<-cmd_D2$points[,2]
plot(x,y, type="p", asp=1, xlab=paste0("PCoA 1 (", round(100*cmd_D2$eig[1]/sum(cmd_D2$eig)),"%)"), 
     ylab=paste0("PCoA 2 (", round(100*cmd_D2$eig[2]/sum(cmd_D2$eig)),"%)"), col="black",
     bg= brewer.pal(8,"Accent"), pch=21)
text(x,y, labels=1:8, pos=2)

