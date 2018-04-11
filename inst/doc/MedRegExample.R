### R code from vignette source 'MedRegExample.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: MedRegExample.Rnw:21-22
###################################################
options(width=67)


###################################################
### code chunk number 2: MedRegExample.Rnw:26-27
###################################################
library(vegclust)


###################################################
### code chunk number 3: MedRegExample.Rnw:32-35
###################################################
data(medreg)
class(medreg)
length(medreg)


###################################################
### code chunk number 4: MedRegExample.Rnw:38-39
###################################################
strataUp = c(20,50,100,300,600,1200,2400)


###################################################
### code chunk number 5: MedRegExample.Rnw:42-43
###################################################
strataWidths = c(20,30,50,200,300,600,1200)


###################################################
### code chunk number 6: MedRegExample.Rnw:46-47
###################################################
medreg[[1]]


###################################################
### code chunk number 7: MedRegExample.Rnw:52-53
###################################################
medreg.CAP <- CAP(medreg)


###################################################
### code chunk number 8: MedRegExample.Rnw:56-58
###################################################
class(medreg.CAP)
length(medreg.CAP)


###################################################
### code chunk number 9: MedRegExample.Rnw:61-62
###################################################
medreg.CAP[[1]]


###################################################
### code chunk number 10: MedRegExample.Rnw:66-71
###################################################
plot(medreg.CAP, plots="1", sizes=strataUp, xlab="Height (cm)", 
     ylab="Cumulative percent cover")
legend("topright", col=1:5, lty=1, 
       legend=c("Pines","Oaks","Tall shrubs","Scrubs","Grass"), 
       bty="n")


###################################################
### code chunk number 11: MedRegExample.Rnw:79-84
###################################################
pl = rep(1,100) # All trees in the same plot
sp = ifelse(runif(100)>0.5,1,2) # Random species identity (species 1 or 2)
h=rgamma(100,10,2) # Heights (m)
d = rpois(100, lambda=h^2) # Diameters (cm)
m = data.frame(plot=pl,species=sp, height=h,diameter=d) 


###################################################
### code chunk number 12: MedRegExample.Rnw:87-88
###################################################
m$ba = pi*(m$diameter/200)^2


###################################################
### code chunk number 13: MedRegExample.Rnw:91-92
###################################################
print(head(m))


###################################################
### code chunk number 14: MedRegExample.Rnw:95-97
###################################################
heights = seq(0,4, by=.25)^2 # Quadratic classes
diams = seq(0,130, by=5) # Linear classes


###################################################
### code chunk number 15: MedRegExample.Rnw:100-104
###################################################
tree.S<-stratifyvegdata(m, sizes1=heights, sizes2=diams, 
                   plotColumn = "plot", speciesColumn = "species", 
                   size1Column = "height", size2Column = "diameter", 
                   abundanceColumn = "ba")


###################################################
### code chunk number 16: MedRegExample.Rnw:107-108
###################################################
tree.CAS <- CAS(tree.S)


###################################################
### code chunk number 17: MedRegExample.Rnw:112-119
###################################################
par(mfrow=c(2,1), mar=c(4,5,2,1))
plot(tree.CAS, species=1, sizes1=heights[-1], xlab="height (m)", 
     ylab="diameter (cm)", sizes2=diams[-1], zlab="Basal area (m2)",
     zlim = c(0,6), main="Species 1")
plot(tree.CAS, species=2, sizes1=heights[-1], xlab="height (m)", 
     ylab="diameter (cm)", sizes2=diams[-1], zlab="Basal area (m2)",
     zlim = c(0,6), main = "Species 2")


###################################################
### code chunk number 18: MedRegExample.Rnw:123-124
###################################################
print(CASmargin(tree.CAS, margin=1))


###################################################
### code chunk number 19: MedRegExample.Rnw:127-131
###################################################
tree.S2<-stratifyvegdata(m, sizes1=heights, plotColumn = "plot", 
                         speciesColumn = "species", size1Column = "height", 
                         abundanceColumn = "ba")
print(CAP(tree.S2))


###################################################
### code chunk number 20: MedRegExample.Rnw:134-139
###################################################
par(mfrow=c(2,1), mar=c(4,5,2,1))
plot(CASmargin(tree.CAS,margin=1), plots=1, sizes=heights[-1], 
     xlab="height (m)", ylab="Basal area (m2)", ylim = c(0,7))
plot(CASmargin(tree.CAS,margin=2), plots=1, sizes=diams[-1], 
     xlab="diameter (cm)", ylab="Basal area (m2)", ylim = c(0,7))


###################################################
### code chunk number 21: MedRegExample.Rnw:143-145
###################################################
medreg.D = vegdiststruct(medreg.CAP, method="bray", 
                         classWeights=strataWidths)


###################################################
### code chunk number 22: MedRegExample.Rnw:148-149
###################################################
as.matrix(medreg.D)[1,2]


###################################################
### code chunk number 23: MedRegExample.Rnw:152-154
###################################################
medreg.Dsqrt = vegdiststruct(medreg.CAP, method="bray", 
                         classWeights=strataWidths, transform="sqrt")


###################################################
### code chunk number 24: MedRegExample.Rnw:158-165
###################################################
par(mfrow=c(2,1), mar=c(4,5,2,1))
X<-cmdscale(medreg.D, k=2)
plot(X, xlab="MDS 1", ylab="MDS 2", asp=1,
     main="Cover untransformed", cex=0.5)
Xsqrt<-cmdscale(medreg.Dsqrt, k=2)
plot(Xsqrt, xlab="MDS 1", ylab="MDS 2", asp=1,
     main="Cover sqrt-transformed", cex=0.5)


###################################################
### code chunk number 25: MedRegExample.Rnw:172-174
###################################################
nclusters = 6
dnoise = 0.40


###################################################
### code chunk number 26: MedRegExample.Rnw:177-179
###################################################
vc<-vegclustdist(medreg.Dsqrt, mobileMemb = nclusters, 
                 method="HNCdd", dnoise=dnoise, nstart=100)


###################################################
### code chunk number 27: MedRegExample.Rnw:182-184
###################################################
medoids<-vc$mobileCenters
print(medoids)


###################################################
### code chunk number 28: MedRegExample.Rnw:187-189
###################################################
cluster<-defuzzify(vc)$cluster
table(cluster)


###################################################
### code chunk number 29: MedRegExample.Rnw:193-198
###################################################
clNum = as.numeric(as.factor(cluster))
plot(Xsqrt, xlab="MDS 1", ylab="MDS 2", 
     pch=clNum, col=clNum)
legend("topleft", col=1:(nclusters+1), pch=1:(nclusters+1),
       legend=levels(as.factor(cluster)), bty="n")


###################################################
### code chunk number 30: MedRegExample.Rnw:204-206
###################################################
CAPm = CAPcenters(medreg.CAP, vc)
n = names(CAPm)


###################################################
### code chunk number 31: MedRegExample.Rnw:209-210
###################################################
round(CAPm[[n[4]]], dig=1)


###################################################
### code chunk number 32: MedRegExample.Rnw:214-225
###################################################
par(mfrow=c(3,2), mar=c(4,4,3,0))
plot(CAPm, plots=n[1], sizes = strataWidths, 
     ylab="Percent cover", main="M1")
plot(CAPm, plots=n[2], sizes = strataWidths, main="M2")
plot(CAPm, plots=n[3], sizes = strataWidths,  
     ylab="Percent cover", main="M3")
plot(CAPm, plots=n[4], sizes = strataWidths, main="M4")
plot(CAPm, plots=n[5], sizes = strataWidths, 
     xlab="Height (cm)", ylab="Percent cover", main="M5")
plot(CAPm, plots=n[6], sizes = strataWidths, 
     xlab="Height (cm)",  main="M6")


