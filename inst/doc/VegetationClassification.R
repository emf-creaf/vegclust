### R code from vignette source 'VegetationClassification.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: VegetationClassification.Rnw:21-22
###################################################
options(width=67)


###################################################
### code chunk number 2: VegetationClassification.Rnw:27-28
###################################################
library(vegclust)


###################################################
### code chunk number 3: VegetationClassification.Rnw:33-35
###################################################
data(wetland)
dim(wetland)


###################################################
### code chunk number 4: VegetationClassification.Rnw:38-39
###################################################
wetlandchord = decostand(wetland,"normalize")


###################################################
### code chunk number 5: VegetationClassification.Rnw:42-43
###################################################
dchord = dist(wetlandchord)


###################################################
### code chunk number 6: VegetationClassification.Rnw:178-180
###################################################
wetland.km = vegclust(x = wetlandchord, mobileCenters=3, 
                      method="KM", nstart=20)


###################################################
### code chunk number 7: VegetationClassification.Rnw:183-184
###################################################
names(wetland.km)


###################################################
### code chunk number 8: VegetationClassification.Rnw:187-188
###################################################
t(wetland.km$memb)


###################################################
### code chunk number 9: VegetationClassification.Rnw:191-192
###################################################
round(wetland.km$mobileCenters, dig=3)


###################################################
### code chunk number 10: VegetationClassification.Rnw:195-198
###################################################
wetland.kmdist = vegclustdist(x = dchord, mobileMemb=3, 
                              method="KM", nstart = 20)
names(wetland.kmdist)


###################################################
### code chunk number 11: VegetationClassification.Rnw:201-202
###################################################
wetland.kmdist$mobileCenters


###################################################
### code chunk number 12: VegetationClassification.Rnw:205-206
###################################################
t(wetland.kmdist$memb)


###################################################
### code chunk number 13: VegetationClassification.Rnw:209-211
###################################################
wetland.km$mode
wetland.kmdist$mode


###################################################
### code chunk number 14: VegetationClassification.Rnw:217-218
###################################################
round(t(wetland.km$dist2clusters), dig=2)


###################################################
### code chunk number 15: VegetationClassification.Rnw:221-224
###################################################
wetland.fcm = vegclust(x = wetlandchord, mobileCenters=3, 
                       method="FCM", m=1.2, nstart=20)
round(t(wetland.fcm$memb), dig=3)


###################################################
### code chunk number 16: VegetationClassification.Rnw:229-232
###################################################
groups = defuzzify(wetland.fcm)$cluster
groups
table(groups)


###################################################
### code chunk number 17: VegetationClassification.Rnw:235-238
###################################################
groups = defuzzify(wetland.fcm, method = "cut", alpha = 0.8)$cluster
groups
table(groups, useNA = "always")


###################################################
### code chunk number 18: VegetationClassification.Rnw:241-244
###################################################
wetland.fcm2 = vegclust(x = wetlandchord, mobileCenters=3, 
                       method="FCM", m=10, nstart=20)
round(t(wetland.fcm2$memb), dig=3)


###################################################
### code chunk number 19: VegetationClassification.Rnw:247-249
###################################################
groups2 = defuzzify(wetland.fcm2, method = "cut", alpha = 0.8)$cluster
table(groups2, useNA = "always")


###################################################
### code chunk number 20: VegetationClassification.Rnw:255-258
###################################################
wetland.nc = vegclust(x = wetlandchord, mobileCenters=3,
                       method="NC", m=1.2, dnoise=0.8, nstart=20)
round(t(wetland.nc$memb), dig=2)


###################################################
### code chunk number 21: VegetationClassification.Rnw:261-264
###################################################
groups = defuzzify(wetland.nc)$cluster
groups
table(groups)


###################################################
### code chunk number 22: VegetationClassification.Rnw:267-270
###################################################
groups = defuzzify(wetland.nc, method="cut", alpha=0.8)$cluster
groups
table(groups, useNA = "always")


###################################################
### code chunk number 23: VegetationClassification.Rnw:275-278
###################################################
dist(wetland.km$mobileCenters)
dist(wetland.fcm$mobileCenters)
dist(wetland.nc$mobileCenters)


###################################################
### code chunk number 24: VegetationClassification.Rnw:284-287
###################################################
wetland.kmdd = vegclust(x = wetlandchord, mobileCenters=3, 
                      method="KMdd", nstart=20)
t(wetland.kmdd$memb)


###################################################
### code chunk number 25: VegetationClassification.Rnw:290-291
###################################################
round(wetland.kmdd$mobileCenters, dig=3)


###################################################
### code chunk number 26: VegetationClassification.Rnw:294-297
###################################################
wetland.kmdd = vegclustdist(x = dchord, mobileMemb=3, 
                      method="KMdd", nstart=20)
wetland.kmdd$mobileCenters


###################################################
### code chunk number 27: VegetationClassification.Rnw:302-308
###################################################
wetland.31 = wetlandchord[1:31,]
wetland.31 = wetland.31[,colSums(wetland.31)>0]
dim(wetland.31)
wetland.10 = wetlandchord[-(1:31),]
wetland.10 = wetland.10[,colSums(wetland.10)>0] 
dim(wetland.10)


###################################################
### code chunk number 28: VegetationClassification.Rnw:311-314
###################################################
km = kmeans(wetland.31, 2)
groups = km$cluster
groups


###################################################
### code chunk number 29: VegetationClassification.Rnw:317-318
###################################################
wetland.31.km = as.vegclust(wetland.31, groups)


###################################################
### code chunk number 30: VegetationClassification.Rnw:321-322
###################################################
wetland.31.km$method


###################################################
### code chunk number 31: VegetationClassification.Rnw:325-327
###################################################
wetland.10.km = vegclass(wetland.31.km, wetland.10)
defuzzify(wetland.10.km)$cluster


###################################################
### code chunk number 32: VegetationClassification.Rnw:330-331
###################################################
wetland.31.km.d = as.vegclust(dist(wetland.31), groups)


###################################################
### code chunk number 33: VegetationClassification.Rnw:334-335
###################################################
wetland.d.10.31 = as.data.frame(as.matrix(dchord)[32:41,1:31])


###################################################
### code chunk number 34: VegetationClassification.Rnw:338-340
###################################################
wetland.d.11.km = vegclass(wetland.31.km.d,wetland.d.10.31)
defuzzify(wetland.d.11.km)$cluster


###################################################
### code chunk number 35: VegetationClassification.Rnw:344-348
###################################################
wetland.31.nc = as.vegclust(wetland.31, groups, method="HNC", 
                            dnoise = 0.8)
wetland.10.nc = vegclass(wetland.31.nc, wetland.10)
defuzzify(wetland.10.nc)$cluster


###################################################
### code chunk number 36: VegetationClassification.Rnw:357-362
###################################################
cf = conformveg(wetland.31, wetland.10)
wetland.31.cf<- cf$x
wetland.10.cf<- cf$y
dim(wetland.31.cf)
dim(wetland.10.cf)


###################################################
### code chunk number 37: VegetationClassification.Rnw:368-369
###################################################
fixed = clustcentroid(wetland.31.cf, groups)


###################################################
### code chunk number 38: VegetationClassification.Rnw:375-380
###################################################
wetland.nc = vegclust(wetland.10.cf, mobileCenters=1, 
                      fixedCenters = fixed, 
                      method = wetland.31.nc$method,
                      dnoise=wetland.31.nc$dnoise, nstart=10)
defuzzify(wetland.nc)$cluster


###################################################
### code chunk number 39: VegetationClassification.Rnw:385-390
###################################################
wetland.km = vegclust(wetland.10.cf, mobileCenters=1, 
                      fixedCenters = fixed, 
                      method = "KM",
                      nstart=10)
defuzzify(wetland.km)$cluster


###################################################
### code chunk number 40: VegetationClassification.Rnw:397-402
###################################################
wetland.nc = vegclust(rbind(wetland.31.cf,wetland.10.cf), mobileCenters=1, 
                      fixedCenters = fixed, 
                      method = wetland.31.nc$method,
                      dnoise=wetland.31.nc$dnoise, nstart=10)
defuzzify(wetland.nc)$cluster


###################################################
### code chunk number 41: VegetationClassification.Rnw:408-409
###################################################
fixedDist = wetland.d.11.km$dist2clusters


###################################################
### code chunk number 42: VegetationClassification.Rnw:412-417
###################################################
wetland.km.d = vegclustdist(dist(wetland.10), mobileMemb = 1,
                            fixedDistToCenters=fixedDist, 
                            method = "KM",
                            nstart=10)
defuzzify(wetland.km.d)$cluster


###################################################
### code chunk number 43: VegetationClassification.Rnw:421-422
###################################################
fixedDist = rbind(wetland.31.km.d$dist2clusters, wetland.d.11.km$dist2clusters)


###################################################
### code chunk number 44: VegetationClassification.Rnw:425-430
###################################################
wetland.km.d = vegclustdist(dchord, mobileMemb = 1,
                            fixedDistToCenters=fixedDist, 
                            method = "KM",
                            nstart=10)
defuzzify(wetland.km.d)$cluster


###################################################
### code chunk number 45: VegetationClassification.Rnw:437-438
###################################################
groups = c(rep(1, 17), rep(2, 14), rep(3,10))


###################################################
### code chunk number 46: VegetationClassification.Rnw:445-447
###################################################
centroids = clustcentroid(wetlandchord, groups)
round(centroids, dig=3)


###################################################
### code chunk number 47: VegetationClassification.Rnw:451-453
###################################################
medoids = clustmedoid(wetlandchord, groups)
medoids


###################################################
### code chunk number 48: VegetationClassification.Rnw:464-465
###################################################
clustvar(wetlandchord, groups)


###################################################
### code chunk number 49: VegetationClassification.Rnw:477-478
###################################################
clustvar(dchord, groups)


###################################################
### code chunk number 50: VegetationClassification.Rnw:481-482
###################################################
clustvar(wetlandchord)


###################################################
### code chunk number 51: VegetationClassification.Rnw:487-488
###################################################
as.dist(as.matrix(dchord)[medoids,medoids])


###################################################
### code chunk number 52: VegetationClassification.Rnw:492-493
###################################################
dist(centroids)


###################################################
### code chunk number 53: VegetationClassification.Rnw:500-501
###################################################
interclustdist(dchord,groups)


###################################################
### code chunk number 54: VegetationClassification.Rnw:508-509
###################################################
c = clustconst(wetlandchord, memb = as.memb(groups))


###################################################
### code chunk number 55: VegetationClassification.Rnw:513-514
###################################################
d=summary(c, mode="all")


###################################################
### code chunk number 56: VegetationClassification.Rnw:518-519
###################################################
summary(c, mode="cluster", name=names(c)[1])


