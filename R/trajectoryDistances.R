
.twoSegmentDistance<-function(dmat12, type="directed-segment") {
  ds1e1 = dmat12[1,2]
  ds1s2 = dmat12[1,3]
  ds1e2 = dmat12[1,4]
  de1s2 = dmat12[2,3]
  de1e2 = dmat12[2,4]
  ds2e2 = dmat12[3,4]
  if(type=="Hausdorff" || type == "directed-segment") {
    ps1_2 =.distanceToSegment(ds2e2,ds1s2, ds1e2)
    pe1_2 =.distanceToSegment(ds2e2,de1s2, de1e2)
    ps2_1 =.distanceToSegment(ds1e1,ds1s2, de1s2)
    pe2_1 =.distanceToSegment(ds1e1,ds1e2, de1e2)
    ds1_2 = ps1_2[3]
    de1_2 = pe1_2[3]
    ds2_1 = ps2_1[3]
    de2_1 = pe2_1[3]
    if(type == "directed-segment") { #Modifications for directionality of segments
      if(ps1_2[1]>pe1_2[1]) {
        de1_2 = min(ds1e1+ps1_2[3], ds1e1+pe1_2[3])
      }
      if(ps2_1[1]>pe2_1[1]) {
        de2_1 = min(ds2e2+ps2_1[3],ds2e2+pe2_1[3])
      }
    }
    return(max(ds1_2, de1_2,ds2_1, de2_1))
  } else if (type=="PPA"){ #Perpendicular/Parallel/Angle
    if(ds1e1 > ds2e2) { # Switch roles if longest segment is 1
      ds2e2 = dmat12[1,2]
      ds2s1 = dmat12[1,3]
      ds2e1 = dmat12[1,4]
      de2s1 = dmat12[2,3]
      de2e1 = dmat12[2,4]
      ds1e1 = dmat12[3,4]
    }
    #Assumes longer segment is 2
    ps1_2 =.distanceToSegment(ds2e2,ds1s2, ds1e2)
    pe1_2 =.distanceToSegment(ds2e2,de1s2, de1e2)
    lp1 = ps1_2[3]
    lp2 = pe1_2[3]
    dperpendicular = (lp1^2+lp2^2)/(lp1+lp2)
    lpar1 = min(ps1_2[1],ps1_2[2])
    lpar2 = min(pe1_2[1],pe1_2[2])
    dparallel = min(lpar1,lpar2)
    if(ps1_2[1]>pe1_2[1]) dangle = ds1e1
    else dangle = (max(lp2,lp1)-min(lp2,lp1))
    
    return(dperpendicular+dparallel+dangle)
  } 
  stop("Wrong distance type")
}
segmentDistances<-function(d, surveys, distance.type ="directed-segment") {
  dmat = as.matrix(d)
  n = nrow(dmat)
  nobj = n/surveys
  nseg = nobj*(surveys-1)
  dsegmat = matrix(0, nseg, nseg)
  rownames(dsegmat) =paste0(rep(1:nobj,surveys-1),"_",gl(surveys-1,nobj))
  colnames(dsegmat) =rownames(dsegmat)
  dinisegmat = dsegmat
  dfinsegmat = dsegmat
  dinifinsegmat = dsegmat
  dfininisegmat = dsegmat
  for(t1 in 1:(surveys-1)) {
    for(o1 in 1:nobj) {
      istart1 = (t1-1)*nobj + o1
      ifin1 = (t1)*nobj + o1
      for(t2 in 1:(surveys-1)) {
        for(o2 in 1:nobj) {
          istart2 = (t2-1)*nobj + o2
          ifin2 = (t2)*nobj + o2
          dmat12 = dmat[c(istart1,ifin1,istart2, ifin2),c(istart1,ifin1,istart2, ifin2)]
          # print(dmat12)
          dsegmat[istart1,istart2] <- .twoSegmentDistance(dmat12, type=distance.type)
          dsegmat[istart2,istart1] <- dsegmat[istart1,istart2]
          # print(c(istart1, istart2, dsegmat[istart1,istart2]))
          dinisegmat[istart2,istart1] <- dinisegmat[istart1,istart2]<-dmat[istart1,istart2]
          dfinsegmat[istart2,istart1] <- dfinsegmat[istart1,istart2] <- dmat[ifin1,ifin2]
          dinifinsegmat[istart1,istart2] <- dmat[istart1,ifin2]
          dinifinsegmat[istart2,istart1] <- dmat[istart2,ifin1]
          dfininisegmat[istart1,istart2] <- dmat[ifin1,istart2]
          dfininisegmat[istart2,istart1] <- dmat[ifin2,istart1]
          # print(c(istart1,ifin1,istart2, ifin2, dinifinsegmat[istart1,istart2],dfininisegmat[istart1,istart2]))
        }
      }
    }
  }
  return(list(Dseg = as.dist(dsegmat), Dini=as.dist(dinisegmat), Dfin = as.dist(dfinsegmat),
              Dinifin=dinifinsegmat, Dfinini = dfininisegmat))
}
trajectoryDistances<-function(d, surveys, distance.type="DSPD") {
  if(distance.type=="DSPD"){
    lsd = segmentDistances(d,surveys,distance.type="directed-segment")
  } else if(distance.type=="SPD" || distance.type=="Hausdorff") {
    lsd = segmentDistances(d,surveys,distance.type="Hausdorff")
  }
  dmat = as.matrix(lsd$Dseg)
  nseg = nrow(dmat)
  nobj = nseg/(surveys-1)
  dtraj = matrix(0, nrow=nobj, ncol = nobj)
  if(distance.type=="SPD" || distance.type=="DSPD") {
    for(i1 in 1:nobj) {
      for(i2 in i1:nobj) {
        dt12 = 0
        for(t1 in 1:(surveys-1)) {
          dt12ivec = numeric(0)
          iseg1 = i1+(t1-1)*nobj
          for(t2 in 1:(surveys-1)) {
            iseg2 = i2+(t2-1)*nobj
            dt12ivec = c(dt12ivec, dmat[iseg1, iseg2])
          }
          dt12 = dt12 + min(dt12ivec)
        }
        dt21 = 0 
        for(t2 in 1:(surveys-1)) {
          dt21ivec = numeric(0)
          iseg2 = i2+(t2-1)*nobj
          for(t1 in 1:(surveys-1)) {
            iseg1 = i1+(t1-1)*nobj
            dt21ivec = c(dt21ivec, dmat[iseg1, iseg2])
          }
          dt21 = dt21 + min(dt21ivec)
        }
        dtraj[i1,i2] = (dt12+dt21)/2
        dtraj[i2,i1] = dtraj[i1,i2]
      }
    }
  } else if(distance.type=="Hausdorff") {
    for(i1 in 1:nobj) {
      for(i2 in i1:nobj) {
        dt12ivec = numeric(0)
        for(t1 in 1:(surveys-1)) {
          iseg1 = i1+(t1-1)*nobj
          for(t2 in 1:(surveys-1)) {
            iseg2 = i2+(t2-1)*nobj
            dt12ivec = c(dt12ivec, dmat[iseg1, iseg2])
          }
        }
        dtraj[i1,i2] = max(dt12ivec)
        dtraj[i2,i1] = dtraj[i1,i2]
      }
    }
    
  }
  return(as.dist(dtraj))
}
