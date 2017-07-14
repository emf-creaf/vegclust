
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

segmentDistances<-function(d, sites, surveys=NULL, distance.type ="directed-segment") {
  distance.type <- match.arg(distance.type, c("directed-segment", "Hausdorff", "PPA"))
  if(length(sites)!=nrow(as.matrix(d))) stop("'sites' needs to be of length equal to the number of rows/columns in d")
  if(!is.null(surveys)) if(length(sites)!=length(surveys)) stop("'sites' and 'surveys' need to be of the same length")
  
  siteIDs = unique(sites)
  nsite = length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) {
    nsurveysite[i] = sum(sites==siteIDs[i])
  }
  if(sum(nsurveysite==1)>0) stop("All sites need to be surveyed at least twice")
  dmat = as.matrix(d)
  n = nrow(dmat)
  nseg = sum(nsurveysite)-nsite
  segnames = character(nseg)
  cnt=1
  for(i in 1:nsite) {
    if(!is.null(surveys)) surv = surveys[sites==siteIDs[i]]
    else surv = 1:nsurveysite[i]
    for(j in 1:(nsurveysite[i]-1)) {
      segnames[cnt] = paste0(siteIDs[i],"[",surv[j],"-",surv[j+1],"]")
      cnt = cnt+1
    }
  }
  dsegmat = matrix(0, nseg, nseg)
  rownames(dsegmat) =segnames
  colnames(dsegmat) =segnames
  dinisegmat = dsegmat
  dfinsegmat = dsegmat
  dinifinsegmat = dsegmat

  os1 = 1
  for(i1 in 1:nsite) {
    ind_surv1 = which(sites==siteIDs[i1])
    for(s1 in 1:(nsurveysite[i1]-1)) {
      os2 = 1
      for(i2 in 1:nsite) {
        ind_surv2 = which(sites==siteIDs[i2])
        for(s2 in 1:(nsurveysite[i2]-1)) {
          # os2 = sum(nsurveysite[1:(i2-1)]-1)+s2 #Output index of segment 2
          #Select submatrix from dmat
          ind12 = c(ind_surv1[s1],ind_surv1[s1+1],ind_surv2[s2],ind_surv2[s2+1])
          # print(ind12)
          # print(c(os1, os2))
          dmat12 = dmat[c(ind_surv1[s1],ind_surv1[s1+1],ind_surv2[s2],ind_surv2[s2+1]),
                        c(ind_surv1[s1],ind_surv1[s1+1],ind_surv2[s2],ind_surv2[s2+1])]
          dsegmat[os1,os2] <- .twoSegmentDistance(dmat12, type=distance.type)
          dsegmat[os2,os1] <- dsegmat[os1,os2]
          dinisegmat[os2,os1] <- dinisegmat[os1,os2]<-dmat[ind_surv1[s1],ind_surv2[s2]]
          dfinsegmat[os2,os1] <- dfinsegmat[os1,os2]<-dmat[ind_surv1[s1+1],ind_surv2[s2+1]]
          dinifinsegmat[os1,os2] <- dmat[ind_surv1[s1],ind_surv2[s2+1]]
          dinifinsegmat[os2,os1] <- dmat[ind_surv1[s1+1],ind_surv2[s2]]
          os2 = os2+1
        }
      }
      os1 = os1+1
    }
  }

  return(list(Dseg = as.dist(dsegmat), Dini=as.dist(dinisegmat), Dfin = as.dist(dfinsegmat),
              Dinifin=dinifinsegmat))
}

trajectoryDistances<-function(d, sites, surveys=NULL, distance.type="DSPD") {
  distance.type <- match.arg(distance.type, c("DSPD", "SPD", "Hausdorff"))
  if(length(sites)!=nrow(as.matrix(d))) stop("'sites' needs to be of length equal to the number of rows/columns in d")
  if(!is.null(surveys)) if(length(sites)!=length(surveys)) stop("'sites' and 'surveys' need to be of the same length")
  siteIDs = unique(sites)
  nsite = length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) nsurveysite[i] = sum(sites==siteIDs[i])
  if(sum(nsurveysite==1)>0) stop("All sites need to be surveyed at least twice")
  n = nrow(as.matrix(d))
  nseg = sum(nsurveysite)-nsite
  
  #Init output
  dtraj = matrix(0, nrow=nsite, ncol = nsite)
  rownames(dtraj) = siteIDs
  colnames(dtraj) = siteIDs
  if(distance.type=="DSPD"){
    lsd = segmentDistances(d,sites, surveys,distance.type="directed-segment")
    dsegmat = as.matrix(lsd$Dseg)
    for(i1 in 1:nsite) {
      for(i2 in 1:nsite) {
        dt12 = 0
        for(s1 in 1:(nsurveysite[i1]-1)) {
          dt12ivec = numeric(0)
          iseg1 = sum(nsurveysite[1:i1]-1)-(nsurveysite[i1]-1)+s1
          for(s2 in 1:(nsurveysite[i2]-1)) {
            iseg2 = sum(nsurveysite[1:i2]-1)-(nsurveysite[i2]-1)+s2
            dt12ivec = c(dt12ivec, dsegmat[iseg1, iseg2])
          }
          dt12 = dt12 + min(dt12ivec)
        }
        dt12 = dt12/(nsurveysite[i1]-1) #Average of distances between segments of T1 and trajectory T2
        dt21 = 0 
        for(s2 in 1:(nsurveysite[i2]-1)) {
          dt21ivec = numeric(0)
          iseg2 = sum(nsurveysite[1:i2]-1)-(nsurveysite[i2]-1)+s2
          for(s1 in 1:(nsurveysite[i1]-1)) {
            iseg1 = sum(nsurveysite[1:i1]-1)-(nsurveysite[i1]-1)+s1
            dt21ivec = c(dt21ivec, dsegmat[iseg1, iseg2])
          }
          dt21 = dt21 + min(dt21ivec)
        }
        dt21 = dt21/(nsurveysite[i2]-1) #Average of distances between segments of T2 and trajectory T1
        
        dtraj[i1,i2] = (dt12+dt21)/2 #Symmetrization
        dtraj[i2,i1] = dtraj[i1,i2]
      }
    }
    
  } 
  else if(distance.type=="SPD") {
    dmat = as.matrix(d)
    for(i1 in 1:nsite) {
      ind_surv1 = which(sites==siteIDs[i1])
      for(i2 in 1:nsite) {
        ind_surv2 = which(sites==siteIDs[i2])
        dt12 = 0
        for(p1 in 1:nsurveysite[i1]) {
          dt12ivec = numeric(0)
          ip1 = ind_surv1[p1]
          for(s2 in 1:(nsurveysite[i2]-1)) {
            ipi2 = ind_surv2[s2] #initial point
            ipe2 = ind_surv2[s2+1] #end point
            dt12ivec = c(dt12ivec, .distanceToSegment(dmat[ipi2,ipe2], dmat[ip1, ipi2], dmat[ip1,ipe2])[3])
          }
          dt12 = dt12 + min(dt12ivec)
        }
        dt12 = dt12/nsurveysite[i1] #Average of distances between points of T1 and trajectory T2
        dt21 = 0 
        for(p2 in 1:nsurveysite[i2]) {
          dt21ivec = numeric(0)
          ip2 = ind_surv2[p2]
          for(s1 in 1:(nsurveysite[i1]-1)) {
            ipi1 = ind_surv1[s1] #initial point
            ipe1 = ind_surv1[s1+1] #end point
            dt21ivec = c(dt21ivec, .distanceToSegment(dmat[ipi1,ipe1], dmat[ip2, ipi1], dmat[ip2,ipe1])[3])
          }
          dt21 = dt21 + min(dt21ivec)
        }
        dt21 = dt21/nsurveysite[i2] #Average of distances between points of T2 and trajectory T1
        
        dtraj[i1,i2] = (dt12+dt21)/2 #Symmetrization
        dtraj[i2,i1] = dtraj[i1,i2]
      }
    }
  }
  else if(distance.type=="Hausdorff") {
    dmat = as.matrix(d)
    for(i1 in 1:nsite) {
      ind_surv1 = which(sites==siteIDs[i1])
      for(i2 in 1:nsite) {
        ind_surv2 = which(sites==siteIDs[i2])
        dt12 = 0
        dt12vec = numeric(0)
        for(p1 in 1:nsurveysite[i1]) {
          ip1 = ind_surv1[p1]
          for(s2 in 1:(nsurveysite[i2]-1)) {
            ipi2 = ind_surv2[s2] #initial point
            ipe2 = ind_surv2[s2+1] #end point
            dt12vec = c(dt12vec, .distanceToSegment(dmat[ipi2,ipe2], dmat[ip1, ipi2], dmat[ip1,ipe2])[3])
          }
        }
        dt12 = max(dt12vec) #Maximum of distances between points of T1 and segments of T2
        dt21 = 0 
        dt21vec = numeric(0)
        for(p2 in 1:nsurveysite[i2]) {
          ip2 = ind_surv2[p2]
          for(s1 in 1:(nsurveysite[i1]-1)) {
            ipi1 = ind_surv1[s1] #initial point
            ipe1 = ind_surv1[s1+1] #end point
            dt21vec = c(dt21vec, .distanceToSegment(dmat[ipi1,ipe1], dmat[ip2, ipi1], dmat[ip2,ipe1])[3])
          }
        }
        dt21 = max(dt21vec) #Maximum of distances between points of T2 and segments of T1
        
        dtraj[i1,i2] = max(dt12, dt21) #maximum of maximums
        dtraj[i2,i1] = dtraj[i1,i2]
      }
    }
  } 
  else stop("Wrong distance type")
  return(as.dist(dtraj))
}
