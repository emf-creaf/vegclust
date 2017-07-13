.projection<-function(dref, d1, d2) {
  a1 = (d1^2+dref^2-d2^2)/(2*dref)
  a2 = dref-a1
  s = d1^2-a1^2
  h = ifelse(s>=0, sqrt(s), NA)
  return(c(a1,a2,h))
}
.distanceToSegment<-function(dref, d1,d2) {
  p =.projection(dref,d1, d2)
  if(is.na(p[3]) || p[1]<0 || p[2]<0) {
    if(d1<d2) {
      p[1] = 0
      p[2] = dref
      p[3] = d1
    } else {
      p[1] = dref
      p[2] = 0
      p[3] = d2
    }
  }
  return(p)
}
.triangleinequality<-function(d1,d2,d3, tol=0.0001) {
  if((d1+d2)<d3*(1-tol)) return(FALSE)
  else if((d1+d3)<d2*(1-tol)) return(FALSE)
  else if((d2+d3)<d1*(1-tol)) return(FALSE)
  return(TRUE)
}
.pt<-function(dIT, dXT, dPX, dPI) {
  if(dPI==0) return(dIT)
  else if(dPX==0) return(dXT)
  A = (dXT^2)-(dPX^2)
  B = (dIT^2)-(dPI^2)
  
  ax =(dXT^2)+(dIT^2)
  bx = ((dXT^2)*2*B) + ((dIT^2)*2*A) - (4*(dIT^2)*(dXT^2))
  cx = (dXT^2)*(B^2) + (dIT^2)*(A^2)
  z = (bx^2)-(4*ax*cx)
  if(z<0) {
    cat(dIT,dXT,dPX, dPI," - ",z,"\n")
    stop()
  }
  d2 = ((-1)*bx + sqrt(z))/(2*ax)
  return(sqrt(d2))
}

is.metric<-function(D, tol=0.0001) {
  dmat = as.matrix(D)
  n = nrow(dmat)
  for(i in 1:n) {
    for(j in i:n) {
      for(k in j:n) {
        ti = .triangleinequality(dmat[i,j], dmat[i,k], dmat[j,k],tol)
        if(!ti) {
          return(FALSE)
        }
      }
    }
  }
  return(TRUE)
}

#Calculates Hausdorff distance between two line segments
distanceToTrajectory<-function(dsteps, d2ref, eps) {
  
  nsteps = length(dsteps)
  npoints = nrow(d2ref)
  
  #Cumulative distance between steps
  dstepcum = rep(0,nsteps+1)
  for(i in 2:nsteps) {
    dstepcum[i] = dstepcum[i-1]+dsteps[i-1]
  }
  dstepcum[nsteps+1] = sum(dsteps)
  
  projH = matrix(NA, nrow=npoints, ncol = nsteps)
  projA1 = matrix(NA, nrow=npoints, ncol = nsteps)
  projA2 = matrix(NA, nrow=npoints, ncol = nsteps)
  projIn = matrix(FALSE, nrow=npoints, ncol = nsteps)
  whichstep = rep(NA, npoints)
  dgrad = rep(NA, npoints)
  posgrad = rep(NA, npoints)
  
  for(i in 1:npoints) {
    for(j in 1:nsteps) {
      if(!.triangleinequality(dsteps[j], d2ref[i, j], d2ref[i, j+1])) cat(paste0(i," to [",j,", ",j+1,"] does not meet triangle inequality\n"))
      p <-.projection(dsteps[j], d2ref[i, j], d2ref[i, j+1])
      projA1[i,j] = p[1]
      projA2[i,j] = p[2]
      projH[i,j] = p[3]
      if(!is.na(projH[i,j])) {
        projIn[i,j] = (projH[i,j]< eps) && (projA1[i,j]>0) && (projA2[i,j]>0)
        if(is.na(projIn[i,j])) projIn[i,j] = FALSE
      } else {
        projIn[i,j] = FALSE
      }
    }
    if(sum(projIn[i,])>0) {
      h = projH[i,projIn[i,]]
      whichstep[i] = which(projIn[i,])[which.min(h)]
      dg = dstepcum[whichstep[i]]+projA1[i,whichstep[i]]
      posgrad[i] = dg/sum(dsteps)
      dgrad[i] = min(h)
    } else { #Check if point lies in gradient turns
      jmin = which.min(d2ref[i,])
      dgrad[i] = d2ref[i,jmin]
      if(d2ref[i,jmin]<eps) {
        # cat("in angle\n")
        posgrad[i] = dstepcum[jmin]/sum(dsteps)
      }
    }
  }
  res = cbind(dgrad, posgrad)
  row.names(res)<-row.names(d2ref)
  return(res)
}
