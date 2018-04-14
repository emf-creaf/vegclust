#' Community trajectory analysis
#' 
#' Given a distance matrix between community states, functions \code{segmentDistances} and \code{trajectoryDistances} calculate the distance between pairs of directed segments and community trajectories, respectively. 
#' Function \code{trajectoryLengths} calculates lengths of directed segments and complete trajectories. 
#' Function \code{trajectoryAngles} calculates the angle between consecutive pairs of directed segments.
#' Function \code{trajectoryPCoA} performs principal coordinates analysis (\code{\link{cmdscale}}) and draws trajectories in the ordination scatterplot.
#' Function \code{trajectoryProjection} projects a set of target points onto a specified trajectory and returns the distance to the trajectory (i.e. rejection) and the relative position of the projection point within the trajectory.
#' 
#' These functions consider community dynamics as trajectories in a chosen space of community resemblance and takes trajectories as objects to be compared. 
#' By adapting concepts and procedures used for the analysis of trajectories in space (i.e. movement data) (Besse et al. 2016), the functions allow assessing the resemblance between trajectories. 
#' Details of calculations are given in De \enc{Cáceres}{Caceres} et al (submitted)
#' 
#' @encoding UTF-8
#' @name trajectories
#' @aliases segmentDistances trajectoryDistances trajectoryLengths trajectoryAngles trajectoryPCoA trajectoryProjection trajectoryConvergence trajectoryDirectionality
#' 
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of community states.
#' @param sites A vector indicating the site corresponding to each community state.
#' @param surveys A vector indicating the survey corresponding to each community state (only necessary when surveys are not in order).
#' @param distance.type 
#' The type of distance index to be calculated (Besse et al. 2016; De Cáceres et al. submitted). For \code{segmentDistances} the available indices are:
#'   \itemize{
#'     \item{\code{Hausdorff}: Hausdorff distance between two segments.}
#'     \item{\code{directed-segment}: Directed segment distance (default).}
#'     \item{\code{PPA}: Perpendicular-parallel-angle distance.}
#'   }
#' whereas for \code{trajectoryDistances} the available indices are:
#'   \itemize{
#'     \item{\code{Hausdorff}: Hausdorff distance between two trajectories.}
#'     \item{\code{SPD}: Segment path distance.}
#'     \item{\code{DSPD}: Directed segment path distance (default).}
#'   }
#' @param symmetrization Function used to obtain a symmetric distance, so that DSPD(T1,T2) = DSPD(T2,T1) (e.g., \code{mean} or \code{min}).
#' @param verbose Provides console output informing about process (useful for large dataset).
#' 
#' @return Function \code{trajectoryDistances} returns an object of class \code{\link{dist}} containing the distances between trajectories. Function \code{trajectorySegments} returns a list with the following elements:
#' \itemize{
#'   \item{\code{Dseg}: Distance matrix between segments.}
#'   \item{\code{Dini}: Distance matrix between initial points of segments.}
#'   \item{\code{Dfin}: Distance matrix between final points of segments.}
#'   \item{\code{Dinifin}: Distance matrix between initial points of one segment and the final point of the other.}
#'   \item{\code{Dfinini}: Distance matrix between final points of one segment and the initial point of the other.}
#' }
#' 
#' Function \code{trajectoryLengths} returns a data frame with the length of each segment on each trajectory and the total length of all trajectories. Function \code{trajectoryPCoA} returns the result of calling \code{\link{cmdscale}}.
#' 
#' Function \code{trajectoryAngles} returns a data frame with the angle between each pair of segments on each trajectory and the mean and standard deviation of those angles across each trajectory. 
#' 
#' Function \code{trajectoryPCoA} returns the result of calling \code{\link{cmdscale}}.
#' 
#' Function \code{trajectoryProjection} returns a data frame with the following columns:
#' \itemize{
#'   \item{\code{distanceToTrajectory}: Distances to the trajectory, i.e. rejection (\code{NA} for target points whose projection is outside the trajectory).}
#'   \item{\code{segment}: Segment that includes the projected point (\code{NA} for target points whose projection is outside the trajectory).}
#'   \item{\code{relativePosition}: Relative position of the projected point within the trajectory, i.e. values from 0 to 1 with 0 representing the start of the trajectory and 1 representing the end (\code{NA} for target points whose projection is outside the trajectory).}
#' }
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, Forest Sciences Center of Catalonia
#' 
#' @references
#' Besse, P., Guillouet, B., Loubes, J.-M. & François, R. (2016). Review and perspective for distance based trajectory clustering. IEEE Trans. Intell. Transp. Syst., 17, 3306–3317.
#' 
#' De \enc{Cáceres}{Caceres} M, Coll L, Legendre P, Allen RB, Wiser SK, Fortin MJ, Condit R & Hubbell S. (submitted). Trajectory analysis in community ecology.
#' 
#' @seealso \code{\link{cmdscale}}
#' 
#' @examples 
#'   #Description of sites and surveys
#'   sites = c(1,1,1,2,2,2)
#'   surveys=c(1,2,3,1,2,3)
#'   
#'   #Raw data table
#'   xy<-matrix(0, nrow=6, ncol=2)
#'   xy[2,2]<-1
#'   xy[3,2]<-2
#'   xy[4:6,1] <- 0.5
#'   xy[4:6,2] <- xy[1:3,2]
#'   xy[6,1]<-1
#'   
#'   #Distance matrix
#'   d = dist(xy)
#'   d
#'   
#'   trajectoryLengths(d, sites, surveys)
#'   trajectoryAngles(d, sites, surveys)
#'   segmentDistances(d, sites, surveys)$Dseg
#'   trajectoryDistances(d, sites, surveys, distance.type = "Hausdorff")
#'   trajectoryDistances(d, sites, surveys, distance.type = "DSPD")
#'   
#'   #Draw trajectories
#'   trajectoryPCoA(d, sites, traj.colors = c("black","red"), lwd = 2)
#'   
#'   
#'   #Should give the same results if surveys are not in order 
#'   #(here we switch surveys for site 2)
#'   temp = xy[5,]
#'   xy[5,] = xy[6,]
#'   xy[6,] = temp
#'   surveys[5] = 3
#'   surveys[6] = 2
#'   trajectoryLengths(dist(xy), sites, surveys)
#'   segmentDistances(dist(xy), sites, surveys)$Dseg
#'   trajectoryDistances(dist(xy), sites, surveys, distance.type = "Hausdorff")
#'   trajectoryDistances(dist(xy), sites, surveys, distance.type = "DSPD")
#'  
segmentDistances<-function(d, sites, surveys=NULL, distance.type ="directed-segment", verbose=FALSE) {
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
    if(!is.null(surveys)) {
      surv = surveys[sites==siteIDs[i]]
      surv = sort(surv) #Surveys may not be in order
    }
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
  if(verbose) {
    cat("\nCalculating segment distances...\n")
    tb = txtProgressBar(1, nsite, style=3)
  }
  for(i1 in 1:nsite) {
    if(verbose) setTxtProgressBar(tb, i1)
    ind_surv1 = which(sites==siteIDs[i1])
    #Surveys may not be in order
    if(!is.null(surveys)) ind_surv1 = ind_surv1[order(surveys[sites==siteIDs[i1]])]
    for(s1 in 1:(nsurveysite[i1]-1)) {
      os2 = 1
      for(i2 in 1:nsite) {
        ind_surv2 = which(sites==siteIDs[i2])
        #Surveys may not be in order
        if(!is.null(surveys)) ind_surv2 = ind_surv2[order(surveys[sites==siteIDs[i2]])]
        for(s2 in 1:(nsurveysite[i2]-1)) {
          # os2 = sum(nsurveysite[1:(i2-1)]-1)+s2 #Output index of segment 2
          #Select submatrix from dmat
          ind12 = c(ind_surv1[s1],ind_surv1[s1+1],ind_surv2[s2],ind_surv2[s2+1])
          # print(ind12)
          # print(c(os1, os2))
          dmat12 = dmat[c(ind_surv1[s1],ind_surv1[s1+1],ind_surv2[s2],ind_surv2[s2+1]),
                        c(ind_surv1[s1],ind_surv1[s1+1],ind_surv2[s2],ind_surv2[s2+1])]
          dsegmat[os1,os2] <- .twoSegmentDistanceC(dmat12, type=distance.type)
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

#' @rdname trajectories
trajectoryDistances<-function(d, sites, surveys=NULL, distance.type="DSPD", symmetrization = "mean" , verbose=FALSE) {
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
    lsd = segmentDistances(d,sites, surveys,distance.type="directed-segment", verbose)
    dsegmat = as.matrix(lsd$Dseg)
    if(verbose) {
      cat("\nCalculating trajectory distances...\n")
      tb = txtProgressBar(1, nsite, style=3)
    }
    for(i1 in 1:nsite) {
      if(verbose) setTxtProgressBar(tb, i1)
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
        
        dtraj[i1,i2] = do.call(symmetrization, list(c(dt12,dt21))) #Symmetrization
        dtraj[i2,i1] = dtraj[i1,i2]
      }
    }
    
  } 
  else if(distance.type=="SPD") {
    dmat = as.matrix(d)
    for(i1 in 1:nsite) {
      ind_surv1 = which(sites==siteIDs[i1])
      #Surveys may not be in order
      if(!is.null(surveys)) ind_surv1 = ind_surv1[order(surveys[sites==siteIDs[i1]])]
      for(i2 in 1:nsite) {
        ind_surv2 = which(sites==siteIDs[i2])
        #Surveys may not be in order
        if(!is.null(surveys)) ind_surv2 = ind_surv2[order(surveys[sites==siteIDs[i2]])]
        dt12 = 0
        for(p1 in 1:nsurveysite[i1]) {
          dt12ivec = numeric(0)
          ip1 = ind_surv1[p1]
          for(s2 in 1:(nsurveysite[i2]-1)) {
            ipi2 = ind_surv2[s2] #initial point
            ipe2 = ind_surv2[s2+1] #end point
            dt12ivec = c(dt12ivec, .distanceToSegmentC(dmat[ipi2,ipe2], dmat[ip1, ipi2], dmat[ip1,ipe2])[3])
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
            dt21ivec = c(dt21ivec, .distanceToSegmentC(dmat[ipi1,ipe1], dmat[ip2, ipi1], dmat[ip2,ipe1])[3])
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
      #Surveys may not be in order
      if(!is.null(surveys)) ind_surv1 = ind_surv1[order(surveys[sites==siteIDs[i1]])]
      for(i2 in 1:nsite) {
        ind_surv2 = which(sites==siteIDs[i2])
        #Surveys may not be in order
        if(!is.null(surveys)) ind_surv2 = ind_surv2[order(surveys[sites==siteIDs[i2]])]
        dt12 = 0
        dt12vec = numeric(0)
        for(p1 in 1:nsurveysite[i1]) {
          ip1 = ind_surv1[p1]
          for(s2 in 1:(nsurveysite[i2]-1)) {
            ipi2 = ind_surv2[s2] #initial point
            ipe2 = ind_surv2[s2+1] #end point
            dt12vec = c(dt12vec, .distanceToSegmentC(dmat[ipi2,ipe2], dmat[ip1, ipi2], dmat[ip1,ipe2])[3])
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
            dt21vec = c(dt21vec, .distanceToSegmentC(dmat[ipi1,ipe1], dmat[ip2, ipi1], dmat[ip2,ipe1])[3])
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

#' @rdname trajectories
trajectoryLengths<-function(d, sites, surveys=NULL, verbose= FALSE) {
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

  maxnsurveys = max(nsurveysite)

  lengths = as.data.frame(matrix(NA, nrow=nsite, ncol=maxnsurveys))
  row.names(lengths)<-siteIDs
  names(lengths)<-c(paste0("S",as.character(1:(maxnsurveys-1))),"Trajectory")
  if(verbose) {
    cat("\nCalculating trajectory lengths...\n")
    tb = txtProgressBar(1, nsite, style=3)
  }
  for(i1 in 1:nsite) {
    if(verbose) setTxtProgressBar(tb, i1)
    ind_surv1 = which(sites==siteIDs[i1])
    #Surveys may not be in order
    if(!is.null(surveys)) ind_surv1 = ind_surv1[order(surveys[sites==siteIDs[i1]])]
    for(s1 in 1:(nsurveysite[i1]-1)) {
      lengths[i1,s1] = dmat[ind_surv1[s1], ind_surv1[s1+1]]
    }
    lengths[i1, maxnsurveys] = sum(lengths[i1,], na.rm=T)
  }
  return(lengths)
}

#' @rdname trajectories
trajectoryAngles<-function(d, sites, surveys=NULL, verbose= FALSE) {
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
  
  maxnsurveys = max(nsurveysite)
  angles = matrix(NA, nrow=nsite, ncol=maxnsurveys)
  if(verbose) {
    cat("\nCalculating trajectory angles...\n")
    tb = txtProgressBar(1, nsite, style=3)
  }
  for(i1 in 1:nsite) {
    if(verbose) setTxtProgressBar(tb, i1)
    ind_surv1 = which(sites==siteIDs[i1])
    # print(ind_surv1)
    #Surveys may not be in order
    if(!is.null(surveys)) ind_surv1 = ind_surv1[order(surveys[sites==siteIDs[i1]])]
    for(s1 in 1:(nsurveysite[i1]-2)) {
      d12 = dmat[ind_surv1[s1], ind_surv1[s1+1]]
      d23 = dmat[ind_surv1[s1+1], ind_surv1[s1+2]]
      d13 = dmat[ind_surv1[s1], ind_surv1[s1+2]]
      angles[i1, s1] = .angleConsecutiveC(d12,d23,d13, TRUE)
      # cat(paste(i1,s1,":", d12,d23,d13,.angleConsecutiveC(d12,d23,d13, TRUE),"\n"))
    }
    angles[i1, maxnsurveys-1] = mean(angles[i1,1:(nsurveysite[i1]-2)], na.rm=T)
    angles[i1, maxnsurveys] = sd(angles[i1,1:(nsurveysite[i1]-2)], na.rm=T)
  }
  angles = as.data.frame(angles)
  row.names(angles)<-siteIDs
  names(angles)<-c(paste0("S",as.character(1:(maxnsurveys-2)),"-S",as.character(2:(maxnsurveys-1))),"mean", "sd")
  return(angles)
}

#' @rdname trajectories
#' @param selection A numeric or logical vector of the same length as \code{sites}, indicating a subset of site trajectories to be plotted.
#' @param traj.colors A vector of colors (one per site). If \code{selection != NULL} the length of the color vector should be equal to the number of sites selected.
#' @param axes The pair of principal coordinates to be plotted.
#' @param ... Additional parameters for function \code{\link{arrows}}.
trajectoryPCoA<-function(d, sites, surveys = NULL, selection = NULL, traj.colors = NULL, axes=c(1,2), ...) {
  siteIDs = unique(sites)
  nsite = length(siteIDs)
  
  #Apply site selection
  
  if(is.null(selection)) selection = 1:nsite 
  selIDs = siteIDs[selection]
  
  D2 =as.dist(as.matrix(d)[sites %in% selIDs, sites %in% selIDs])
  cmd_D2<-cmdscale(D2,eig=TRUE, add=TRUE, k=nrow(as.matrix(D2))-1)
  
  x<-cmd_D2$points[,axes[1]]
  y<-cmd_D2$points[,axes[2]]
  plot(x,y, type="n", asp=1, xlab=paste0("PCoA ",axes[1]," (", round(100*cmd_D2$eig[axes[1]]/sum(cmd_D2$eig)),"%)"), 
       ylab=paste0("PCoA ",axes[2]," (", round(100*cmd_D2$eig[axes[2]]/sum(cmd_D2$eig)),"%)"))
  
  sitesred = sites[sites %in% selIDs]
  if(!is.null(surveys)) surveysred = surveys[sites %in% selIDs]
  else surveysred = NULL
  #Draw arrows
  for(i in 1:length(selIDs)) {
    ind_surv = which(sitesred==selIDs[i])
    #Surveys may not be in order
    if(!is.null(surveysred)) ind_surv = ind_surv[order(surveysred[sitesred==selIDs[i]])]
    for(t in 1:(length(ind_surv)-1)) {
      niini =ind_surv[t]
      nifin =ind_surv[t+1]
      if(!is.null(traj.colors)) arrows(x[niini],y[niini],x[nifin],y[nifin], col = traj.colors[i], ...)
      else arrows(x[niini],y[niini],x[nifin],y[nifin], ...)
    }
  }
  #Draw legend
  invisible(cmd_D2)
}


#' @rdname trajectories
#' @param target An integer vector of the community states to be projected.
#' @param trajectory An integer vector of the trajectory onto which target states are to be projected.
#' @param tol Numerical tolerance value to determine that projection of a point lies within the trajectory.
trajectoryProjection<-function(d, target, trajectory, tol = 0.000001) {
  if(length(trajectory)<2) stop("Trajectory needs to include at least two states")
  dmat = as.matrix(d)
  npoints = length(target)
  nsteps = length(trajectory) -1
  #Distance betwen target points and trajectory points
  d2ref = dmat[target, trajectory, drop=FALSE]
  #Distance between trajectory steps
  dsteps = diag(dmat[trajectory[1:(length(trajectory)-1)], trajectory[2:length(trajectory)]])
  #Cumulative distance between steps
  dstepcum = rep(0,nsteps+1)
  if(nsteps>1) {
    for(i in 2:nsteps) {
      dstepcum[i] = dstepcum[i-1]+dsteps[i-1]
    }
  }
  dstepcum[nsteps+1] = sum(dsteps)
  
  projH = matrix(NA, nrow=npoints, ncol = nsteps)
  projA1 = matrix(NA, nrow=npoints, ncol = nsteps)
  projA2 = matrix(NA, nrow=npoints, ncol = nsteps)
  whichstep = rep(NA, npoints)
  dgrad = rep(NA, npoints)
  posgrad = rep(NA, npoints)
  
  for(i in 1:npoints) {
    for(j in 1:nsteps) {
      if(!.triangleinequalityC(dsteps[j], d2ref[i, j], d2ref[i, j+1])) warning(paste0(i," to [",j,", ",j+1,"] does not meet triangle inequality\n"))
      p <-.projectionC(dsteps[j], d2ref[i, j], d2ref[i, j+1])
      if((!is.na(p[3])) & (p[1]>-tol) & (p[2]>-tol)) {
        projA1[i,j] = p[1]
        projA2[i,j] = p[2]
        projH[i,j] = p[3]
        if(is.na(dgrad[i])) {
          dgrad[i] = p[3]
          whichstep[i] = j
        } else {
          if(p[3]<dgrad[i]) {
            dgrad[i] = p[3]
            whichstep[i] = j
          }
        }
      }
    }
    if(!is.na(whichstep[i])) {
      dg = dstepcum[whichstep[i]]+projA1[i,whichstep[i]]
      posgrad[i] = dg/sum(dsteps)
    }
  }
  res = data.frame(distanceToTrajectory=dgrad, segment = whichstep, relativePosition = posgrad)
  row.names(res)<-row.names(d2ref)
  return(res)
}


#' @rdname trajectories
#' @param symmetric A logical flag to indicate a symmetric convergence comparison of trajectories.
trajectoryConvergence<-function(d, sites, surveys, symmetric = FALSE, verbose = FALSE){
  if(length(sites)!=nrow(as.matrix(d))) stop("'sites' needs to be of length equal to the number of rows/columns in d")
  if(!is.null(surveys)) if(length(sites)!=length(surveys)) stop("'sites' and 'surveys' need to be of the same length")
  siteIDs = unique(sites)
  nsite = length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) nsurveysite[i] = sum(sites==siteIDs[i])
  if(sum(nsurveysite<3)>0) stop("All sites need to be surveyed at least three times")
  n = nrow(as.matrix(d))

  #Init output
  tau = matrix(NA, nrow=nsite, ncol = nsite)
  rownames(tau) = siteIDs
  colnames(tau) = siteIDs
  p.value = tau
  dmat = as.matrix(d)
  if(verbose) {
    cat("\nCalculating trajectory convergence...\n")
    tb = txtProgressBar(1, nsite, style=3)
  }
  for(i1 in 1:(nsite-1)) {
    if(verbose) setTxtProgressBar(tb, i1)
    ind_surv1 = which(sites==siteIDs[i1])
    #Surveys may not be in order
    if(!is.null(surveys)) ind_surv1 = ind_surv1[order(surveys[sites==siteIDs[i1]])]
    for(i2 in (i1+1):nsite) {
      ind_surv2 = which(sites==siteIDs[i2])
      #Surveys may not be in order
      if(!is.null(surveys)) ind_surv2 = ind_surv2[order(surveys[sites==siteIDs[i2]])]
      if(!symmetric) {
        trajectory = ind_surv2
        target = ind_surv1
        trajProj = trajectoryProjection(d,target, trajectory)
        dT = trajProj$distanceToTrajectory
        mk.test = MannKendall(dT)
        tau[i1,i2] = mk.test$tau
        p.value[i1,i2] = mk.test$sl
        trajectory = ind_surv1
        target = ind_surv2
        trajProj = trajectoryProjection(d,target, trajectory)
        dT = trajProj$distanceToTrajectory
        mk.test = MannKendall(dT)
        tau[i2,i1] = mk.test$tau
        p.value[i2,i1] = mk.test$sl
      } 
      else {
        if(length(ind_surv1)==length(ind_surv2)) {
          dT = numeric(length(ind_surv1))
          for(j in 1:length(ind_surv1)) dT[j] = dmat[ind_surv1[j], ind_surv2[j]]
          mk.test = MannKendall(dT)
          tau[i1,i2] = mk.test$tau
          p.value[i1,i2] = mk.test$sl
          tau[i2,i1] = mk.test$tau
          p.value[i2,i1] = mk.test$sl
        } else {
          warning(paste0("sites ",i1, " and ",i2," do not have the same number of surveys."))
        }
      }
    }
  }
  return(list(tau = tau, p.value = p.value))
}


#' @rdname trajectories
#' @param nperm Number of permutations in the Mantel test
trajectoryDirectionality<-function(d, sites, surveys, nperm = 999, verbose = FALSE) {
  if(length(sites)!=nrow(as.matrix(d))) stop("'sites' needs to be of length equal to the number of rows/columns in d")
  if(!is.null(surveys)) if(length(sites)!=length(surveys)) stop("'sites' and 'surveys' need to be of the same length")
  siteIDs = unique(sites)
  nsite = length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) nsurveysite[i] = sum(sites==siteIDs[i])
  if(sum(nsurveysite<3)>0) stop("All sites need to be surveyed at least three times")
  n = nrow(as.matrix(d))
  
  dmat = as.matrix(d)
  #Init output
  res = data.frame(r= rep(NA,nsite), p.value = rep(NA,nsite), row.names = siteIDs)
  if(verbose) {
    cat("\nCalculating trajectory directionality...\n")
    tb = txtProgressBar(1, nsite, style=3)
  }
  for(i1 in 1:nsite) {
    if(verbose) setTxtProgressBar(tb, i1)
    ind_surv1 = which(sites==siteIDs[i1])
    #Surveys may not be in order
    if(!is.null(surveys)) ind_surv1 = ind_surv1[order(surveys[sites==siteIDs[i1]])]
    d1 = as.dist(dmat[ind_surv1, ind_surv1])
    d2 = dist(ind_surv1)
    m = mantel(d1,d2, permutations = nperm)
    res[i1, "r"] = m$statistic
    res[i1, "p.value"] = m$signif
  }
  return(res)
}