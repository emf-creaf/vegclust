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
    if(!is.null(surveysred)) ind_surv = ind_surv[order(surveysred[sitesred==siteIDs[i]])]
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
