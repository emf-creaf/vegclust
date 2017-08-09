trajectoryPCoA<-function(D, sites, selsites, clsites=NULL, cl_colors=NULL, axes=c(1,2), lwd=1.5, length=0.05) {
  siteIDs = unique(sites)
  nsite = length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) {
    nsurveysite[i] = sum(sites==siteIDs[i])
  }
  selIDs = siteIDs[selsites]
  D2 =as.dist(as.matrix(D)[sites %in% selIDs, sites %in% selIDs])
  cmd_D2<-cmdscale(D2,eig=TRUE, add=TRUE, k=nrow(as.matrix(D2))-1)
  # R2 = sum(cmd_D2$eig[axes])/sum(cmd_D2$eig)
  x<-cmd_D2$points[,axes[1]]
  y<-cmd_D2$points[,axes[2]]
  plot(x,y, type="n", asp=1, xlab=paste0("PCoA ",axes[1]," (", round(100*cmd_D2$eig[axes[1]]/sum(cmd_D2$eig)),"%)"), 
       ylab=paste0("PCoA ",axes[2]," (", round(100*cmd_D2$eig[axes[2]]/sum(cmd_D2$eig)),"%)"))
  
  #Draw arrows
  cum=0
  sitesred = sites[sites %in% selIDs]
  for(i in 1:nsite) {
    if(selsites[i]) {
      ind_surv = which(sitesred==siteIDs[i])
      for(t in 1:(nsurveysite[i]-1)) {
        niini =ind_surv[t]
        nifin =ind_surv[t+1]
        if(!is.null(clsites)) {
          arrows(x[niini],y[niini],x[nifin],y[nifin], length=length,
                 col=cl_colors[as.numeric(clsites)][i], lwd=lwd)
        } else {
          arrows(x[niini],y[niini],x[nifin],y[nifin], length=length, lwd=lwd)
        }
      }
      cum = cum + nsurveysite[i]
    }
  }
  #Draw legend
  a = table(as.factor(clsites)[selsites])
  if(!is.null(clsites)) legend("topright", legend= names(a)[a>0], bty="n", lty=1, col=cl_colors[a>0], lwd=lwd)
  invisible(cmd_D2)
}
