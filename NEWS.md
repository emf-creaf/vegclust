-------------------------------
 NEWS for R Package "vegclust"
-------------------------------

# Version 2.0.0
- IMPORTANT: Functions for Community Trajectory Analysis were moved to package 'ecotraj', available at GitHub (https://github.com/emf-creaf/ecotraj/) and soon also in CRAN.

# Version 1.8.0
- Update of function 'trajectoryAngles2D' by Anthony Sturbois

# Version 1.7.9
- New option to 'trajectoryLengths' and 'trajectoryAngles' to calculate distances and angles relative to initial survey
- New functions 'trajectoryLengths2D' and 'trajectoryAngles2D' by Anthony Sturbois
- Bug correction in vegdiststruct
- Update of CTA vignette (suggestions P. Legendre)

# Version 1.7.8
- New option to 'trajectoryPCoA' and 'trajectoryPlot' to draw survey labels

# Version 1.7.7
- Bug correction 'stratifyvegdata'

# Version 1.7.6
- Adapt to Rcpp changes

# Version 1.7.5
- Now trajectory angles are measured between the direction of the first segment and the direction of the second segment.
- Now mean and sd of trajectory angles are calculated using functions of package 'circular'
- Parameter 'add' added to community trajectory analysis to deactivate constant addition in triplets violating the triangle inequality.
- Update "centerTrajectories" according to Anderson (2017)
- Improved CTA vignette

# Version 1.7.4
- Improvement of parameter checking in function 'stratifyvegdata'
- Documentation of avoca data set
- New function 'trajectoryPlot'

# Version 1.7.3
- New function 'trajectoryConvergence'
- New function 'trajectoryDirectionality'
- New function 'centerTrajectories'

# Version 1.7.2
- New function 'trajectoryAngles'.
- New function 'trajectoryProjection'.
- New vignette for trajectory analysis.

# Version 1.7.1
- Package encoding. Changed package and author description.

# Version 1.7.0
- NEW feature: Trajectory analysis. Functions 'trajectoryDistances', 'trajectoryLengths' and 'trajectoryPCoA.
- Added option to 'trajectoryDistances' to allow different symmetrization functions (e.g. min) in addition to 'mean'.
