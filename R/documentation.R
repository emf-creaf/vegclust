#'  Wetland vegetation data set
#'  
#'  Vegetation of the Adelaide river alluvial plain (Australia). This data set was published by Bowman & Wilson (1987) and used in Dale (1988) to compare fuzzy classification approaches.
#'  
#' @name wetland
#' @docType data
#' @encoding UTF-8
#' @format A data frame with 41 sites (rows) and 33 species (columns). Abundance values are represented in abundance classes.
#' @source Bowman, D. M. J. S. and B. A. Wilson. 1986. Wetland vegetation pattern on the Adelaide River flood plain, Northern Territory, Australia. Proceedings of the Royal Society of Queensland 97:69-77.
#' @keywords data
#' 
#' @references 
#' Dale, M. B. 1988. Some fuzzy approaches to phytosociology. Ideals and instances. Folia geobotanica et phytotaxonomica 23:239-274.
#' 
#' @examples
#' data(wetland)
#' 
NULL

#' Synthetic vegetation data set with tree data
#' 
#' A synthetic data set used to illustrate the stratification of data originally collected on an individual basis (e.g. forest inventory).
#' 
#' @name treedata
#' @docType data
#' @encoding UTF-8
#' @format A data frame where each row corresponds to a different tree. Columns are plot code, species identity, tree height, tree diameter and cover value.
#' @keywords data
#' @seealso \code{\link{stratifyvegdata}}
#' 
NULL

#' Regeneration of Mediterranean vegetation data set
#' 
#' A stratified vegetation data set containing with several plot records laid to assess vegetation recovery three years after a wildfire. Collected in 2012 by Miquel De Caceres and Albert Petit in Horta de Sant Joan (Catalonia, Spain).
#' 
#' @name medreg
#' @docType data
#' @encoding UTF-8
#' @format An object of class \code{stratifiedvegdata} with 96 elements (plots), each of them consisting of a data.frame where rows correspond to species groups and columns correspond to vegetation strata. Abundance values are percentage cover.
#' @keywords data
#' @seealso \code{\link{CAP}}, \code{\link{plot.CAP}}, \code{\link{stratifyvegdata}}
#' 
NULL

