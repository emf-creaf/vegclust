library(vegclust)

## Loads data  
data(wetland)

## This equals the chord transformation 
wetland.chord <- as.data.frame(sweep(as.matrix(wetland), 1, 
                                     sqrt(rowSums(as.matrix(wetland)^2)), "/"))

test_that("Paritioning can be performed",{
  expect_s3_class(vegclust(wetland.chord, mobileCenters=3, m = 1.2, dnoise=0.75, method="NC", nstart=2),
                  "vegclust")
  expect_s3_class(vegclust(wetland.chord, mobileCenters=3, m = 1.2, method="FCM", nstart=2),
                  "vegclust")
  expect_s3_class(vegclust(wetland.chord, mobileCenters=3, method="KM", nstart=2),
                  "vegclust")
})