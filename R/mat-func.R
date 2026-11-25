##### Libraries

library(tidyverse)
library(popbio)
library(gridExtra)
library(Rcompadre)
library(tidyverse)
library(patchwork)
library(ggsci)



##### OPTIMIZATION function for recruitment

# EXAMPLE
# mat <- matfnc("Acropora cf. digitifera")
# opt <- optimise(f=recfnc, interval=c(0, 10000), mat=mat, maximum = FALSE, lam=1)
# mat[1, ncol(mat)] <- opt$minimum
# lambda(mat)

recfnc <- function(x, mat,lam=1, i=1, j=ncol(mat)) {
  mat[i, j] <- x
  return((lam - lambda(mat))^2)
}

mfnc <- function(mat, x, lam=1, i=1, j=ncol(mat)) {
  mat[i, j] <- x
  (lam - lambda(mat))^2
}

afnc <- function(mat, x, stg=1, lam=1.1) {
  ss <- stable.stage(mat) * 1
  ss[stg] <- ss[stg] + x
  (lam - sum(mat %*% ss))^2
}

mfnc_prop <- function(mat, x, lam=1, i=1, j=ncol(mat), prop=1) {
  mat[i, j] <- x * prop + mat[i, j] * (1 - prop)
  (lam - lambda(mat))^2
}

##### MATRIX MAKING function for continuous data

# EXAMPLE
# matfnc("Acropora cf. digitifera")

matfnc <- function(spp) {

  mat <- matrix(0, 3, 3)
  
  gdat <- gall[gall$species == spp,]
  sdat <- sall[sall$species == spp,]
  
  rec_juv <- params$rec_juv[params$species==spp]
  juv_adu <- params$juv_adu[params$species==spp]
  
  r_s <- mean(sdat$surv[sdat$area_cm2 < rec_juv])
  
  mat[1, 1] <- nrow(gdat[gdat$area_cm2 < rec_juv & gdat$area_cm2_next < rec_juv,]) / nrow(gdat[gdat$area_cm2 < rec_juv,]) * r_s
  mat[2, 1] <- nrow(gdat[gdat$area_cm2 < rec_juv & gdat$area_cm2_next >= rec_juv & gdat$area_cm2_next < juv_adu,]) / nrow(gdat[gdat$area_cm2 < rec_juv,]) * r_s
  mat[3, 1] <- nrow(gdat[gdat$area_cm2 < rec_juv & gdat$area_cm2_next >= juv_adu,]) / nrow(gdat[gdat$area_cm2 < rec_juv,]) * r_s
  
  j_s <- mean(sdat$surv[sdat$area_cm2 >= rec_juv & sdat$area_cm2 < juv_adu])
  
  mat[1, 2] <- nrow(gdat[gdat$area_cm2 >= rec_juv & gdat$area_cm2 < juv_adu & gdat$area_cm2_next < rec_juv,]) / nrow(gdat[gdat$area_cm2 >= rec_juv & gdat$area_cm2 < juv_adu,]) * j_s
  mat[2, 2] <- nrow(gdat[gdat$area_cm2 >= rec_juv & gdat$area_cm2 < juv_adu & gdat$area_cm2_next >= rec_juv & gdat$area_cm2_next < juv_adu,]) / nrow(gdat[gdat$area_cm2 >= rec_juv & gdat$area_cm2 < juv_adu,]) * j_s
  mat[3, 2] <- nrow(gdat[gdat$area_cm2 >= rec_juv & gdat$area_cm2 < juv_adu & gdat$area_cm2_next >= juv_adu,]) / nrow(gdat[gdat$area_cm2 >= rec_juv & gdat$area_cm2 < juv_adu,]) * j_s
  
  a_s <- mean(sdat$surv[sdat$area_cm2 >= juv_adu])
  
  mat[1, 3] <- nrow(gdat[gdat$area_cm2 >= juv_adu & gdat$area_cm2_next < rec_juv,]) / nrow(gdat[gdat$area_cm2 >= juv_adu,]) * a_s
  mat[2, 3] <- nrow(gdat[gdat$area_cm2 >= juv_adu & gdat$area_cm2_next >= rec_juv & gdat$area_cm2_next < juv_adu,]) / nrow(gdat[gdat$area_cm2 >= juv_adu,]) * a_s
  mat[3, 3] <- nrow(gdat[gdat$area_cm2 >= juv_adu & gdat$area_cm2_next >= juv_adu,]) / nrow(gdat[gdat$area_cm2 >= juv_adu,]) * a_s

  if (any(mat == "NaN")) {
    return(FALSE)
  } else {
    return(mat)
  }
}
