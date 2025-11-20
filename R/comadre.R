# Comadre

library(Rcompadre)
library(popbio)

# optimization function for recruitment
fnc <- function(x, mat,lam=1, i=1, j=ncol(mat)) {
  mat[i, j] <- x
  return((lam - lambda(mat))^2)
}

# Comadre start

comadre <- cdb_fetch("data/comadre/COMADRE_v.4.23.3.1.RData")

x <- subset(comadre, Order == "Scleractinia")
x <- subset(x, MatrixTreatment=="Unmanipulated")

load("output/mats_tri_haw.RData")
mats <- store
# mats <- list()

## Matrix navigator

for (i in 1:nrow(x)) {
  mat <- x$mat[[i]]@matA
  print(i)
  print(mat)
  print("-------------------------------")
}

store_size <- data.frame()
meta <- data.frame()

for (i in 1:nrow(x)) {
  spp <- x$SpeciesAccepted[i]
  if (spp %in% names(mats)) {
    spp <- paste0(x$SpeciesAccepted[i], " 2")
    if (spp %in% names(mats)) {
      spp <- paste0(x$SpeciesAccepted[i], " 3")
      if (spp %in% names(mats)) {
        spp <- paste0(x$SpeciesAccepted[i], " 4")
        if (spp %in% names(mats)) {
          spp <- paste0(x$SpeciesAccepted[i], " 5")
        }
      }
    }
  }
  
  mat <- x$mat[[i]]@matA
  siz <- x$mat[[i]]@matrixClass$MatrixClassAuthor
  
  if (nrow(mat) == 3) {
    if (all(mat[cbind(c(1, 2, 2, 3, 3), c(1, 1, 2, 2, 3))] > 0)) {
      if (all(apply(mat, 2, sum) <= 1)) {
        if (!any(is.na(mat))) {
          print(paste0(i, "----", spp))
          opt <- optimise(f=fnc, interval=c(0, 10000), mat=mat)
          mat[1, ncol(mat)] <- opt$minimum
          print(lambda(mat))
        
          mats[[spp]] <- mat
          
          store_size <- rbind(store_size, data.frame(species=spp, rec=siz[1], juv=siz[2], adult=siz[3]))
          meta <- rbind(meta, data.frame(i=i, species=spp, authors=x$Authors[i], year=x$YearPublication[i], doi=x$DOI_ISBN[i]))  
          
        }
      }
    }
  }
}

write.csv(store_size, "output/camadre_size.csv", row.names = FALSE)
write.csv(meta, "output/camadre_meta.csv", row.names = FALSE)

save(mats, file="output/mats.RData")
