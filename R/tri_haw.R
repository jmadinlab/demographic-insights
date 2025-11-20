# TRIMODAL

#library load
# library(stats)
library(tidyverse)
library(popbio)
library(gridExtra)

# optimization function for recruitment
fnc <- function(x, mat,lam=1, i=1, j=ncol(mat)) {
  mat[i, j] <- x
  return((lam - lambda(mat))^2)
}

# Load Trimodal Data, Prep for loop
gtri <- read.csv("data/trimodal/growth.csv")
stri <- read.csv("data/trimodal/survival.csv")
stri <- stri[!is.na(stri$surv),]

# Add Hawaii Taxa ---------------------------------------------------------
# ADD and morph Hawaii Taxa
load("data/hawaii/Hawaii_NCRMP/Colony_Data_20210917_edited.rdata")
load("data/hawaii/Hawaii_NCRMP/HA_MA_AS_ColonyTransition_withAnnSurv.rdata")

#Species Names
SPLU <- c("Montipora capitata", "Montipora flabelata", "Montipora sp.", "Montipora patula", "Pocillopora grandis", "Porites lichen", "Porites lobata", "Porites lutea", "Pocillopora meandrina", "Pocillopora sp.", "Porites sp.")
names(SPLU) <- sort(unique(ColonyLevel$Spec_Code))

#Date Formats
ColTrans$StartingDate <- ymd(ColTrans$StartingDate)
ColonyLevel$StartingDate <- ymd(ColonyLevel$StartingDate)
ColTrans$EndingDate <- ymd(ColTrans$EndingDate)
ColonyLevel$EndingDate <- ymd(ColonyLevel$EndingDate)

ColTrans_Sel <- ColTrans %>% filter(REGION%in%c("MHI","NWHI")) %>%
  select(Site_Genet,StartingDate,EndingDate,AnnSurvProbs)
ColonyLevel_ADD <- ColonyLevel %>% 
  left_join(ColTrans_Sel,
            by=join_by(ColonyID==Site_Genet,StartingDate,EndingDate))

gshaw <- ColonyLevel_ADD %>% 
  filter(DataOrError=="DATA",StartingSize>pi,# minimum 2 cm colony size
         Spec_Code%in%c("MCAP","MPAT","PLIC","PLOB","PLUT","PMEA")) %>%
  mutate(Spec_Name=SPLU[Spec_Code]) %>% 
  select(ColonyID,Spec_Code,Spec_Name,StartingYear,StartingSize,AnnualEndingSize_E,Survival,AnnSurvProbs) %>% 
  rename(c(colony_id=ColonyID,
           spp=Spec_Code,
           species=Spec_Name,
           year=StartingYear,
           area_cm2=StartingSize,area_cm2_next=AnnualEndingSize_E,surv=Survival,Pasurv=AnnSurvProbs))
head(gshaw)

#Now split into the growth and surv datasets
ghaw <- gshaw %>% filter(surv==1) %>% select(!c(surv,Pasurv))
head(ghaw)

# Here I'm changing to use the continuous probability for annual survival estimates
shaw <- gshaw %>%
  select(!c(area_cm2_next,surv)) %>% 
  rename(surv=Pasurv)
head(shaw)

# Merge Hawaii and GBR data for growth and survival
gall <- rbind(gtri, ghaw)
sall <- rbind(stri, shaw)

# Move on with analysis ---------------------------------------------------
spp <- unique(sort(gall$species))
nspp <- length(spp)
nspp

# read params with assumed HA taxa params
params <- read.csv("data/params.csv")

# Calculate size mid points
params$rec <- params$rec_juv / 2
params$juv <- (params$rec_juv + params$juv_adu) / 2
params$adu <- params$juv_adu + (params$juv_adu - params$juv)
head(params)

### FECUNDITY TRIMODAL
ftri <- read.csv("data/trimodal/fecundity.csv")
pp <- params[params$species %in% unique(ftri$species),]
pp <- pp[c("species", "juv_adu")]

ftri <- merge(ftri, pp)

tail(ftri)

ggplot(ftri, aes(x=area_cm2, y=reproductive)) +
  geom_jitter(height=0.02, width=0.05, color="grey") +
  geom_smooth(method="glm",method.args=list(family="binomial")) +
  scale_x_log10() +
  geom_hline(yintercept=0.75, lty=2) +
  facet_wrap(~species) +
  ylim(c(0, 1))  +
  geom_vline(data = . %>% group_by(species) %>% summarise(vl=mean(juv_adu)), 
             aes(xintercept=vl), lty=2)

ggsave("figs/figureS1-trifecund.png", width = 8, height = 2.5)

# LOOP
plot_list <- list() 
store <- list()

for (i in spp) {

  mat <- matrix(0, 3, 3)

  gdat <- gall[gall$species == i,]
  sdat <- sall[sall$species == i,]

  rec_juv <- params$rec_juv[params$species==i]
  juv_adu <- params$juv_adu[params$species==i]

  ###
  plot_list[[i]] <- ggplot(sdat, aes(y = surv, x = area_cm2)) +
    geom_jitter(height=0.05, color="grey") +
    geom_smooth(se = FALSE) +
    scale_x_log10() +
    geom_vline(xintercept=rec_juv, linetype=1) +
    geom_vline(xintercept=juv_adu, linetype=2) +
    ggtitle(i, ) + 
    theme(plot.title = element_text(size = 10))
  ###  
  

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
  
  print(lambda(mat))
  
  opt <- optimise(f=fnc, interval=c(0, 10000), mat=mat, maximum = FALSE)
  mat[1, ncol(mat)] <- opt$minimum
  print(lambda(mat))
  
  store[[i]] <- mat

}

gg <- grid.arrange(grobs=plot_list, ncol=3)

ggsave("figs/figureS2-allsurv.png", gg, width = 8, height = 7)

save(store, file="output/mats_tri_haw.RData")

