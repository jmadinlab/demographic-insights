# TRIMODAL & HAWAII
params <- read.csv("data/params.csv")

# Load Trimodal Data
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
spp <- unique(sort(gall$species))
nspp <- length(spp)
nspp

### FECUNDITY TRIMODAL
ftri <- read.csv("data/trimodal/fecundity.csv")
pp <- params[params$species %in% unique(ftri$species),]
pp <- pp[c("species", "juv_adu")]

ftri <- merge(ftri, pp)

ggplot(ftri, aes(x=area_cm2, y=reproductive)) +
  geom_jitter(height=0.02, width=0.05, color="grey") +
  geom_smooth(method="glm",method.args=list(family="binomial")) +
  scale_x_log10() +
  geom_hline(yintercept=0.75, lty=2) +
  facet_wrap(~species, ncol=3) +
  ylim(c(0, 1))  +
  labs(x="Area (cm^2)", y="Probability of being reproductive") +
  geom_vline(data = . %>% group_by(species) %>% summarise(vl=mean(juv_adu)), 
             aes(xintercept=vl), lty=2)

ggsave("figs/figureS1-trifecund.png", width = 8, height = 5)

### ALL SURZVIVAL

temp <- left_join(sall, params[c("species", "rec_juv", "juv_adu")], by = "species")
ggplot(temp, aes(y = surv, x = area_cm2)) +
  geom_jitter(height=0.05, color="grey") +
  geom_smooth(se = FALSE) +
  scale_x_log10() +
  geom_vline(aes(xintercept=rec_juv), linetype=1) +
  geom_vline(aes(xintercept=juv_adu), linetype=2) +
  labs(x="Area (cm^2)", y="Probability of survival") +
  facet_wrap(~species, ncol=3) +
  theme(plot.title = element_text(size = 10))

ggsave("figs/figureS2-allsurv.png", width = 8, height = 7)

# NORMAL RECRUITMENT

spp <- unique(sort(gall$species))

store_norm <- lapply(setNames(spp, spp), matfnc)
length(store_norm)

params$rec_juv <- 100 # pi * (2.5^2) # 5 cm diameter coral as per reviewer
store_100 <- lapply(setNames(spp, spp), matfnc)
names(store_100) <- paste0(names(store_100), "_100")

params$rec_juv <- 50 # pi * (2.5^2) # 5 cm diameter coral as per reviewer
store_50 <- lapply(setNames(spp, spp), matfnc)
names(store_50) <- paste0(names(store_50), "_50")

params$rec_juv <- 20 # pi * (2.5^2) # 5 cm diameter coral as per reviewer
store_20 <- lapply(setNames(spp, spp), matfnc)
names(store_20) <- paste0(names(store_20), "_20")

params$rec_juv <- 10 # pi * (2.5^2) # 5 cm diameter coral as per reviewer
store_10 <- lapply(setNames(spp, spp), matfnc)
names(store_10) <- paste0(names(store_10), "_10")

params$rec_juv <- 5 # pi * (2.5^2) # 5 cm diameter coral as per reviewer
store_5 <- lapply(setNames(spp, spp), matfnc)
names(store_5) <- paste0(names(store_5), "_5")



# 5 cm RECRUITMRNT
# gall <- gall[!(gall$species %in% c("Acropora cytherea", "Acropora humilis", "Acropora hyacinthus", "Acropora intermedia", "Acropora millepora", "Acropora nasuta", "Acropora robusta", "Acropora spathulata")),]
# sall <- sall[!(sall$species %in% c("Acropora cytherea", "Acropora humilis", "Acropora hyacinthus", "Acropora intermedia", "Acropora millepora", "Acropora nasuta", "Acropora robusta", "Acropora spathulata")),]

# spp <- unique(sort(gall$species))
# length(spp)
# params$rec_juv <- 20 # pi * (2.5^2) # 5 cm diameter coral as per reviewer
# 
# store_5cm <- lapply(setNames(spp, spp), matfnc)
# length(store_5cm)
# 
# names(store_5cm) <- paste0(names(store_5cm), "_5cm")
