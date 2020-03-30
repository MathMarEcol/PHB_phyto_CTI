## Phenology comparisons ajani/imos/dakin&colefax
## Claire with Ant's comments and feedback
## March 2020

suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
  library(reshape)
  library(ggplot2)
  library(visreg)
  library(vegan)
  library(lubridate)
  library(patchwork)})

source("Harm.R")

## open and manage Ajani (1998-2009) data
ajani <- read.csv("Ajani_data.csv") %>% 
  dplyr::rename("Survey" = "SURVEY", "SampleDate" = "SAMPLE_DATE", "Taxon" = "TAXON_NAME",
                "AbundanceM3"="ABUNDANCE_M3", "TaxonGroup" = "TAXON_GROUP") %>%
  mutate(Date = ymd_hms(SampleDate),
         mon = month(SampleDate))  %>%
  mutate(Taxon = recode(Taxon, "Ditylum brightwellii > 40 µm width" = "Ditylum brightwellii",
                        "Ditylum brightwellii < 40 µm width" = "Ditylum brightwellii",
                        "Scrippsiella trochoidea" = "Scrippsiella spp.", #Penny not confident about this id
                        "Bacteriastrum furcatum" =	"Bacteriastrum spp." #Penny not confident about this id
                        )) 

# remove taxon where there are no counts
atax <- ajani %>% group_by(Taxon) %>% summarise(sums = sum(AbundanceM3, na.rm=TRUE)) %>%
  filter(sums>0) %>% droplevels()

ajani <- ajani %>% filter(Taxon %in% levels(atax$Taxon)) %>% droplevels() 

## open and manage IMOS data
imos <- read.csv("imos_data.csv", header = TRUE) %>% select(c(2:5)) %>%
  cast(SAMPLE_DATE ~ TAXON_NAME + TAXON_GROUP, mean, fill=0) %>%
  gather("Taxon", "AbundanceM3", -SAMPLE_DATE) %>%
  mutate(Survey = 'IMOS',
         TaxonGroup = as.factor(sub(".*_", "", Taxon)),
         Taxon = as.factor(sub("_.*", "", Taxon))) %>%
  dplyr::rename("SampleDate" = "SAMPLE_DATE") %>%
  mutate(Date= ymd_hms(SampleDate),
         mon = month(SampleDate)) %>% 
  mutate(Taxon = recode(Taxon, "Guinardia flaccida <150 µm length" = "Guinardia flaccida",
                                                      "Guinardia flaccida >150 µm length" = "Guinardia flaccida",
                                                       "Guinardia striata (with Richelia)" = "Guinardia striata",
                                                       "Lauderia annulata < 70 µm length" = "Lauderia annulata",
                                                       "Lauderia annulata > 70 µm length" = "Lauderia annulata",
                                                       "Ditylum brightwellii > 40 µm width" = "Ditylum brightwellii",
                                                       "Ditylum brightwellii < 40 µm width" = "Ditylum brightwellii",
                                                       "Dactyliosolen fragillissimus >150 µm length" = "Dactyliosolen fragillissimus",
                                                       "Dactyliosolen fragillissimus <150 µm length" = "Dactyliosolen fragillissimus",  
                                                       "Leptocylindrus mediterraneus (no flagellates)" = "Leptocylindrus mediterraneus", 
                                                       "Leptocylindrus mediterraneus (with flagellates)" = "Leptocylindrus mediterraneus",
                                                       "Proboscia alata ~ 50 µm cell width" = "Proboscia alata" ,
                                                       "Asterionellopsis spp. ~50 µm length" = "Asterionellopsis glacialis",
                        "Nitzschia cf. bicapitata" = 'Nitzschia bicapitata',
                        "Chaetoceros peruvianus < 40 µm cell width" = "Chaetoceros peruvianus",
                        "Chaetoceros peruvianus > 40 µm cell width" = "Chaetoceros peruvianus"))  %>% droplevels()


## open and manage Dakin and Colefax (DC) data
dc <- read.csv("dc_data.csv", header = TRUE) %>%
    dplyr::rename("Survey" = "SURVEY", "SampleDate" = "SAMPLE_DATE", "Taxon" = "TAXON_NAME",
                       "AbundanceM3"="ABUNDANCE_M3", "TaxonGroup" = "TAXON_GROUP") %>%
  mutate(Date= ymd_hms(SampleDate),
         mon = month(SampleDate)) %>%
  mutate(Taxon = recode(Taxon, "Ditylum brightwellii > 40 µm width" = "Ditylum brightwellii",
                        "Ditylum brightwellii < 40 µm width" = "Ditylum brightwellii"))

## open and manage OSL, Hallegraeff, Ajani97 (Ajani2) data
ph <- read.csv("ph_data.csv", header = TRUE) %>%
  dplyr::rename("Survey" = "SURVEY", "SampleDate" = "SAMPLE_DATE_UTC", "Taxon" = "TAXON_NAME",
                "AbundanceM3"="ABUNDANCE_M3", "TaxonGroup" = "TAXON_GROUP") %>%
  mutate(Date= ymd_hms(SampleDate),
         mon = month(SampleDate)) %>%
  mutate(Taxon = recode(Taxon, "Ditylum brightwellii > 40 µm width" = "Ditylum brightwellii",
                        "Ditylum brightwellii < 40 µm width" = "Ditylum brightwellii",
                        "Scrippsiella trochoidea" = "Scrippsiella spp.")) #Penny not confident about this id

# DC: 37um net, surface, 1930-1931
# OSL: Dorn sampler, 0-100m, 1965-1968
# Hallegraeff: 6L water sampler 0-100m, 1978-1979
# Ajani97: 37um net, 0-100m, 1997-1998
# Ajani98: 20um net, 0-50m, 1998-2009
# IMOS: niskin, 0-50m, 2009-2013
# see paper for different sampling methods and depths

## combine all the datasets
aid <- rbind(ajani, imos, dc, ph)
#write.csv(aid, "all_ph_data.csv")

#################################################################
## sampling frequency

dates <- aid[,c(1,6,7)] %>% 
  mutate(year = year(Date))
dates$Survey <- factor(dates$Survey, levels(dates$Survey)[c(3, 6, 5, 4, 1, 2)])

samps <- ggplot(dates, aes(year, mon, colour=Survey)) + geom_point() +
  theme_bw(base_size = 12) + labs(x="Year", y="Month") + 
  scale_y_reverse(breaks = seq(1,12,length.out = 12), label = c("J","F","M","A","M","J","J","A","S","O","N","D")) 
#x11(width=11, height=8)
samps 

ggsave("PH_samples.png", samps, dpi = 1200)

#################################################################################################

# how similiar are these sample compositions, look with an MDS
fg <- aid %>% group_by(Date, Survey, Taxon) %>% 
  summarise(abund = sum(AbundanceM3, na.rm = TRUE)) %>%
  cast(Date + Survey ~ Taxon, mean, fill=0) 

fgm <- fg %>% select(3:326) %>% as.matrix(fgm)

taxaSum <- colSums(fgm) # Take the sum of each column (taxa). There should be no 0 for a species now.
min(taxaSum)

summary(fgm)

fg_mds <- metaMDS(fgm, distance = "bray", autotransform = TRUE, k= 2, maxit=1000, 
                  try = 20, trymax = 50)

Pal <- c(rainbow(n=6)) # Define the palette
x11(width=11, height=8)
grp <- fg$Survey
plot(fg_mds$points, pch = 16, col=Pal[grp]) #plot so 3 outlying samples don't show
#text(fg_mds, display = 'sites')
legend('topleft', legend = levels(fg$Survey), pch = 16, col = Pal,  y.intersp = 0.7)
# all pretty different according to survey
# 1998-03-03 and 1998-03-07, 2019-04-09, 2015-07-23 outlying samples

ggsave("PH_mds.png", dpi = 1200)

#############################################################################################

# psti_id can be generated from STI.R or used from the saved file.
psti_id <- read.csv("pSTI.csv") 

# remove groups of taxa with suspect id's and silicos that DC did not count
mintemp <- 6
aidSTI <- aid  %>% 
  filter(!grepl("Thalassiosira s", Taxon)  & !grepl("Corethron", Taxon)  & !grepl("type", Taxon)  &
           !grepl("Thalassionema", Taxon) & !grepl("Cerat", Taxon)  &
           !grepl("Silicof", TaxonGroup) & !grepl("Chaetoceros spp", Taxon) &
           !grepl(".group", Taxon) & !grepl("Pseudo", Taxon) & !grepl("Thalassiosira g", Taxon)  & 
           !grepl("Lauder", Taxon) & !grepl("Detonu", Taxon)  & # these id's are confused but the STI's are similiar, remove anyway.
           !grepl("Eucamp", Taxon) & !grepl("Tripos p", Taxon)  & 
           !grepl("Nitzschia bicapitata", Taxon)) # STI is highly skewed by SO data, not sure it is the same species everywhere 

# assign species an STI
aidSTI <- left_join(aidSTI, psti_id, by = "Taxon") %>% 
  filter(!is.na(STI) & STI > mintemp) %>% 
  mutate(Taxon = as.factor(Taxon),
         Period = ifelse(Survey=='Ajani', "1997-2009", 
                         ifelse(Survey=='Halle', "1978-1979", 
                                ifelse(Survey=='OSL', "1965-1966", 
                                       ifelse(Survey=='Ajani2', "1997-2009",  # works better if the Ajani studies are combined
                                              ifelse(Survey == 'IMOS', "2009-2018", '1931-1932'))))),
  Sampler = ifelse(Survey=='Ajani', "net", 
                  ifelse(Survey=='Halle', "net", 
                         ifelse(Survey=='OSL', "bottle", 
                                ifelse(Survey=='Ajani2', "net",  # works better if the Ajani studies are combined
                                       ifelse(Survey == 'IMOS', "bottle", 'net')))))) %>%
  subset(AbundanceM3 > 0) %>%
  mutate(Period=as.factor(Period)) %>%
  subset(Survey != 'Halle') %>% # so few species in this dataset that it is non representative
  droplevels()

## Calculate CTI for each sample

# Case 1- use species that occur in every survey, too few n, so occur in all but one of the surveys

n=3 #no of surveys taxa to be found in 

taxa <- distinct(aidSTI, Period, Taxon) %>% group_by(Taxon) %>% summarise(freq = n()) %>%
  filter(freq>=n)  %>% droplevels() #find taxa present in n surveys
aidSTIt <- aidSTI %>% filter(Taxon %in% levels(taxa$Taxon)) %>% droplevels() 
noSTI <- distinct(subset(aidSTIt[,c(3,8)], is.na(aidSTIt$STI))) # check what has no STI , should be higher taxonomic groups or rare species

aidctit <- aidSTIt %>% filter(!is.na(STI))  %>%
  group_by(Period, SampleDate) %>% summarise(cti = sum(AbundanceM3*STI, na.rm=TRUE)/sum(AbundanceM3, na.rm=TRUE),
                                                                                        n = n()) %>%
  filter(n>1) %>% # remove samples where CTI only made up of 1 species
  mutate(Dates = ymd_hms(SampleDate),
         Month = month(Dates),
         Monh = Month/12*2*base::pi,
         Year = year(Dates)) 

# Use lm to see how CTI changes with years after removing the seasonal cycle

aidctit$Period <- as.factor(aidctit$Period)
lmCTIaidt <- lm(cti ~ Period + Harm(Monh, k=1), data=aidctit) #k=1 not significant, k=2 unrealistic
summary(lmCTIaidt)
anova(lmCTIaidt)
#x11(width=11, height=8)
par(mfrow=c(2,2))
plot(lmCTIaidt) # not too bad

Pert <- visreg(lmCTIaidt, "Period", rug = FALSE, gg = TRUE) + theme_bw(base_size = 14)+
  ylab(bquote("Community Temperature Index ("* degree *"C)")) +
  xlab("Period")
mont <- visreg(lmCTIaidt, "Monh", rug = FALSE, gg = TRUE) + theme_bw(base_size = 14)+
  ylab(bquote("Community Temperature Index ("* degree *"C)")) +
  xlab("Month") + scale_x_continuous(breaks = seq(0.52,6.28,length.out = 12), label = c("J","F","M","A","M","J","J","A","S","O","N","D")) 

#x11(width = 12, height = 5)
pt <- Pert + mont
pt
ggsave("PHB_CTI_SurveysPhyto.png", pt, dpi = 1200)

## look at relative abundances over periods for species occuring in all surveys
addAbun <- aidSTIt  %>% filter(!is.na(STI)) %>% 
  group_by(Period) %>% summarise(sums = sum(AbundanceM3))
pie <- inner_join(aidSTIt, addAbun, by=c("Period")) %>% filter(!is.na(STI)) %>% 
  mutate(pies = AbundanceM3/sums) %>% group_by(Period, Taxon) %>% summarise(sums = sum(pies))
pie2 <- pie %>% mutate(TaxonT = ifelse(sums>0.01, as.character(Taxon), "Other")) %>% 
  group_by(TaxonT, Period) %>% summarise(tots = sum(sums)) %>% complete(TaxonT, Period) %>%
  mutate(tots = ifelse(is.na(tots), 0, tots)) 

pie2$TaxonT <- as.factor(pie2$TaxonT)

#x11(width = 11, height = 8)
percPlot <- ggplot(data=pie2, aes(Period, tots)) + geom_col(aes(fill=TaxonT)) + theme_bw()
percPlot
ggsave("PHB_prop_3surveysphyto.png", percPlot, dpi = 1200)
##biggest change is relative decrease in A. glacilis and increase in L. danicus


#####################################################################################
# Case 2 - use species that occur in IMOS more than once and at least one other survey

taxa4 <- aidSTI %>% 
  mutate(Survey = ifelse(Survey == 'IMOS', 'IMOS', 'Other')) %>%
  group_by(Survey, STI, Taxon) %>% summarise(freq = n()) %>%
  filter(freq>1) %>% group_by(STI, Taxon) %>% summarise(freq = n()) %>%
  filter(freq>1) %>% droplevels()

#taxa4 <- aidSTI %>% group_by(Survey, STI, Taxon) %>% summarise(freq = n()) %>%
#  filter(freq>1) %>% filter(Survey == 'IMOS') %>% droplevels()
aidSTI4 <- aidSTI %>% filter(Taxon %in% levels(taxa4$Taxon)) %>% droplevels() 
noSTI4 <- distinct(subset(aidSTI4[,c(3,8)], is.na(aidSTI4$STI))) # check what has no STI , should be higher taxonomic groups or rare species

aidcti4 <- aidSTI4 %>% filter(!is.na(STI))  %>% 
  group_by(Period, Sampler, SampleDate) %>% summarise(cti = sum(AbundanceM3*STI, na.rm=TRUE)/sum(AbundanceM3, na.rm=TRUE),
                                             n = n()) %>%
  filter(n>1) %>% # remove samples where CTI only made up of 1 species
  mutate(Dates = ymd_hms(SampleDate),
         Month = month(Dates),
         Monh = Month/12*2*base::pi,
         Year = year(Dates)) 

# Use lm to see how CTI changes with years after removing the seasonal cycle (with SAMPLER, SAMPLER IS NOT SIGNIFICANT)

aidcti4$Period <- as.factor(aidcti4$Period)
lmCTIaid4 <- lm(cti ~ Period + Harm(Monh, k=1), data=aidcti4) #k2 makes month wobbly, k=1 insignificant
summary(lmCTIaid4)
anova(lmCTIaid4)
#x11(width=11, height=8)
par(mfrow=c(2,2))
plot(lmCTIaid4) # Not real good

# month term is barely significant
# IMOS and Ajani are significantly different from D&C but not from each other
Per4 <- visreg(lmCTIaid4, "Period", rug = FALSE, gg = TRUE) + theme_bw(base_size = 14)+
  ylab(bquote("Community Temperature Index ("* degree *"C)")) + 
  xlab("Period") #+ scale_y_continuous(limits=c(15, 25))
mon4 <- visreg(lmCTIaid4, "Monh", rug = FALSE, gg = TRUE) + theme_bw(base_size = 14)+
  ylab(bquote("Community Temperature Index ("* degree *"C)")) +
  xlab("Month") + scale_x_continuous(breaks = seq(0.52,6.28,length.out = 12), label = c("J","F","M","A","M","J","J","A","S","O","N","D")) 
#Sam4 <- visreg(lmCTIaid4, "Sampler", rug = FALSE, gg = TRUE) + theme_bw(base_size = 14)+
#  ylab(bquote("Community Temperature Index ("* degree *"C)")) + 
#  xlab("sampler") 

#x11(width = 12, height = 5)
p4 <- Per4 + mon4 #+ Sam4 
p4
ggsave("PHB_CTI_IMOSphyto.png", p4, dpi = 1200)


###########################################################################
## look at relative abundances over periods for IMOS species
addAbun4 <- aidSTI4  %>% filter(!is.na(STI))  %>% 
          group_by(Period) %>% summarise(sums = sum(AbundanceM3)) 
pie4 <- inner_join(aidSTI4, addAbun4, by=c("Period")) %>% filter(!is.na(STI))  %>% 
  mutate(pies = AbundanceM3/sums) %>% group_by(Period, Taxon) %>% summarise(sums = sum(pies)) 
pie4 <- pie4 %>% mutate(TaxonT = ifelse(sums>0.01, as.character(Taxon), "Other")) %>% 
  group_by(TaxonT, Period) %>% summarise(tots = sum(sums)) %>% complete(TaxonT, Period) %>%
  mutate(tots = ifelse(is.na(tots), 0, tots)) 

pie4$TaxonT <- as.factor(pie4$TaxonT)

#x11(width = 11, height = 8)
percPlot4 <- ggplot(data=pie4, aes(Period, tots)) + geom_col(aes(fill=TaxonT)) + theme_bw()
percPlot4
ggsave("PHB_prop_IMOSsurveysphyto.png", percPlot4, dpi = 1200)
##biggest change is decrease in A. glacilis and increase in L. danicus

###########################################################################
## look at relative abundances over periods for all species counted in a sample (true relative proportions)
addAbunAll <- aidSTI  %>% filter(!is.na(STI))  %>% 
  group_by(Period) %>% summarise(sums = sum(AbundanceM3)) 
pieAll <- inner_join(aidSTI, addAbunAll, by=c("Period")) %>% filter(!is.na(STI))  %>% 
  mutate(pies = AbundanceM3/sums) %>% group_by(Period, Taxon) %>% summarise(sums = sum(pies)) 
pieAll <- pieAll %>% mutate(TaxonT = ifelse(sums>0.01, as.character(Taxon), "Other")) %>% 
  group_by(TaxonT, Period) %>% summarise(tots = sum(sums)) %>% complete(TaxonT, Period) %>%
  mutate(tots = ifelse(is.na(tots), 0, tots)) 

pieAll$TaxonT <- as.factor(pieAll$TaxonT)

#x11(width = 11, height = 8)
percPlotAll <- ggplot(data=pieAll, aes(Period, tots)) + geom_col(aes(fill=TaxonT)) + theme_bw()
percPlotAll
ggsave("PHB_prop_Allphyto.png", percPlotAll, dpi = 1200)
##biggest change is decrease in A. glacilis and increase in L. danicus



