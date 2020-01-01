## Phenology comparisons ajani/imos/dakin&colefax
## Claire with Ant's comments and feedback
## Dec 2019

suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
  library(reshape)
  library(ggplot2)
  library(visreg)
  library(lubridate)
  library(patchwork)})

source("Harm.R")

## open and manage Ajani data
ajani <- read.csv("Ajani_data.csv") %>% 
  dplyr::rename("Survey" = "SURVEY", "SampleDate" = "SAMPLE_DATE", "Taxon" = "TAXON_NAME",
                "AbundanceM3"="ABUNDANCE_M3", "TaxonGroup" = "TAXON_GROUP") %>%
  mutate(Date = ymd_hms(SampleDate),
         mon = month(SampleDate))  %>%
  mutate(Taxon = recode(Taxon, "Ditylum brightwellii > 40 µm width" = "Ditylum brightwellii",
                        "Ditylum brightwellii < 40 µm width" = "Ditylum brightwellii",
                        "Scrippsiella trochoidea" = "Scrippsiella spp.", #Penny not confident about this id
                        "Bacteriastrum furcatum" =	"Bacteriastrum spp.", #Penny not confident about this id
                        #"Nitzschia bicapitata" = 'Nitzschia spp.' # STI is highly skewed by SO data, not sure it is the same id 
                        )) 

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

## open and manage DC data
dc <- read.csv("dc_data.csv", header = TRUE) %>%
    dplyr::rename("Survey" = "SURVEY", "SampleDate" = "SAMPLE_DATE", "Taxon" = "TAXON_NAME",
                       "AbundanceM3"="ABUNDANCE_M3", "TaxonGroup" = "TAXON_GROUP") %>%
  mutate(Date= ymd_hms(SampleDate),
         mon = month(SampleDate)) %>%
  mutate(Taxon = recode(Taxon, "Ditylum brightwellii > 40 µm width" = "Ditylum brightwellii",
                        "Ditylum brightwellii < 40 µm width" = "Ditylum brightwellii"))

## combine all the datasets
aid <- rbind(ajani, imos, dc)
#write.csv(aid, "all_ph_data.csv")

#Try CTI approach
aidSTI <- aid  %>% # remove groups of taxa with suspect id's and silicos that DC did not count
  filter(!grepl("Thalassiosira s", Taxon)  & !grepl("Corethron", Taxon)  &
           !grepl("Thalassionema", Taxon) & 
           !grepl("Silicof", TaxonGroup) & !grepl("Chaetoceros spp", Taxon) &
           !grepl(".group", Taxon) & !grepl("Pseudo", Taxon)  ) 

psti_id <- read.csv("pSTI.csv") 

aidSTI <- left_join(aidSTI, psti_id, by = "Taxon") %>% mutate(Taxon = as.factor(Taxon)) %>% 
  subset(AbundanceM3 > 0) 

## Calculate CTI for each sample

# Case 1 - use species that occur in every survey

taxa <- distinct(aidSTI, Survey, Taxon) %>% group_by(Taxon) %>% summarise(freq = n()) %>%
  subset(freq>2)  %>% droplevels() #find taxa present in all surveys
aidSTIt <- aidSTI %>% filter(Taxon %in% levels(taxa$Taxon)) %>% droplevels() 
noSTI <- distinct(subset(aidSTIt[,c(3,8)], is.na(aidSTIt$STI))) # check what has no STI , should be higher taxonomic groups or rare species

aidctit <- aidSTIt %>% filter(!is.na(STI))  %>%
  group_by(Survey, SampleDate) %>% summarise(cti = sum(AbundanceM3*STI, na.rm=TRUE)/sum(AbundanceM3, na.rm=TRUE)) %>%
  mutate(Dates = ymd_hms(SampleDate),
         Month = month(Dates),
         Monh = Month/12*2*base::pi,
         Period = ifelse(Survey=='Ajani', "1998-2009", 
                         ifelse(Survey=='IMOS', "2009-2018", '1931-1932')),
         Year = year(Dates)) 

# Use lm to see how CTI changes with years after removing the seasonal cycle

aidctit$Period <- as.factor(aidctit$Period)
lmCTIaidt <- lm(cti ~ Period + Harm(Monh, k=1), data=aidctit) #k=1 not significant, k=2 unrealistic
summary(lmCTIaidt)
anova(lmCTIaidt)
#x11(width=11, height=8)
par(mfrow=c(2,2))
plot(lmCTIaidt) # Residuals aren't good, not normal
# many samples only have 1 species defining CTI so not very accurate, overlap of species isn't great

Pert <- visreg(lmCTIaidt, "Period", rug = FALSE, gg = TRUE) + theme_bw(base_size = 14)+
  ylab(bquote("Community Temperature Index ("* degree *"C)")) +
  xlab("Period")
mont <- visreg(lmCTIaidt, "Monh", rug = FALSE, gg = TRUE) + theme_bw(base_size = 14)+
  ylab(bquote("Community Temperature Index ("* degree *"C)")) +
  xlab("Month") + scale_x_continuous(breaks = seq(0.52,6.28,length.out = 12), label = c("J","F","M","A","M","J","J","A","S","O","N","D")) 

#x11(width = 12, height = 5)
pt <- Pert + mont
pt
ggsave("PHB_CTI_3SurveysPhyto.png", pt, dpi = 1200)

################################################
# Case 2 - use species that occur over two of the surveys

taxa2 <- distinct(aidSTI, Survey, Taxon) %>% group_by(Taxon) %>% summarise(freq = n()) %>%
  subset(freq>1)  %>% droplevels() 
aidSTI2 <- aidSTI %>% filter(Taxon %in% levels(taxa2$Taxon)) %>% droplevels() 
noSTI2 <- distinct(subset(aidSTI2[,c(3,8)], is.na(aidSTI2$STI))) # check what has no STI , should be higher taxonomic groups or rare species

aidcti2 <- aidSTI2 %>% filter(!is.na(STI))  %>%
  group_by(Survey, SampleDate) %>% summarise(cti = sum(AbundanceM3*STI, na.rm=TRUE)/sum(AbundanceM3, na.rm=TRUE),
                                             n = n()) %>%
  #filter(n>2) %>% # remove samples where CTI only made up of 2 or less species
  mutate(Dates = ymd_hms(SampleDate),
         Month = month(Dates),
         Monh = Month/12*2*base::pi,
         Period = ifelse(Survey=='Ajani', "1998-2009", 
                         ifelse(Survey=='IMOS', "2009-2018", '1931-1932')),
         Year = year(Dates)) 

# Use lm to see how CTI changes with years after removing the seasonal cycle

aidcti2$Period <- as.factor(aidcti2$Period)
lmCTIaid2 <- lm(cti ~ Period + Harm(Monh, k=1), data=aidcti2) 
summary(lmCTIaid2)
anova(lmCTIaid2)
#x11(width=11, height=8)
par(mfrow=c(2,2))
plot(lmCTIaid2) # Residuals are still dodgy, not normal

# month term is better, more believable
# periods are starting to be significantly different
Per2 <- visreg(lmCTIaid2, "Period", rug = FALSE, gg = TRUE) + theme_bw(base_size = 14)+
  ylab(bquote("Community Temperature Index ("* degree *"C)")) +
  xlab("Period")
mon2 <- visreg(lmCTIaid2, "Monh", rug = FALSE, gg = TRUE) + theme_bw(base_size = 14)+
  ylab(bquote("Community Temperature Index ("* degree *"C)")) +
  xlab("Month") + scale_x_continuous(breaks = seq(0.52,6.28,length.out = 12), label = c("J","F","M","A","M","J","J","A","S","O","N","D")) 

#x11(width = 12, height = 5)
p2 <- Per2 + mon2
p2
ggsave("PHB_CTI_2SurveysPhyto.png", p2, dpi = 1200)


################################################
# Case 3 - use species that occur in any of the surveys

aidcti3 <- aidSTI %>% filter(!is.na(STI))  %>% 
    group_by(Survey, SampleDate) %>% summarise(cti = sum(AbundanceM3*STI, na.rm=TRUE)/sum(AbundanceM3, na.rm=TRUE),
                                             n = n()) %>%
  filter(n>2) %>% # remove samples where CTI only made up of 2 or less species
  mutate(Dates = ymd_hms(SampleDate),
         Month = month(Dates),
         Monh = Month/12*2*base::pi,
         Period = ifelse(Survey=='Ajani', "1998-2009", 
                         ifelse(Survey=='IMOS', "2009-2018", '1931-1932')),
         Year = year(Dates)) 

# Use lm to see how CTI changes with years after removing the seasonal cycle

aidcti3$Period <- as.factor(aidcti3$Period)
lmCTIaid3 <- lm(cti ~ Period + Harm(Monh, k=1), data=aidcti3) #k2 makes month wobbly, k=1 insignificant
summary(lmCTIaid3)
anova(lmCTIaid3)
#x11(width=11, height=8)
par(mfrow=c(2,2))
plot(lmCTIaid3) # Residuals are slightly better, not normal

# month term is barely significant
# IMOS and Ajani are significantly different from D&C but not from each other
Per3 <- visreg(lmCTIaid3, "Period", rug = FALSE, gg = TRUE) + theme_bw(base_size = 14)+
  ylab(bquote("Community Temperature Index ("* degree *"C)")) +
  xlab("Period")
mon3 <- visreg(lmCTIaid3, "Monh", rug = FALSE, gg = TRUE) + theme_bw(base_size = 14)+
  ylab(bquote("Community Temperature Index ("* degree *"C)")) +
  xlab("Month") + scale_x_continuous(breaks = seq(0.52,6.28,length.out = 12), label = c("J","F","M","A","M","J","J","A","S","O","N","D")) 

#x11(width = 12, height = 5)
p3 <- Per3 + mon3
p3
ggsave("PHB_CTI_1phyto.png", p3, dpi = 1200)

################################################
# Case 4 - use species that occur in IMOS

taxa4 <- aidSTI %>% filter(Survey == 'IMOS') %>% distinct(Taxon) %>% droplevels() 
aidSTI4 <- aidSTI %>% filter(Taxon %in% levels(taxa4$Taxon)) %>% droplevels() 
noSTI4 <- distinct(subset(aidSTI4[,c(3,8)], is.na(aidSTI4$STI))) # check what has no STI , should be higher taxonomic groups or rare species

aidcti4 <- aidSTI4 %>% filter(!is.na(STI))  %>% 
  group_by(Survey, SampleDate) %>% summarise(cti = sum(AbundanceM3*STI, na.rm=TRUE)/sum(AbundanceM3, na.rm=TRUE),
                                             n = n()) %>%
  #filter(n>2) %>% # remove samples where CTI only made up of 2 or less species
  mutate(Dates = ymd_hms(SampleDate),
         Month = month(Dates),
         Monh = Month/12*2*base::pi,
         Period = ifelse(Survey=='Ajani', "1998-2009", 
                         ifelse(Survey=='IMOS', "2009-2018", '1931-1932')),
         Year = year(Dates)) 

# Use lm to see how CTI changes with years after removing the seasonal cycle

aidcti4$Period <- as.factor(aidcti4$Period)
lmCTIaid4 <- lm(cti ~ Period + Harm(Monh, k=1), data=aidcti4) #k2 makes month wobbly, k=1 insignificant
summary(lmCTIaid4)
anova(lmCTIaid4)
#x11(width=11, height=8)
par(mfrow=c(2,2))
plot(lmCTIaid4) # Residuals are slightly better, not normal

# month term is barely significant
# IMOS and Ajani are significantly different from D&C but not from each other
Per4 <- visreg(lmCTIaid4, "Period", rug = FALSE, gg = TRUE) + theme_bw(base_size = 14)+
  ylab(bquote("Community Temperature Index ("* degree *"C)")) +
  xlab("Period")
mon4 <- visreg(lmCTIaid4, "Monh", rug = FALSE, gg = TRUE) + theme_bw(base_size = 14)+
  ylab(bquote("Community Temperature Index ("* degree *"C)")) +
  xlab("Month") + scale_x_continuous(breaks = seq(0.52,6.28,length.out = 12), label = c("J","F","M","A","M","J","J","A","S","O","N","D")) 

#x11(width = 12, height = 5)
p4 <- Per4 + mon4
p4
ggsave("PHB_CTI_IMOSphyto.png", p3, dpi = 1200)


###########################################################################
## look at relative abundances over periods for IMOS species
addAbun4 <- aidSTI4  %>% filter(!is.na(STI))  %>% 
          group_by(Survey) %>% summarise(sums = sum(AbundanceM3))
pie4 <- inner_join(aidSTI4, addAbun4, by=c("Survey")) %>% filter(!is.na(STI))  %>% 
  mutate(pies = AbundanceM3/sums) %>% group_by(Survey, Taxon) %>% summarise(sums = sum(pies)) %>%
  mutate(Period = ifelse(Survey=='Ajani', "1998-2009", ifelse(Survey=='IMOS', "2009-2018", '1931-1932')))
pie4 <- pie4 %>% mutate(TaxonT = ifelse(sums>0.01, as.character(Taxon), "Other")) %>% 
  group_by(Survey, TaxonT, Period) %>% summarise(tots = sum(sums)) %>% complete(TaxonT, Period) %>%
  mutate(tots = ifelse(is.na(tots), 0, tots)) 

pie4$TaxonT <- as.factor(pie4$TaxonT)

#x11(width = 11, height = 8)
percPlot4 <- ggplot(data=pie4, aes(Period, tots)) + geom_col(aes(fill=TaxonT)) + theme_bw()
percPlot4
ggsave("PHB_prop_IMOSsurveysphyto.png", percPlot4, dpi = 1200)
##biggest change is decrease in A. glacilis and increase in L. danicus

##########################################################
## look at relative abundances over periods for species occuring in all 3 surveys
addAbun <- aidSTIt  %>% filter(!is.na(STI)) %>% 
  group_by(Survey) %>% summarise(sums = sum(AbundanceM3))
pie <- inner_join(aidSTIt, addAbun, by=c("Survey")) %>% filter(!is.na(STI)) %>% 
  mutate(pies = AbundanceM3/sums) %>% group_by(Survey, Taxon) %>% summarise(sums = sum(pies)) %>%
  mutate(Period = ifelse(Survey=='Ajani', "1998-2009", ifelse(Survey=='IMOS', "2009-2018", '1931-1932')))
pie2 <- pie %>% mutate(TaxonT = ifelse(sums>0.01, as.character(Taxon), "Other")) %>% 
  group_by(Survey, TaxonT, Period) %>% summarise(tots = sum(sums)) %>% complete(TaxonT, Period) %>%
  mutate(tots = ifelse(is.na(tots), 0, tots)) 

pie2$TaxonT <- as.factor(pie2$TaxonT)

#x11(width = 11, height = 8)
percPlot <- ggplot(data=pie2, aes(Period, tots)) + geom_col(aes(fill=TaxonT)) + theme_bw()
percPlot
ggsave("PHB_prop_3surveysphyto.png", percPlot, dpi = 1200)
##biggest change is decrease in A. glacilis and increase in L. danicus

############################################################################
## look at relative abundances over periods for species in at least 2 surveys
addAbun <- aidSTI2  %>% filter(!is.na(STI)) %>% filter(!grepl("imbricata", Taxon) &
                                                         !grepl("Thalassionema", Taxon)) %>% 
  group_by(Survey) %>% summarise(sums = sum(AbundanceM3))
pie3 <- inner_join(aidSTI2, addAbun, by=c("Survey")) %>% filter(!is.na(STI)) %>% filter(!grepl("imbricata", Taxon) &
                                                                                          !grepl("Thalassionema", Taxon)) %>% 
  mutate(pies = AbundanceM3/sums) %>% group_by(Survey, Taxon) %>% summarise(sums = sum(pies)) %>%
  mutate(Period = ifelse(Survey=='Ajani', "1998-2009", ifelse(Survey=='IMOS', "2009-2018", '1931-1932')))
pie4 <- pie3 %>% mutate(TaxonT = ifelse(sums>0.01, as.character(Taxon), "Other")) %>% 
  group_by(Survey, TaxonT, Period) %>% summarise(tots = sum(sums)) %>% complete(TaxonT, Period) %>%
  mutate(tots = ifelse(is.na(tots), 0, tots)) 

pie4$TaxonT <- as.factor(pie4$TaxonT)

#x11(width = 11, height = 8)
percPlot2 <- ggplot(data=pie4, aes(Period, tots)) + geom_col(aes(fill=TaxonT)) + theme_bw()
percPlot2
ggsave("PHB_prop_2surveysphyto.png", percPlot2, dpi = 1200)


#########################################################################################
## look at relative abundances over periods for species occuring in at least 1 survey

aidSTI_pies <- aidSTI %>% filter(!is.na(STI))  %>%
  filter(!grepl("Chaetoceros", Taxon) & !grepl("Thalassiosira s", Taxon) &
           !grepl("Thalassionema", Taxon))

addAbun <- aidSTI_pies  %>%
  group_by(Survey) %>% summarise(sums = sum(AbundanceM3))
pie5 <- inner_join(aidSTI_pies, addAbun, by=c("Survey")) %>% filter(!is.na(STI)) %>% 
  filter(!grepl("Chaetoceros", Taxon) & !grepl("Thalassiosira s", Taxon) &
           !grepl("Thalassionema", Taxon)) %>%
  mutate(pies = AbundanceM3/sums) %>% group_by(Survey, Taxon) %>% summarise(sums = sum(pies)) %>%
  mutate(Period = ifelse(Survey=='Ajani', "1998-2009", ifelse(Survey=='IMOS', "2009-2018", '1931-1932')))
pie6 <- pie5 %>% mutate(TaxonT = ifelse(sums>0.01, as.character(Taxon), "Other")) %>% 
  group_by(Survey, TaxonT, Period) %>% summarise(tots = sum(sums)) %>% complete(TaxonT, Period) %>%
  mutate(tots = ifelse(is.na(tots), 0, tots)) 

pie6$TaxonT <- as.factor(pie6$TaxonT)

x11(width = 11, height = 8)
percPlot3 <- ggplot(data=pie6, aes(Period, tots)) + geom_col(aes(fill=TaxonT)) + theme_bw()
percPlot3
ggsave("PHB_prop_1phyto.png", percPlot3, dpi = 1200)

