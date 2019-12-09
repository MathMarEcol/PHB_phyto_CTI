# Calculate STI and NRS CTI for Phytoplankton
# CTI time series graphs at end of page
# Claire: march 2019

## THis version makes assumptions only valid for the Ajani / DC / IMOS comparison work

library(ggplot2)
library(dplyr)

####################################
# for phyto
psti <- read.csv("DataForSTIphyto.csv", header= TRUE, na.strings="(null)") %>% select(c(1,3:7)) %>%
  subset(!is.na(SST) & ABUNDANCE >0) %>% 
  mutate(sst = round(SST/0.5) * 0.5) %>%
  filter(!grepl("/", TAXON)) %>%
  mutate(TAXON = recode(TAXON, "Guinardia flaccida <150 µm length" = "Guinardia flaccida",
                        "Guinardia flaccida >150 µm length" = "Guinardia flaccida",
                        "Guinardia striata (with Richelia)" = "Guinardia striata",
                        "Lauderia annulata < 70 µm length" = "Lauderia annulata",
                        "Lauderia annulata > 70 µm length" = "Lauderia annulata",
                        "Ditylum brightwellii < 40 µm width" = "Ditylum brightwellii",
                        "Dactyliosolen fragillissimus >150 µm length" = "Dactyliosolen fragillissimus",
                        "Dactyliosolen fragillissimus <150 µm length" = "Dactyliosolen fragillissimus",  
                        "Leptocylindrus mediterraneus (no flagellates)" = "Leptocylindrus mediterraneus", 
                        "Leptocylindrus mediterraneus (with flagellates)" = "Leptocylindrus mediterraneus",
                        "Proboscia alata ~ 50 µm cell width" = "Proboscia alata" ,
                        "Asterionellopsis spp. ~50 µm length" = "Asterionellopsis glacialis",
                        "Chaetoceros peruvianus < 40 µm cell width" = "Chaetoceros peruvianus",
                        "Chaetoceros peruvianus > 40 µm cell width" = "Chaetoceros peruvianus",
                        "Nitzschia cf. bicapitata" = "Nitzschia bicapitata"))  %>% droplevels()
  
nrsp <- subset(psti, PROJECT == "nrs" )
mean_nrsp <- mean(log10(nrsp$ABUNDANCE))
nrsz <- subset(psti, PROJECT == "nrsz" )
mean_nrsz <- mean((nrsz$ABUNDANCE))
cprp <- subset(psti, PROJECT == "cpr")
mean_cprp <- mean(cprp$ABUNDANCE)

#means are so different so log nrs data as the abundance scale is so wide

psti <- psti %>% mutate(mean = ifelse(PROJECT=='cpr', mean_cprp, 
                                      ifelse(PROJECT=='nrs', mean_nrsp, mean_nrsz))) %>%
  mutate(mean_abun =ifelse(PROJECT == "nrs", log10(ABUNDANCE+1)/mean, ABUNDANCE/mean))

#rel_mean <- mean(psti$mean_abun)
psti$relab <- psti$mean_abun #/rel_mean # this is pointless as I am not including the assumption that Wayne made yet
psti$absst <- psti$relab*psti$SST

p <- psti[,c(6,4,10)] #SST, TAXON_NAME, relab

dfp <- p %>% group_by(SST, TAXON) %>%
  summarize(relab = sum(relab), freq = n(), a = sum(relab)/n()) 

dfpn <- dfp %>% group_by(TAXON) %>% summarize(ntax = n())
dfp <- inner_join(dfp, dfpn, by="TAXON")
notaxp <- 5
dfp <- subset(dfp, ntax>notaxp)

dfp$TAXON <- factor(dfp$TAXON)
psti_id <- matrix(0, nrow=nlevels(dfp$TAXON), ncol=2)
colnames(psti_id)=c("Taxon", "STI")

## STI via kernel density

kernStep <- 0.1
kernMin <- min(dfp$SST) - 3
kernMax <- max(dfp$SST) + 3
kernN <- round((kernMax - kernMin) / kernStep + 1)
kernTemps <- seq(kernMin, kernMax, length.out=kernN)
kernBw <- 2

kern_yp <- matrix(0, nrow = kernN, ncol = nlevels(dfp$TAXON))
kypout <- matrix(0, nrow = kernN, ncol = nlevels(dfp$TAXON))

for (i in 1:nlevels(dfp$TAXON)) {
  taxon <- levels(dfp$TAXON)[i]
  kernData <- subset(dfp[,1:3], TAXON == taxon)
  kernData$TAXON <- factor(kernData$TAXON)
  kernData$weight <- with(kernData, abs(relab) / sum(relab))
  kernOut <- with(kernData,
                  density(SST, weight=weight,
                          bw=kernBw,
                          from=kernMin,
                          to=kernMax,
                          n=kernN))
  
  
  z <- cbind(kernTemps, kernOut$y) 
  z <- z[which.max(z[,2]),]
  psti_id[i,1] <- paste0(taxon)
  psti_id[i,2] <- as.numeric(z[1])  
  
  kern_yp[,i] <- kernOut$y/sum(kernOut$y) * 100 * mean(kernData$relab)
  kypout[,i] <- kernOut$y
}

psti_id <- as.data.frame(psti_id)
colnames(psti_id) <- c("Taxon", "STI")
psti_id$STI <- as.numeric(as.character(psti_id$STI))

write.csv(psti_id, "pSTI.csv", row.names = FALSE)


#####################################################################
psti <- subset(psti, TAXON != 'Thalassiosira spp. 10-20 µm' | PROJECT != 'cpr' ) #just for this case to make Penny's thalassisira equivalent
psti <- subset(psti, TAXON != 'Thalassiosira spp. < 10 µm' | PROJECT != 'cpr'  ) #just for this case to make Penny's thalassisira equivalent
psti <- subset(psti, TAXON != 'Thalassiosira spp. 20-40 µm' | PROJECT != 'cpr'  ) #just for this case to make Penny's thalassisira equivalent
psti <- subset(psti, TAXON != 'Thalassiosira spp. 40-60 µm' | PROJECT != 'cpr'  ) #just for this case to make Penny's thalassisira equivalent
psti <- subset(psti, TAXON != 'Pseudo-nitzschia delicatisima complex <=3 µm' | PROJECT != 'cpr'  ) #just for this case to make Penny's thalassisira equivalent
psti <- subset(psti, TAXON != 'Pseudo-nitzschia seriata complex >3 µm' | PROJECT != 'cpr'  ) #just for this case to make Penny's thalassisira equivalent


