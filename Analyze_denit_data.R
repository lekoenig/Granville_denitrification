## The objective of this script is to analyze BIWMA salt marsh denitrification data collected in 2018 by Kayleigh Granville
## Last updated 17 September 2019 by LE Koenig


## Load packages:
library(dplyr)     # general data cleaning and manipulation
library(ggplot)    # create plots
library(cowplot)   # formatting plots created using ggplot
library(rnoaa)     # r package for interfacing with NCDC climate datasets (need API token to access noaa data) 
library(nlme)      # linear mixed effects models
library(lsmeans)   # posthoc multiple comparisons on linear mixed effects models


## ------------------------------ Read in raw data ----------------------------- ##

data <- read.csv("./data/LISS.KG.MasterMatrix.7.24.csv",header=TRUE,stringsAsFactors = FALSE)
str(data)
data$Month.N <- as.factor(data$Month.N)
data$VegZone <- as.factor(data$VegZone)
data$Site <- as.factor(data$Site)

cols <- c("Month.N","VegZone","Site")
data[cols] <- sapply(data[cols],as.factor)

veg.dat <- read.csv("./data/BIWMA_VegetationSurvey.csv",header=TRUE,stringsAsFactors = FALSE)



## ------------------------------ Inspect raw data ----------------------------- ##

# How many sites? Should be samples from 3 different vegetation zones on each of 5 days (n=15):
length(unique(data$Site))


## --------------- Repeated measures ANOVA on covariate datasets --------------- ##

# 1. denitrification rate 
  denit.data <- data[which(!is.na(data$Denit.Rate)),]
  denit=lme(Denit.Rate ~ Month.N*VegZone, data=denit.data, random=~1|Site, method="ML")
  anova(denit) 
  ## denitrification rate varies by month, but not by vegetation zone
  
  ## pairwise comparisons by month (averaged over vegetation zone, since veg. zone was not significant):
  posthocsm=lsmeans(denit, pairwise~Month.N, adjust="tukey")   
  summary(posthocsm)
  emmeans(denit,pairwise ~ Month.N)
  emmip(denit, VegZone ~ Month.N)
  ## denitrification rate differs between May-June; May-July; May-October (denitrification rate was significantly greater in May compared to June, July, and October)
  
  
  emmeans(denit, pairwise ~ Month.N)
  
  res.sm<-residuals(soilmoist)
  qqnorm(res.sm)
  hist(res.sm)
  shapiro.test(res.sm)


## --------------- Repeated measures ANOVA on covariate datasets --------------- ##

# 1. soil moisture

  soilmoist=lme(Soil.Moist ~ Month.N*VegZone, data=data, random=~1|Site, method="ML")
  anova(soilmoist)
  posthocsm=lsmeans(soilmoist, pairwise~Month.N, adjust="tukey")   # pairwise comparisons grouped by vegetation zone, since there was no significant relationship between vegetation zone and denitrification
  summary(posthocsm)
  
  res.sm<-residuals(soilmoist)
  qqnorm(res.sm)
  hist(res.sm)
  shapiro.test(res.sm)
  

# 2. ammonium

  nh4=lme(nh4 ~ Month.N*VegZone, data=data, random=~1|Site, method="ML")
  anova(nh4)
  posthoc.nh4=lsmeans(nh4, pairwise~Month.N, adjust="tukey")        
  summary(posthoc.nh4)
  summary(posthoc.nh4$contrasts)[summary(posthoc.nh4$contrasts)$p.value < 0.05,]
  res.nh4<-residuals(nh4)
  qqnorm(res.nh4)
  hist(res.nh4)
  shapiro.test(res.nh4)
  
  





# regressions by group:
  
  byveg = try3 %>% 
          group_by("Month.N") %>%
          do(fit = lm(Denit.Rate ~ Percent_Bare,data= .))




  emmip(noise.lm, type ~ size | side)



## --------------- Bring in NOAA NCDC climate data for BIWMA --------------- ##

# get GHCND stations:

BIWMA <- data.frame(lat = 41.340194, lon = -71.877144)
stations <- ghcnd_stations() %>%
            filter(state == "CT") 

BI.precip <- meteo_distance(stations, BIWMA$lat, BIWMA$lon, units = "deg", radius = 3,limit = NULL) # radius in km



