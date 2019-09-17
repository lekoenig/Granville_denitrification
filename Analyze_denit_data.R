## The objective of this script is to analyze BIWMA salt marsh denitrification data collected in 2018 by Kayleigh Granville
## Last updated 17 September 2019 by LE Koenig


## Load packages:
library(dplyr)     # general data cleaning and manipulation
library(ggplot)    # create plots
library(cowplot)   # formatting plots created using ggplot
library(rnoaa)     # r package for interfacing with NCDC climate datasets (need API token to access noaa data) 
library(nlme)      # linear mixed effects models
library(lsmeans)   # posthoc multiple comparisons on linear mixed effects models

## --------------- Read in raw data --------------- ##

data <- read.csv("./data/LISS.KG.MasterMatrix.7.24.csv",header=TRUE,stringsAsFactors = FALSE)
veg.dat <- read.csv("./data/BIWMA_VegetationSurvey.csv",header=TRUE,stringsAsFactors = FALSE)



## --------------- Inspect raw data --------------- ##

# How many sites? Should be samples from 3 different vegetation zones on each of 5 days (n=15):
length(unique(data$Site))


## --------------- Repeated measures ANOVA --------------- ##

# 1. soil moisture

  soilmoist=lme(Soil.Moist ~ Month.N*VegZone, data=data, random=~1|Site, method="ML")
  anova(soilmoist)
  posthocsm=lsmeans(soilmoist, pairwise~Month.N, adjust="tukey")
  summary(posthocsm)
  
  res.sm<-residuals(soilmoist)
  qqnorm(res.sm)
  hist(res.sm)
  shapiro.test(res.sm)
  

# 2. ammonium

  nh4=lme(nh4 ~ Month.N*VegZone, data=Dataframe, random=~1|Site, method="ML")
  anova(nh4)
  posthoc.nh4=lsmeans(nh4, pairwise~Month.N, adjust="tukey")        
  summary(posthoc.nh4)
  res.nh4<-residuals(nh4)
  qqnorm(res.nh4)
  hist(res.nh4)
  shapiro.test(res.nh4)
  
  
  log.nh4=lme(log.nh4 ~ Month.N*VegZone, data=Dataframe, random=~1|Site, method="ML")
  anova(log.nh4)
  ph.nh4.log=lsmeans(log.nh4, pairwise~Month.N, adjust="tukey")
  summary(ph.nh4.log)
  res.ph.nh4.log<-residuals(log.nh4)
  qqnorm(res.ph.nh4.log)
  hist(res.ph.nh4.log)
  shapiro.test(res.ph.nh4.log)


















## --------------- Bring in NOAA NCDC climate data for BIWMA --------------- ##

# get GHCND stations:

BIWMA <- data.frame(lat = 41.340194, lon = -71.877144)
stations <- ghcnd_stations() %>%
            filter(state == "CT") 

BI.precip <- meteo_distance(stations, BIWMA$lat, BIWMA$lon, units = "deg", radius = 3,limit = NULL) # radius in km



