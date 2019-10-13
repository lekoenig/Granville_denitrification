## The objective of this script is to analyze BIWMA salt marsh denitrification data collected in 2018 (Kayleigh Granville)
## Last updated 10 October 2019 by LE Koenig


## Load packages:
library(dplyr)         # general data cleaning and manipulation
library(ggplot2)       # create plots
library(cowplot)       # formatting plots created using ggplot
library(rnoaa)         # r package for interfacing with NCDC climate datasets (need API token to access noaa data) 
library(nlme)          # linear mixed effects models
library(lsmeans)       # posthoc multiple comparisons on linear mixed effects models
library(forcats)       # package for categorical variables (used to order factor variables)
library(car)           # used to check ANOVA assumptions
library(lubridate)     # deal with dates and times
library(RcppRoll)      # calculate rolling summary stats
library(dataRetrieval) # r package for interfacing with USGS streamflow data

## ------------------------------ Read in raw data ----------------------------- ##
  
  # Load master data file containing denitrification rates and anion covariate data:
  data <- read.csv("./data/Compiled_Final_Data_20191002.csv",header=TRUE,stringsAsFactors = FALSE)
  data$Month <- as.factor(data$Date)
  data$VegZone <- as.factor(substring(data$SampleID,first=1,last=2))
  
  #str(data)
  data$Month <- factor(data$Month, levels = c("May","June","July","August","October"))
  data$VegZone <- factor(data$VegZone, levels = c("LM","HM","PH"))

  # Inspect raw data - there should be samples from 3 different vegetation zones x 5 sampling days (n = 15 unique VegZone*Month samples):
  length(unique(data$SampleID))

  # Load data file containing data on percent vegetative cover (measured in August 2018):
  veg.dat <- read.csv("./data/BIWMA_VegetationSurvey.csv",header=TRUE,stringsAsFactors = FALSE)


## --------------- Bring in NOAA NCDC climate data for BIWMA --------------- ##

  ## NOTE THAT YOU WILL NEED TO REQUEST A UNIQUE NCDC WEB SERVICE KEY TOKEN TO ACCESS THE NOAA API ##
  
  # Get GHCND stations near BIWMA:
  BIWMA <- data.frame(lat = 41.340194, lon = -71.877144)
  stations <- ghcnd_stations() %>%
              filter(state == "CT" | state == "RI") 

  # Find stations within 4 km of Barn Island lat/lon:
  BI.precip <- meteo_distance(stations, BIWMA$lat, BIWMA$lon, units = "deg", radius = 4,limit = NULL) # radius in km
  BI.precip <- BI.precip[which(BI.precip$element == "PRCP"),]

  # Download precip data for stations within 4 km of Barn Island:
  out1 <- ncdc(datasetid='GHCND', stationid='GHCND:US1CTNL0018', datatypeid='PRCP', startdate = '2018-03-01', enddate = '2018-10-31', token=token,limit=600)
  out1b <- ncdc(datasetid='GHCND', stationid='GHCND:US1CTNL0040', datatypeid='PRCP', startdate = '2018-03-01', enddate = '2018-10-31', token=token,limit=600)
  out1c <- ncdc(datasetid='GHCND', stationid='GHCND:US1CTNL0024', datatypeid='PRCP', startdate = '2018-03-01', enddate = '2018-10-31', token=token,limit=600)

  # Download temperature data for station closest to Barn Island (precip stations above don't seem to record temperature):
  out2 <- ncdc(datasetid='GHCND', stationid='GHCND:USW00014794', datatypeid='TMAX', startdate = '2018-03-01', enddate = '2018-10-31', token=token,limit=600)
  out3 <- ncdc(datasetid='GHCND', stationid='GHCND:USW00014794', datatypeid='TMIN', startdate = '2018-03-01', enddate = '2018-10-31', token=token,limit=600)

  # Precip is stored in tenths of an inch - convert units to mm for all three data records:
  out1$data$precip_mm <- out1$data$value/10
  out1$data$precip_in <- out1$data$precip_mm*0.0393701
  out1$data$date <- as.Date(out1$data$date)

  out1b$data$precip_mm <- out1b$data$value/10
  out1b$data$precip_in <- out1b$data$precip_mm*0.0393701
  out1b$data$date <- as.Date(out1b$data$date)

  out1c$data$precip_mm <- out1c$data$value/10
  out1c$data$precip_in <- out1c$data$precip_mm*0.0393701
  out1c$data$date <- as.Date(out1c$data$date)

  # Join precip data records to calculate average precipitation for the local "region" around Barn Island:
  out <- left_join(out1c$data[,c("date","datatype","station","precip_mm")], out1b$data[,c("date","datatype","station","precip_mm")], by='date') %>%
          left_join(., out1$data[,c("date","datatype","station","precip_mm")], by='date') 
  out$avg_precip_mm <- apply(out[,c("precip_mm.x","precip_mm.y","precip_mm")],1,function(x) mean(x,na.rm=T))

  # Calculate 7-day, 10-day, 14-day, 21-day, and 30-day average antecedent rainfall for Barn Island:
  StoningtonCT.precip <- out %>% mutate(roll_2day = RcppRoll::roll_sum(avg_precip_mm,2,align="right",fill=NA),
                                        roll_7day = RcppRoll::roll_sum(avg_precip_mm,7,align="right",fill=NA),
                                        roll_10day = RcppRoll::roll_sum(avg_precip_mm,10,align="right",fill=NA),
                                        roll_14day = RcppRoll::roll_sum(avg_precip_mm,14,align="right",fill=NA),
                                        roll_21day = RcppRoll::roll_sum(avg_precip_mm,21,align="right",fill=NA),
                                        roll_30day = RcppRoll::roll_sum(avg_precip_mm,30,align="right",fill=NA))
  
  # Convert temperature units to degrees celsius for both the max and min daily temperature records:
  out2$data$MaxTemp_C <- out2$data$value/10
  out2$data$MaxTemp_F <- (out2$data$MaxTemp_C*9/5)+32
  out2$data$date <- as.Date(out2$data$date)
  
  out3$data$MinTemp_C <- out3$data$value/10
  out3$data$MinTemp_F <- (out3$data$MinTemp_C*9/5)+32
  out3$data$date <- as.Date(out3$data$date)

  # Join climate data:
  Table1 <- data.frame(Month.N = c("May","June","July","August","October"),
                     Date = c(as.Date("05/11/2018","%m/%d/%Y"),
                              as.Date("06/18/2018","%m/%d/%Y"),
                              as.Date("07/24/2018","%m/%d/%Y"),
                              as.Date("08/13/2018","%m/%d/%Y"),
                              as.Date("10/09/2018","%m/%d/%Y")))
  Table1 <- left_join(Table1,StoningtonCT.precip[,c("date","datatype.x","avg_precip_mm","roll_2day","roll_7day","roll_10day","roll_14day","roll_21day","roll_30day")],by=c("Date"="date"))
  names(Table1)[names(Table1) == 'datatype.x'] <- 'datatype'

  Table1 <- left_join(Table1,out3$data[,c("date","datatype","station","MinTemp_C","MinTemp_F")],by=c("Date"="date"))
  Table1 <- left_join(Table1,out2$data[,c("date","datatype","MaxTemp_C","MaxTemp_F")],by=c("Date"="date"))
  Table1$Month.N <- factor(Table1$Month.N,levels = c("May","June","July","August","October"))

  # Download tidal height data from NOAA Tides and Currents station 8461490: New London, Thames River, CT:
  tide.dat <- coops_search(begin_date = 20180401, end_date = 20181031, station_name = 8461490,
               datum = "MLLW", units = "metric", 
               #product = "monthly_mean",
               product = "hourly_height",
               time_zone = "lst",
               application = "rnoaa")
  tide.metadata <- do.call("rbind",tide.dat$metadata)
  tide.dat <- tide.dat$data
  tide.dat$date <- as.Date(tide.dat$t)
  tide.dat$month <- month(tide.dat$date)
  
  # Convert tidal data to daily means/max:
  tide.dat <- tide.dat %>% group_by(date) %>% summarize(meantide = mean(v,na.rm=T),
                                                     maxtide = max(v,na.rm=T))
  
  # Calculate 7-day mean tide:
  tide.dat <- tide.dat %>% mutate(meantide = RcppRoll::roll_mean(meantide,1,align="right",fill=NA),
                                  meantide_7day = RcppRoll::roll_mean(meantide,7,align="right",fill=NA),
                                  meantide_14day = RcppRoll::roll_mean(meantide,14,align="right",fill=NA),
                                  maxtide_7day = RcppRoll::roll_max(maxtide,7,align="right",fill=NA),
                                  maxtide_14day = RcppRoll::roll_max(maxtide,14,align="right",fill=NA))
  
  # Join tidal data with other environmental data:
  Table1 <- left_join(Table1,tide.dat[,c("date","meantide","maxtide","meantide_7day","meantide_14day","maxtide_7day","maxtide_14day")],by=c("Date"="date"))

  # Save climate data:
  write.csv(Table1,"./output/Table1_EnvironmentalConditions.csv",row.names = FALSE)
  saveRDS(Table1,"./output/Table1_EnvironmentalConditions.rds")
  
  
  # Download USGS streamflow data for site 01118500, PAWCATUCK RIVER AT WESTERLY, RI:
  siteNumber  <-'01118500'
  siteINFO <- readNWISsite(siteNumber)
  parameterCd <-c('00060','00065')
  startDate   <-"2018-05-01"
  endDate     <-"2018-10-31"
  flow.dat   <- readNWISuv(siteNumber, parameterCd,startDate, endDate,tz="America/New_York") # USGS data stored in Eastern U.S. local time
  flow.dat   <- renameNWISColumns(flow.dat)
  flow.dat$date <- as.Date(flow.dat$dateTime)
  
  # convert flow data from cfs to L/s:
  flow.dat$Flow_m3s <- flow.dat$Flow_Inst*0.0283168
  
  # calculate mean daily flow:
  flow.dat.daily <- flow.dat %>% group_by(date) %>% summarize(totalvolume_m3 = (sum(Flow_m3s)*15*60))
  
  # convert flow (m3/s) to runoff (m/d)
  flow.dat.daily$runoff_mmd <- (flow.dat.daily$totalvolume_m3/(siteINFO$drain_area_va*2.58999*1000*1000))*1000
    

 
  
## ------------------------------ Plot the vegetation cover by vegetation zone ----------------------------- ##
  
  # Rename veg.dat columns:
  for(i in 2:length(veg.dat)){
    colnames(veg.dat)[i] <- substring(colnames(veg.dat[i]), first=9)
  }
  
  veg.dat$VegZone <- as.factor(substring(veg.dat$Site,first=1,last=2))
  
  # Keep only the 6 most dominant vegetation types:
  veg.dat$other <- apply(veg.dat[,c("Salicornia","Juncus","Baccharius","SeaLavender","Sueda","Iva","Arrowleaf")],1,sum)
  veg.dat <- veg.dat[,c("Site","Bare","Litter","Alter","Distichlis","Patens","Phrag","other","VegZone")]
  
  # Calculate mean percent composition across sub-plots:
  veg.dat.summarize <- veg.dat %>% group_by(VegZone) %>%
    summarize(Bare = mean(Bare),
              Litter = mean(Litter),
              Alter = mean(Alter),
              Distichlis = mean(Distichlis),
              Patens = mean(Patens),
              Phrag = mean(Phrag),
              other = mean(other))
  
  # Convert data format from wide to long (with Species as a factor):
  veg.dat.long <- veg.dat %>% tidyr::gather(Species,Percent_cover, Bare:other)
  veg.dat.long$VegZone <- factor(veg.dat.long$VegZone, levels = c("LM","HM","PH"))
  veg.dat.long$Site <- factor(veg.dat.long$Site, levels = c("LM1","LM2","LM3","LM4","LM5",
                                                            "HM1","HM2","HM3","HM4","HM5",
                                                            "PH1","PH2","PH3","PH4","PH5"))
  veg.sum.dat.long <- veg.dat.summarize %>% tidyr::gather(Species,Percent_cover, Bare:other)
  veg.sum.dat.long$VegZone <- factor(veg.sum.dat.long$VegZone, levels = c("LM","HM","PH"))
  
  # Create stacked bar charts of vegetation composition across all of the sub-plots:
  my.palette <- c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5")
  
  plot_grid(
    veg.sum.dat.long %>% ggplot() + geom_bar(aes(y = Percent_cover, x = VegZone, fill = Species),stat="identity")+
      scale_fill_manual(values = my.palette)+
      labs(y=expression(Average~percent~cover~across~the~"5"~plots))+
      NULL,
    
    veg.dat.long %>% ggplot() + geom_bar(aes(y = Percent_cover, x = Site, fill = Species),stat="identity")+
      scale_fill_manual(values = my.palette)+
      labs(y=expression(Percent~cover))+
      theme(legend.position="none"),
    ncol=2,align="v")

    
## ------------------------------ Potential denitrification rates versus N2O Production rates ----------------------------- ##
  
  # Plot the overall relationship (across month x vegetation zone combinations):
  data %>% ggplot() + geom_point(aes(x=DenitRate,y=N2ORate,color=VegZone),size=1.3)  +
    scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3")) +
    geom_smooth(aes(x=DenitRate,y=N2ORate,color=VegZone),method="lm",se=FALSE)
  
  # Overall relationship (across month x vegetation zone combinations):
  Nreg.all <- lm(log10(data$N2ORate+1)~log10(data$DenitRate))
  summary(Nreg.all)
  
  # Regressions between potential denitrification and N2O production rates by vegetation zone:
  Nreg.byveg = data %>% 
    group_by(VegZone) %>%
    do(fit = lm(log10(N2ORate+1)~log10(DenitRate),data= .))
  
  Nreg.byveg.coef = Nreg.byveg %>% tidy(fit)
  print(Nreg.byveg.coef)
  
  
## ------------------------------ How do N cycling rates vary seasonally and across vegetation zones in the BIWMA? ----------------------------- ##
  
## 1. POTENTIAL DENITRIFICATION
  
  denit.data <- data[which(!is.na(data$DenitRate)),]
  
  # Plot potential denitrification rates, including possible main (Month, VegZone), random (SampleID), and interactive effects. 
  # Note y-axis is presented on a log-scale.  
  plot_grid(
    data %>% ggplot() + geom_boxplot(aes(x=Month,y=DenitRate))+scale_y_log10()+theme_bw(),
    data %>% ggplot() + geom_boxplot(aes(x=VegZone,y=DenitRate))+scale_y_log10()+theme_bw(),
    data %>% ggplot() + geom_boxplot(aes(x=SampleID,y=DenitRate))+scale_y_log10()+theme_bw(),
    # interaction plot:
    ggplot(denit.data, aes(x = Month, y = log10(DenitRate), color = VegZone)) + geom_boxplot(size=0.5) +
      stat_summary(fun.y = mean, geom = "line", aes(group = VegZone), size = 1.2) + theme_bw() +
      scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3")),
    ncol=2,nrow=2,align="h")
  
  # Check - is the experimental design balanced?
  table(denit.data$VegZone, denit.data$Month)
  # October-LM is missing one sample, but aside from data loss in August, the sample size is uniform across treatment levels.  

  # Two-way (factorial) ANOVA with SampleID (i.e., sub-plot) as a random effect, assuming that samples from the same sampleID over time are correlated:
  denit=lme(log10(DenitRate) ~ Month*VegZone, data=denit.data, random=~1|SampleID, method="ML")
  anova(denit) 
  
  # Pairwise comparisons by month: 
  #emmeans(denit,pairwise ~ Month)
  posthoc.denit=lsmeans(denit, pairwise~Month, adjust="tukey",conf=.95)   
  summary(posthoc.denit)
  
  # Check model assumptions:
  par(mfrow=c(2,2),mar=c(3.0,1.8,1.8,1.0),oma=c(2,1.5,1.5,1.5))
  res.denit <- residuals(denit)
  qqnorm(res.denit)
  qqline(res.denit)
  xfit <- seq(min(res.denit), max(res.denit), length = 40) 
  yfit <- dnorm(xfit, mean = mean(res.denit), sd = sd(res.denit)) 
  hist(res.denit,freq=FALSE)
  lines(xfit, yfit, col = "darkblue", lwd = 1.8)
  plot(fitted(denit),residuals(denit))
  shapiro.test(res.denit)                                         # test assumption of normality
  car::leveneTest(log10(DenitRate) ~ Month, data = denit.data)    # test assumption of homogeneity of variance

    
## 2. N2O PRODUCTION
  
  # Plot potential N2O Production rates, including possible main (Month, VegZone), random (SampleID), and interactive effects. 
  # Note y-axis is presented on a log-scale.  
  plot_grid(
    data %>% ggplot() + geom_boxplot(aes(x=Month,y=N2ORate))+scale_y_log10()+theme_bw(),
    data %>% ggplot() + geom_boxplot(aes(x=VegZone,y=N2ORate))+scale_y_log10()+theme_bw(),
    data %>% ggplot() + geom_boxplot(aes(x=SampleID,y=N2ORate))+scale_y_log10()+theme_bw(),
    # interaction plot:
    ggplot(denit.data, aes(x = Month, y = log10(N2ORate+1), color = VegZone)) + geom_boxplot(size=0.5) +
      stat_summary(fun.y = mean, geom = "line", aes(group = VegZone), size = 1.2) + theme_bw() +
      scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3")),
    ncol=2,nrow=2,align="h")
  
  # Two-way (factorial) ANOVA with SampleID (i.e., sub-plot) as a random effect, assuming that samples from the same sampleID over time are correlated:  
  n2o=lme(log10(N2ORate+1) ~ Month*VegZone, data=data, random=~1|SampleID, method="ML")
  anova(n2o) 
  
  # Pairwise comparisons by month*vegetation zone: 
  #emmeans(n2o,pairwise ~ Month*VegZone)
  posthoc.n2o=lsmeans(n2o, pairwise~Month*VegZone, adjust="tukey",conf=.95)   
  summary(posthoc.n2o$contrasts)[summary(posthoc.n2o$contrasts)$p.value < 0.05,]
  
  # Create vector of unique interaction combinations of Month x VegZone:
  data$unique.trt <- interaction(data$Month,data$VegZone)
  
  # Check model assumptions:
  par(mfrow=c(2,2),mar=c(3.0,1.8,1.8,1.0),oma=c(2,1.5,1.5,1.5))
  res.n2o <- residuals(n2o)
  qqnorm(res.n2o)
  qqline(res.n2o)
  xfit <- seq(min(res.n2o), max(res.n2o), length = 40) 
  yfit <- dnorm(xfit, mean = mean(res.n2o), sd = sd(res.n2o)) 
  hist(res.n2o,freq=FALSE)
  lines(xfit, yfit, col = "darkblue", lwd = 1.8)
  plot(fitted(n2o),residuals(n2o))
  shapiro.test(res.n2o)                                             # test assumption of normality
  car::leveneTest(log10(N2ORate+1) ~ unique.trt, data = data)    # test assumption of homogeneity of variance
  

## ------------------------------ Climate and other covariate data ----------------------------- ##
  
  climate.dat <- readRDS("./output/Table1_EnvironmentalConditions.rds")
  climate.dat.long <- climate.dat %>% tidyr::gather(Antcdt_rainfall,value_mm, roll_7day:roll_30day)

  # Plot climate and covariate trends across months:
  plot_grid(
    climate.dat.long %>% ggplot() + 
      geom_boxplot(aes(x=Month.N,y=value_mm))+
      geom_point(aes(x=Month.N,y=value_mm,color=Antcdt_rainfall),size=2)+
      theme_bw() + 
      labs(x="Month",y=expression(Antecedent~precip~(mm))),
    
    data %>% ggplot() + geom_boxplot(aes(x=Month,y=LOI))+geom_point(aes(x=Month,y=LOI,color=VegZone))+theme_bw(),
    data %>% ggplot() + geom_boxplot(aes(x=Month,y=SO4_normalized))+geom_point(aes(x=Month,y=SO4_normalized,color=VegZone))+theme_bw()+theme(legend.position="none"),
    data %>% ggplot() + geom_boxplot(aes(x=Month,y=NH4_normalized))+geom_point(aes(x=Month,y=NH4_normalized,color=VegZone))+theme_bw()+theme(legend.position="none"),
    data %>% ggplot() + geom_boxplot(aes(x=Month,y=SoilMoistureFraction))+geom_point(aes(x=Month,y=SoilMoistureFraction,color=VegZone))+theme_bw()+theme(legend.position="none"),
    ncol=2
    
  )
  

## ------------------------------ Scatterplots between environmental/soil covariates and N cycling rates ----------------------------- ##
  
## 1. POTENTIAL DENITRIFICATION
  plot_grid(
    data %>% ggplot() + geom_point(aes(x=SO4_normalized,y=DenitRate))+theme_bw(),
    data %>% ggplot() + geom_point(aes(x=NH4_normalized,y=DenitRate))+theme_bw(),
    data %>% ggplot() + geom_point(aes(x=LOI,y=DenitRate))+theme_bw(),
    data %>% ggplot() + geom_point(aes(x=SoilMoistureFraction,y=DenitRate))+theme_bw(),
    ncol=2)
  
## 2. N2O PRODUCTION
  plot_grid(
    data %>% ggplot() + geom_point(aes(x=SO4_normalized,y=N2ORate))+theme_bw(),
    data %>% ggplot() + geom_point(aes(x=NH4_normalized,y=N2ORate))+theme_bw(),
    data %>% ggplot() + geom_point(aes(x=LOI,y=N2ORate))+theme_bw(),
    data %>% ggplot() + geom_point(aes(x=SoilMoistureFraction,y=N2ORate))+theme_bw(),
    ncol=2)
  
## 3. N2O YIELD
  plot_grid(
    data %>% ggplot() + geom_point(aes(x=SO4_normalized,y=N2OYield))+theme_bw(),
    data %>% ggplot() + geom_point(aes(x=NH4_normalized,y=N2OYield))+theme_bw(),
    data %>% ggplot() + geom_point(aes(x=LOI,y=N2OYield))+theme_bw(),
    data %>% ggplot() + geom_point(aes(x=SoilMoistureFraction,y=N2OYield))+theme_bw(),
    ncol=2)
  
  