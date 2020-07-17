## The objective of this script is to analyze BIWMA salt marsh denitrification data collected in 2018 (Kayleigh Granville)
## Last updated 3 April 2020 by LE Koenig


## Load packages:
library(dplyr)         # general data cleaning and manipulation
library(ggplot2)       # create plots
library(ggpubr)        # format some plots created with ggplot
library(cowplot)       # formatting plots created using ggplot (and creating multi-panel figures)
library(rnoaa)         # r package for interfacing with NCDC climate datasets (need API token to access noaa data) 
library(lme4)          # linear mixed effects models
library(lmerTest)      # linear mixed effects models
library(lsmeans)       # posthoc multiple comparisons on linear mixed effects models
library(forcats)       # package for categorical variables (used to order factor variables)
library(car)           # used to check ANOVA assumptions
library(lubridate)     # deal with dates and times
library(RcppRoll)      # calculate rolling summary stats
library(dataRetrieval) # r package for interfacing with USGS streamflow data
library(broom)         # tidy data and run linear models by factor groupings
library(piecewiseSEM)  # estimate marginal and conditional r2 for linear mixed effects models
library(patchwork)     # create multi-panel figures


## NOTE THAT YOU WILL NEED TO REQUEST A UNIQUE NCDC WEB SERVICE KEY TOKEN TO ACCESS THE NOAA API (copy in place of X below) ##
token <- "X"


## ============================================================================= ## 
##                     Read in raw data and format column names                  ##
## ============================================================================= ## 

  # Load master data file containing denitrification rates and anion covariate data:
  data <- read.csv("./data/Compiled_Final_Data_20200128.csv",header=TRUE,stringsAsFactors = FALSE)
  
  # Add sampling dates column to the data (these are the known sampling dates that correspond with the months data that is recorded in the raw data):
  data <- data %>% mutate(date = case_when(
    Date == "May" ~ "5/11/2018",
    Date == "June" ~ "6/18/2018",
    Date == "July" ~ "7/24/2018",
    Date == "August" ~ "8/13/2018",
    Date == "October" ~ "10/9/2018"))
  data$date <- mdy(data$date)
  
  # Additional formatting
  # Create a factor variable that represents the salt marsh zone (LM, HM, PH):
  data$VegZone <- as.factor(substring(data$Sample.ID,first=1,last=2))
  
  # Adjust factor levels for month and salt marsh zone:
  data$Date <- factor(data$Date, levels = c("May","June","July","August","October"))
  data$VegZone <- factor(data$VegZone, levels = c("LM","HM","PH"))
  data$VegZone2 <- factor(data$VegZone,levels=c("LM","HM","PH"),labels=c("Low Marsh","High Marsh", "Phragmites"))
  
  # Simplify N cycling column names for analysis:
  names(data)[names(data) == 'Denit..ng.N...hr...g.dry.soil.'] <- 'DenitRate'
  names(data)[names(data) == 'N2O..ng...hr...g.dry.soil.'] <- 'N2ORate'
  names(data)[names(data) == 'LOI..g.ashed...g.dry.'] <- 'LOI'
  names(data)[names(data) == 'SO4.mg.per.g.wet.soil'] <- 'SO4_wet_normalized'
  names(data)[names(data) == 'NH4.mg.N.per.g.wet.soil'] <- 'NH4_wet_normalized'
  
  # Inspect raw data - there should be samples from 3 different vegetation zones x 5 sampling days (n = 15 unique VegZone*Month samples):
  length(unique(data$Sample.ID))
  
  # Load data file containing data on percent vegetative cover (measured in August 2018):
  veg.dat <- read.csv("./data/BIWMA_VegetationSurvey.csv",header=TRUE,stringsAsFactors = FALSE)



## ============================================================================= ## 
##                     Bring in NOAA NCDC climate data for BIWMA                 ##
## ============================================================================= ## 

## NOTE THAT YOU WILL NEED TO REQUEST A UNIQUE NCDC WEB SERVICE KEY TOKEN TO ACCESS THE NOAA API (replace token <- "X" at top of script) ##

# Download climate data:
  
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
  
  # Precip is stored in tenths of mm - convert units to mm for all three data records:
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
  out <- left_join(out1c$data[,c("date","datatype","station","precip_mm")], 
                   out1b$data[,c("date","datatype","station","precip_mm")], by='date') %>%
         left_join(., out1$data[,c("date","datatype","station","precip_mm")], by='date') 
  out$avg_precip_mm <- apply(out[,c("precip_mm.x","precip_mm.y","precip_mm")],1,function(x) mean(x,na.rm=T))
  
  # Calculate 7-day, 10-day, 14-day, 21-day, and 30-day average antecedent rainfall for Barn Island:
  StoningtonCT.precip <- out %>% mutate(roll_2day = RcppRoll::roll_sum(avg_precip_mm,2,align="right",fill=NA),
                                        roll_7day = RcppRoll::roll_sum(avg_precip_mm,7,align="right",fill=NA),
                                        roll_10day = RcppRoll::roll_sum(avg_precip_mm,10,align="right",fill=NA),
                                        roll_14day = RcppRoll::roll_sum(avg_precip_mm,14,align="right",fill=NA),
                                        roll_21day = RcppRoll::roll_sum(avg_precip_mm,21,align="right",fill=NA),
                                        roll_30day = RcppRoll::roll_sum(avg_precip_mm,30,align="right",fill=NA))
    
  # Temp is stored in tenths of degrees celsius. Convert temperature units to degrees celsius for both the max and min daily temperature records:
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
  Table1 <- left_join(Table1,out2$data[,c("date","datatype","station","MaxTemp_C","MaxTemp_F")],by=c("Date"="date"))
  Table1$Month.N <- factor(Table1$Month.N,levels = c("May","June","July","August","October"))
  
  
# Download local tidal data:

  # Download tidal height data from NOAA Tides and Currents station 8461490: New London, Thames River, CT:
  tide.dat <- coops_search(begin_date = 20180401, end_date = 20181031, station_name = 8461490,
               datum = "MLLW", units = "metric", 
               #product = "monthly_mean",
               product = "hourly_height",
               time_zone = "lst",
               application = "rnoaa")
  tide.metadata <- do.call("rbind",tide.dat$metadata)
  tide.dat <- tide.dat$data
  tide.dat$date <- as.Date(as.character(tide.dat$t))
  tide.dat$month <- month(tide.dat$date)
    
  # Convert tidal data to daily means/max:
  tide.dat2 <- tide.dat %>% group_by(date) %>% summarize(meantide = mean(v,na.rm=T),
                                                        maxtide = max(v,na.rm=T))
  
  # Calculate rolling 7-day mean tide:
  tide.dat3 <- tide.dat2 %>% mutate(meantide = RcppRoll::roll_mean(meantide,1,align="right",fill=NA),
                                  meantide_7day = RcppRoll::roll_mean(meantide,7,align="right",fill=NA),
                                  meantide_14day = RcppRoll::roll_mean(meantide,14,align="right",fill=NA),
                                  maxtide_7day = RcppRoll::roll_max(maxtide,7,align="right",fill=NA),
                                  maxtide_14day = RcppRoll::roll_max(maxtide,14,align="right",fill=NA))

  
# Join tidal data with other environmental data:
  Table1 <- left_join(Table1,tide.dat3[,c("date","meantide","maxtide","meantide_7day","meantide_14day","maxtide_7day","maxtide_14day")],by=c("Date"="date"))
  
# Save climate data:
  write.csv(Table1,"./output/data_processed/Table1_EnvironmentalConditions.csv",row.names = FALSE)
  saveRDS(Table1,"./output/data_processed/Table1_EnvironmentalConditions.rds")
  


## ============================================================================= ## 
##                     Plot vegetation cover by vegetation zone                  ##
## ============================================================================= ## 
  
  # Rename veg.dat columns to remove "Percent" from all of the names:
  for(i in 2:length(veg.dat)){
    colnames(veg.dat)[i] <- substring(colnames(veg.dat[i]), first=9)
  }
  veg.dat$VegZone <- as.factor(substring(veg.dat$Site,first=1,last=2))
  
  # Keep only the 6 most dominant vegetation types from the vegetation survey:
  veg.dat$Other <- apply(veg.dat[,c("Salicornia","Juncus","Baccharius","SeaLavender","Sueda","Iva","Arrowleaf")],1,sum)
  veg.dat <- veg.dat[,c("Site","Bare","Litter","Alter","Distichlis","Patens","Phrag","Other","VegZone")]

  # Calculate total percentages for each vegetation type:
  veg.dat$sum <- apply(veg.dat[,c(2:8)],1,sum)
  
  # Scale each vegetation zone so that all plots sum to 100% coverage:
  for(i in 1:length(veg.dat$Site)){
    veg.dat$Bare_scaled[i] <- round((veg.dat$Bare[i]*100)/veg.dat$sum[i],3)
    veg.dat$Litter_scaled[i] <- round((veg.dat$Litter[i]*100)/veg.dat$sum[i],3)
    veg.dat$Alter_scaled[i] <- round((veg.dat$Alter[i]*100)/veg.dat$sum[i],3)
    veg.dat$Distichlis_scaled[i] <- round((veg.dat$Distichlis[i]*100)/veg.dat$sum[i],3)
    veg.dat$Patens_scaled[i] <- round((veg.dat$Patens[i]*100)/veg.dat$sum[i],3)
    veg.dat$Phrag_scaled[i] <- round((veg.dat$Phrag[i]*100)/veg.dat$sum[i],3)
    veg.dat$Other_scaled[i] <- round((veg.dat$Other[i]*100)/veg.dat$sum[i],3)
  }
  veg.dat$sum_scaled <- apply(veg.dat[,c(11:17)],1,sum)
  
  # Convert data format from wide to long (with Species as a factor):
  veg.dat.long <- veg.dat %>% tidyr::gather(Species,Percent_cover, Bare_scaled:Other_scaled)
  veg.dat.long$VegZone <- factor(veg.dat.long$VegZone, levels = c("LM","HM","PH"))
  veg.dat.long$Site <- factor(veg.dat.long$Site, levels = c("LM1","LM2","LM3","LM4","LM5",
                                                            "HM1","HM2","HM3","HM4","HM5",
                                                            "PH1","PH2","PH3","PH4","PH5"))
  veg.dat.long$Species <- factor(veg.dat.long$Species,levels = c("Bare_scaled","Litter_scaled","Other_scaled","Distichlis_scaled","Alter_scaled","Patens_scaled","Phrag_scaled"),
                                                      labels = c("Bare","Litter","Other","Distichlis","Alter","Patens","Phrag"))
  veg.dat.long$Site2 <- substring(veg.dat.long$Site, first=3)
  veg.dat.long <- veg.dat.long %>% 
                  mutate(VegZone2 = case_when(
                         VegZone == "LM" ~ "Low Marsh",
                         VegZone == "HM" ~ "High Marsh",
                         VegZone == "PH" ~ "Phragmites"))
  veg.dat.long$VegZone2 <- factor(veg.dat.long$VegZone2,levels=c("Low Marsh","High Marsh", "Phragmites"))
  
  # Plot the vegetation composition across each of the sub-plots:
    my.veg.palette <-  c("#999999", "#56B4E9", "#E69F00", "#F0E442","#009E73", "#0072B2", "#D55E00",  "#009E73")
    
    # open jpeg file
    jpeg("./output/figures/Figure2.jpg", width = 8, height = 4.7,units = "in",res = 350)
    
    # create plot
    vegplot.indiv <- veg.dat.long %>% ggplot() + 
      geom_bar(aes(y = Percent_cover, x = Site2, fill = Species),stat="identity")+
      scale_fill_manual(values = my.veg.palette)+    
      labs(y=expression(Percent~cover),x="Plot")+
      facet_grid(. ~ VegZone2)+
      theme_cowplot()+
      theme(legend.position="right",
            axis.title = element_text(size=16),
            axis.text = element_text(size=14),
            legend.text = element_text(size=13),
            legend.key.size = unit(1.5,"line"),
            legend.title = element_text(size=15),
            axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
            #strip.background=element_blank(),          # un-comment this line to get rid of gray boxes around panel titles
            strip.text.x = element_text(size=14,margin = margin(.1, 0, .1, 0, "cm")))
    print(vegplot.indiv) 
    
    # close jpeg file
    dev.off()
  
  # Plot mean vegetation composition across each of the zones (averaged across all of the sub-plots)
    
    # open jpeg file
    jpeg("./output/figures/Figure2_subplot_means.jpg", width = 5.5, height = 4.7,units = "in",res = 350)
    
    # create plot
    vegplot.means <-  veg.dat.long %>% group_by(VegZone, Species) %>% summarize(Mean_cover = mean(Percent_cover,na.rm=T)) %>% 
      ggplot() + geom_bar(aes(y = Mean_cover, x = VegZone, fill = Species),stat="identity")+
      scale_fill_manual(values = my.veg.palette)+
      labs(y=expression(Average~percent~cover~across~the~"5"~plots))+
      labs(y=expression(Mean~percent~cover),x="Salt marsh zone")+
      theme_cowplot()+
      theme(legend.position="right",
            axis.title = element_text(size=16),
            axis.text = element_text(size=14),
            legend.text = element_text(size=13),
            legend.key.size = unit(1.5,"line"),
            legend.title = element_text(size=15),
            axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
            strip.text.x = element_text(size=14,margin = margin(.1, 0, .1, 0, "cm")))
    print(vegplot.means) 
    
    # close jpeg file
    dev.off()
    
    
## ============================================================================= ## 
##                   N CYCLING ACROSS SEASONS AND SALT MARSH ZONES               ##
## ============================================================================= ## 
    
     
## ============================================================================= ## 
##                           1. POTENTIAL DENITRIFICATION                        ##
## ============================================================================= ## 
  
## Plot data
    
  # Plot potential denitrification rates, including possible main (Month, VegZone), random (SampleID), and interactive effects. 
  # Note y-axis is presented on a log-scale.  
  plot_grid(
    data %>% ggplot() + geom_boxplot(aes(x=Date,y=DenitRate))+scale_y_log10()+theme_bw(),
    data %>% ggplot() + geom_boxplot(aes(x=VegZone,y=DenitRate))+scale_y_log10()+theme_bw(),
    data %>% ggplot() + geom_boxplot(aes(x=Sample.ID,y=DenitRate,color=VegZone))+scale_y_log10()+theme_bw()+scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3"))+theme(legend.position="none"),
    # interaction plot:
    ggplot(data, aes(x = Date, y = log10(DenitRate), color = VegZone)) + geom_boxplot(size=0.5) +
      stat_summary(fun = mean, geom = "line", aes(group = VegZone), size = 1.2) + theme_bw() +
      scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3")),
      ncol=2,nrow=2,align="h")
    
## Check - is the experimental design balanced?
  table(data[which(!is.na(data$DenitRate)),]$VegZone, data[which(!is.na(data$DenitRate)),]$Date)

## Analyze data across season and salt marsh zone:
  
  #Two-way ANOVA with SampleID (i.e. sub-plot) as a random effect, assuming that samples from the same plot (sampleID) are correlated over time. 
  #Note that potential denitrification data were natural log-transformed prior to analysis.
  denit = lme4::lmer(log(DenitRate) ~ Date*VegZone + (1|Sample.ID),data=data,REML=TRUE)
  anova(denit)
  # run denit model again using lmerTest Satterthwaite approximation of degrees of freedom in lmerTest:
  options(contrasts = c("contr.sum", "contr.poly"))
  denit <- lmerTest::lmer(log(DenitRate) ~ Date*VegZone + (1|Sample.ID), data = data,REML = TRUE)
  anova(denit)
      
  # Pairwise comparisons by month: 
  posthoc.denit.emm <- emmeans(denit,c("Date"),adjust="tukey",conf=.95)
  pairs(posthoc.denit.emm,type="response")
  
  # Compare effect size:
  summary(posthoc.denit.emm,type="response")

  # Calculate marginal and conditional r2:
  piecewiseSEM::rsquared(denit)
  #r2glmm::r2beta(model=denit,partial=TRUE,method='nsj')
  
  # Calculate relative AIC weights to assess variable importance:
  #require(MuMIn)
  #denit1 = lme4::lmer(log(DenitRate) ~ Date + (1|Sample.ID),data=data,REML=TRUE)
  #denit2 = lme4::lmer(log(DenitRate) ~ Date*VegZone + (1|Sample.ID),data=data,REML=TRUE)
  #denit3 = lme4::lmer(log(DenitRate) ~ VegZone + (1|Sample.ID),data=data,REML=TRUE)
  #denit.aic = c(MuMIn::AICc(denit1),MuMIn::AICc(denit2),MuMIn::AICc(denit3))
  #delAIC <- denit.aic - min(denit.aic)
  #relLik <- exp(-0.5 * delAIC)
  #aicweight <- relLik/sum(relLik)
  #aic.table <- data.frame(AICc = denit.aic, delAIC = delAIC, relLik = relLik, 
  #                        weight = aicweight)
  #rownames(aic.table) <- c("date", "date.vegzone", "vegzone")

## Check model assumptions:
  par(mfrow=c(2,2),mar=c(3.0,1.8,1.8,1.0),oma=c(2,1.5,1.5,1.5))
  res.denit <- residuals(denit)
  qqnorm(res.denit)
  qqline(res.denit)
  hist(res.denit,freq=FALSE)
  plot(fitted(denit),residuals(denit))
  shapiro.test(res.denit)                                         # test assumption of normality
  car::leveneTest(log(DenitRate) ~ Date, data = data)             # test assumption of homogeneity of variance
  
  # Log transformation improves the normality and homoscedasticity of the model residuals
  # Q-Q plot looks OK, and there's no obvious patterning in the residuals plot. 
    
  
## ============================================================================= ## 
##                             2. N2O PRODUCTION RATE                            ##
## ============================================================================= ## 
  
## Plot data
  
  # Plot N2O production rates, including possible main (Month, VegZone), random (SampleID), and interactive effects. 
  # Note y-axis is presented on a log-scale
  plot_grid(
    data %>% ggplot() + geom_boxplot(aes(x=Date,y=(N2ORate+1)))+scale_y_log10()+theme_bw(),
    data %>% ggplot() + geom_boxplot(aes(x=VegZone,y=(N2ORate+1)))+scale_y_log10()+theme_bw(),
    data %>% ggplot() + geom_boxplot(aes(x=Sample.ID,y=(N2ORate+1),color = VegZone))+scale_y_log10()+theme_bw()+scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3"))+theme(legend.position="none"),
    # interaction plot:
    ggplot(data, aes(x = Date, y = log10(N2ORate+1), color = VegZone)) + geom_boxplot(size=0.5) +
      stat_summary(fun = mean, geom = "line", aes(group = VegZone), size = 1.2) + theme_bw() +
      scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3")),
    ncol=2,nrow=2,align="h")
  
## Check - is the experimental design balanced?
  table(data[which(!is.na(data$N2ORate)),]$VegZone, data[which(!is.na(data$N2ORate)),]$Date)
  
## Analyze data across season and salt marsh zone:
  
  #Two-way ANOVA with SampleID (i.e. sub-plot) as a random effect, assuming that samples from the same sampleID are correlated over time. 
  #Note that N2O production data were natural log(x+1)-transformed prior to analysis.
  n2o = lme4::lmer(log(N2ORate+1) ~ Date*VegZone + (1|Sample.ID),data=data,REML=TRUE)
  anova(n2o)
  # run n2o model again using lmerTest Satterthwaite approximation of degrees of freedom in lmerTest:
  options(contrasts = c("contr.sum", "contr.poly"))
  n2o = lmerTest::lmer(log(N2ORate+1) ~ Date*VegZone + (1|Sample.ID),data=data,REML=TRUE)
  anova(n2o)
  
  # Pairwise comparisons by month*vegetation zone: 
  posthoc.n2o.emm <- emmeans(n2o,c("Date","VegZone"),adjust="tukey",conf=.95)
  pairs(posthoc.n2o.emm,type="response")
  
  # Compare effect size:
  summary(posthoc.n2o.emm,type="response")

  # Calculate marginal and conditional r2:
  piecewiseSEM::rsquared(n2o)

## Check model assumptions:
  
  # Create vector of unique interaction combinations of Month x VegZone:
  data$unique.trt <- interaction(data$Date,data$VegZone)
  
  # Check model assumptions:
  par(mfrow=c(2,2),mar=c(3.0,1.8,1.8,1.0),oma=c(2,1.5,1.5,1.5))
  res.n2o <- residuals(n2o)
  qqnorm(res.n2o)
  qqline(res.n2o)
  hist(res.n2o,freq=FALSE)
  plot(fitted(n2o),residuals(n2o))
  shapiro.test(res.n2o)                                             # test assumption of normality
  car::leveneTest(log(N2ORate+1) ~ unique.trt, data = data)         # test assumption of homogeneity of variance
  
  # Log transformation improves the normality and homoscedasticity of the model residuals
  # Q-Q plot looks OK, and there's no obvious patterning in the residuals plot (although points with N2ORate=0 have leverage). 
  

## ============================================================================= ## 
##                                3. N2O YIELD                                   ##
## ============================================================================= ## 
  
## Plot data
  
  # Plot N2O yield, including possible main (Month, VegZone), random (SampleID), and interactive effects. 
  # Note y-axis is presented on a log-scale
  plot_grid(
    data %>% ggplot() + geom_boxplot(aes(x=Date,y=(N2O.yield+0.1)))+scale_y_log10()+theme_bw(),
    data %>% ggplot() + geom_boxplot(aes(x=VegZone,y=(N2O.yield+0.1)))+scale_y_log10()+theme_bw(),
    data %>% ggplot() + geom_boxplot(aes(x=Sample.ID,y=(N2O.yield+0.1),color=VegZone))+scale_y_log10()+theme_bw()+scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3"))+theme(legend.position="none"),
    # interaction plot:
    ggplot(data, aes(x = Date, y = log(N2O.yield+0.1), color = VegZone)) + geom_boxplot(size=0.5) +
      stat_summary(fun = mean, geom = "line", aes(group = VegZone), size = 1.2) + theme_bw() +
      scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3")),
    ncol=2,nrow=2,align="h")
  
## Check - is the experimental design balanced?
  table(data[which(!is.na(data$N2O.yield)),]$VegZone, data[which(!is.na(data$N2O.yield)),]$Date)
  
## Analyze data across season and salt marsh zone:
  
  #Two-way ANOVA with SampleID (i.e. sub-plot) as a random effect, assuming that samples from the same sampleID are correlated over time. 
  #Note that N2O yield data were natural log(x+0.1)-transformed prior to analysis.
  n2oyield = lme4::lmer(log(N2O.yield+0.1) ~ Date*VegZone + (1|Sample.ID),data=data,REML=TRUE)
  anova(n2oyield)
  # run n2oyield model again using lmerTest Satterthwaite approximation of degrees of freedom in lmerTest:
  options(contrasts = c("contr.sum", "contr.poly"))
  n2oyield = lmerTest::lmer(log(N2O.yield+0.1) ~ Date*VegZone + (1|Sample.ID),data=data,REML=TRUE)
  anova(n2oyield)
  
  # Pairwise comparisons by month*vegetation zone: 
  posthoc.n2oyield.emm <- emmeans(n2oyield,c("Date","VegZone"),adjust="tukey",conf=.95)
  pairs(posthoc.n2oyield.emm,type="response")
  
  # Compare effect size:
  summary(posthoc.n2oyield.emm,type="response")

  # Calculate marginal and conditional r2:
  piecewiseSEM::rsquared(n2oyield)
  
## Check model assumptions:
  par(mfrow=c(2,2),mar=c(3.0,1.8,1.8,1.0),oma=c(2,1.5,1.5,1.5))
  res.n2oyield <- residuals(n2oyield)
  qqnorm(res.n2oyield)
  qqline(res.n2oyield)
  hist(res.n2oyield,freq=FALSE)
  plot(fitted(n2oyield),residuals(n2oyield))
  shapiro.test(res.n2oyield)                                       # test assumption of normality
  car::leveneTest(log(N2O.yield+0.1) ~ unique.trt, data = data)    # test assumption of homogeneity of variance
  
  # Log transformation improves the normality and homoscedasticity of the model residuals
  
## Look at N2O yields as scatterplots between potential denitrification rate and N2O production rate:
  
  # Plot the overall relationship (across month x vegetation zone combinations):
  data[which(data$N2ORate>0),] %>% ggplot() + geom_point(aes(x=DenitRate,y=N2ORate,color=VegZone),size=1.3)  +
    geom_smooth(aes(x=DenitRate,y=N2ORate,color=VegZone),method="lm",se=FALSE)+
    theme_cowplot()+
    scale_x_log10()+
    scale_y_log10()+
    NULL
  
  # Overall relationship (across month x vegetation zone combinations):
  Nreg.all <- lm(log(N2ORate)~log(DenitRate),data[which(data$N2ORate>0),])
  par(mfrow = c(2, 2))
  plot(Nreg.all)
  summary(Nreg.all)
  
  par(mfrow=c(2,2),mar=c(3.0,1.8,1.8,1.0),oma=c(2,1.5,1.5,1.5))
  res.Nreg <- residuals(Nreg.all)
  qqnorm(res.Nreg)
  qqline(res.Nreg)
  hist(res.Nreg,freq=FALSE)
  plot(fitted(Nreg.all),residuals(Nreg.all))
  shapiro.test(res.Nreg)                        # test assumption of normality

  # Regressions between potential denitrification and N2O production rates by vegetation zone:
  Nreg.byveg = data[which(data$N2ORate>0),] %>% 
               group_by(VegZone) %>%
               do(fit = lm(log(N2ORate)~log(DenitRate),data= .))
  #Nreg.byveg.coef = Nreg.byveg %>% tidy(fit)
  Nreg.byveg.coef = Nreg.byveg %>% glance(fit)  
  # get r-squared values:
  print(Nreg.byveg.coef)
  # get slopes:
  print(Nreg.byveg$fit)
  
  # ANCOVA with potential denitrification as covariate and salt marsh zone as categorical independent variable:
  options(contrasts = c("contr.treatment", "contr.poly"))
  Nreg.ancova=aov(log(N2ORate)~log(DenitRate)*VegZone, data=data[which(data$N2ORate>0),])
  anova(Nreg.ancova) # interaction not significant, indicating no heterogeneity in slopes
  
  par(mfrow=c(2,2),mar=c(3.0,1.8,1.8,1.0),oma=c(2,1.5,1.5,1.5))
  res.Nrega <- residuals(Nreg.ancova)
  qqnorm(res.Nrega)
  qqline(res.Nrega)
  hist(res.Nrega,freq=FALSE)
  plot(fitted(Nreg.ancova),residuals(Nreg.ancova))
  shapiro.test(res.Nrega)                          # test assumption of normality
  
  
## ============================================================================= ## 
##                  CLIMATE AND ENVIRONMENTAL COVARIATE DATA                     ##
## ============================================================================= ## 
  
## Read in data (generated from code chunk above)
  climate.dat <- readRDS("./output/data_processed/Table1_EnvironmentalConditions.rds")
  climate.dat.long <- climate.dat %>% tidyr::gather(Antcdt_rainfall,value_mm, roll_7day:roll_30day)
  
## Plot environmental covariate data
  
  # Boxplots of environmental variables over seasons:
  plot_grid(climate.dat.long %>% ggplot() + 
      geom_boxplot(aes(x=Month.N,y=value_mm))+
      geom_point(aes(x=Month.N,y=value_mm,color=Antcdt_rainfall),size=2)+
      theme_bw() + 
      labs(x="Month",y=expression(Antecedent~precip~(mm))),
    data %>% ggplot() + geom_boxplot(aes(x=Date,y=LOI))+geom_point(aes(x=Date,y=LOI,color=VegZone))+theme_bw(),
    data %>% ggplot() + geom_boxplot(aes(x=Date,y=SO4_wet_normalized))+geom_point(aes(x=Date,y=SO4_wet_normalized,color=VegZone))+theme_bw()+theme(legend.position="none"),
    data %>% ggplot() + geom_boxplot(aes(x=Date,y=NH4.mg.N.per.g.dry.soil))+geom_point(aes(x=Date,y=NH4.mg.N.per.g.dry.soil,color=VegZone))+theme_bw()+theme(legend.position="none"),
    data %>% ggplot() + geom_boxplot(aes(x=Date,y=Soil.Moisture.Fraction))+geom_point(aes(x=Date,y=Soil.Moisture.Fraction,color=VegZone))+theme_bw()+theme(legend.position="none"),
    data %>% ggplot() + geom_boxplot(aes(x=Date,y=Belowground.Biomass.Weight..g.))+geom_point(aes(x=Date,y=Belowground.Biomass.Weight..g.,color=VegZone))+labs(y=expression(Belowground~biomass~(g)))+theme_bw()+theme(legend.position="none"),
    ncol=2
  )
  
  # Scatterplots between environmental variables and potential denitrification rates:
  plot_grid(
    data %>% ggplot() + geom_point(aes(x=SO4_wet_normalized,y=log(DenitRate),color=VegZone))+theme_bw()+theme(legend.position="none"),
    data %>% ggplot() + geom_point(aes(x=NH4_wet_normalized,y=log(DenitRate),color=VegZone))+theme_bw()+theme(legend.position="none"),
    data %>% ggplot() + geom_point(aes(x=LOI,y=log(DenitRate),color=VegZone))+theme_bw()+theme(legend.position="none"),
    data %>% ggplot() + geom_point(aes(x=Soil.Moisture.Fraction,y=log(DenitRate),color=VegZone))+theme_bw()+theme(legend.position="none"),
    data %>% ggplot() + geom_point(aes(x=Belowground.Biomass.Weight..g.,y=log(DenitRate),color=VegZone))+theme_bw()+labs(x=expression(Belowground~biomass~(g))),
    ncol=3)
  
  # Scatterplots between environmental variables and N2O production rates:
  plot_grid(
    data %>% ggplot() + geom_point(aes(x=SO4_wet_normalized,y=log(N2ORate+1), color=VegZone))+theme_bw()+theme(legend.position="none"),
    data %>% ggplot() + geom_point(aes(x=NH4_wet_normalized,y=log(N2ORate+1),color=VegZone))+theme_bw()+theme(legend.position="none"),
    data %>% ggplot() + geom_point(aes(x=LOI,y=log(N2ORate+1),color=VegZone))+theme_bw()+theme(legend.position="none"),
    data %>% ggplot() + geom_point(aes(x=Soil.Moisture.Fraction,y=log(N2ORate+1),color=VegZone))+theme_bw()+theme(legend.position="none"),
    data %>% ggplot() + geom_point(aes(x=Belowground.Biomass.Weight..g.,y=log(N2ORate+1),color=VegZone))+theme_bw()+labs(x=expression(Belowground~biomass~(g))),
    ncol=3)
  
  # Scatterplots between environmental variables and N2O yield:
  plot_grid(
    data %>% ggplot() + geom_point(aes(x=SO4_wet_normalized,y=log(N2O.yield+0.1),color=VegZone))+theme_bw()+theme(legend.position="none"),
    data %>% ggplot() + geom_point(aes(x=NH4_wet_normalized,y=log(N2O.yield+0.1),color=VegZone))+theme_bw()+theme(legend.position="none"),
    data %>% ggplot() + geom_point(aes(x=LOI,y=log(N2O.yield+0.1),color=VegZone))+theme_bw()+theme(legend.position="none"),
    data %>% ggplot() + geom_point(aes(x=Soil.Moisture.Fraction,y=log(N2O.yield+0.1),color=VegZone))+theme_bw()+theme(legend.position="none"),
    data %>% ggplot() + geom_point(aes(x=Belowground.Biomass.Weight..g.,y=log(N2O.yield+0.1),color=VegZone))+theme_bw()+labs(x=expression(Belowground~biomass~(g))),
    ncol=3)
  
  
## ============================================================================= ## 
##                              1. SOIL SULFATE                                  ##
## ============================================================================= ## 
  
## Plot data
  
  # Plot soil SO4 (mg SO4 per g wet soil), including possible main (Month, VegZone), random (SampleID), and interactive effects. 
  plot_grid(
    data %>% ggplot() + geom_boxplot(aes(x=Date,y=SO4_wet_normalized))+theme_bw(),
    data %>% ggplot() + geom_boxplot(aes(x=VegZone,y=SO4_wet_normalized))+theme_bw(),
    data %>% ggplot() + geom_boxplot(aes(x=Sample.ID,y=SO4_wet_normalized,color=VegZone))+theme_bw()+scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3"))+theme(legend.position="none"),
    # interaction plot:
    ggplot(data, aes(x = Date, y = SO4_wet_normalized, color = VegZone)) + geom_boxplot(size=0.5) +
      stat_summary(fun = mean, geom = "line", aes(group = VegZone), size = 1.2) + theme_bw() +
      scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3")),
    ncol=2,nrow=2,align="h")
  
## Check - is the experimental design balanced?
  table(data[which(!is.na(data$SO4_wet_normalized)),]$VegZone, data[which(!is.na(data$SO4_wet_normalized)),]$Date)
  
## Analyze data across season and salt marsh zone:
  
  #Two-way ANOVA with SampleID (i.e. sub-plot) as a random effect, assuming that samples from the same sampleID are correlated over time. 
  so4 = lme4::lmer(SO4_wet_normalized ~ Date*VegZone + (1|Sample.ID),data=data,REML=TRUE,na.action = na.exclude)
  anova(so4)
  # run so4 model again using lmerTest Satterthwaite approximation of degrees of freedom in lmerTest:
  options(contrasts = c("contr.sum", "contr.poly"))
  so4 <- lmerTest::lmer(SO4_wet_normalized ~ Date*VegZone + (1|Sample.ID), data = data,REML = TRUE,na.action=na.exclude)
  anova(so4)
  
  # Pairwise comparisons by month*vegetation zone: 
  posthoc.so4.emm <- emmeans(so4,c("Date","VegZone"),adjust="tukey",conf=.95)
  pairs(posthoc.so4.emm,type="response")
  
  # Compare effect sizes: 
  summary(posthoc.so4.emm,type="response")
  
## Check model assumptions:
  par(mfrow=c(2,2),mar=c(3.0,1.8,1.8,1.0),oma=c(2,1.5,1.5,1.5))
  res.so4 <- residuals(so4)
  qqnorm(res.so4)
  qqline(res.so4)
  hist(res.so4,freq=FALSE)
  plot(fitted(so4),residuals(so4))
  shapiro.test(res.so4)                                            # test assumption of normality
  car::leveneTest(SO4_wet_normalized ~ unique.trt, data = data)    # test assumption of homogeneity of variance
  
  # Q-Q plot looks OK, and there's no obvious patterning in the residuals plot.
  
  
## ============================================================================= ## 
##                             2. SOIL AMMONIUM                                  ##
## ============================================================================= ## 
  
## Plot data
  
  # Plot soil NH4 (mg N per g wet soil), including possible main (Month, VegZone), random (SampleID), and interactive effects. 
  # Note y-axis is presented on a log-scale
  plot_grid(
    data %>% ggplot() + geom_boxplot(aes(x=Date,y=NH4_wet_normalized))+scale_y_log10()+theme_bw(),
    data %>% ggplot() + geom_boxplot(aes(x=VegZone,y=NH4_wet_normalized))+scale_y_log10()+theme_bw(),
    data %>% ggplot() + geom_boxplot(aes(x=Sample.ID,y=NH4_wet_normalized,color=VegZone))+scale_y_log10()+theme_bw()+scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3"))+theme(legend.position="none"),
    # interaction plot:
    ggplot(data, aes(x = Date, y = log(NH4_wet_normalized), color = VegZone)) + geom_boxplot(size=0.5) +
      stat_summary(fun = mean, geom = "line", aes(group = VegZone), size = 1.2) + theme_bw() +
      scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3")),
    ncol=2,nrow=2,align="h")
  
## Check - is the experimental design balanced?
  table(data[which(!is.na(data$NH4_wet_normalized)),]$VegZone, data[which(!is.na(data$NH4_wet_normalized)),]$Date)
  
## Analyze data across season and salt marsh zone:
  
  #Two-way ANOVA with SampleID (i.e. sub-plot) as a random effect, assuming that samples from the same sampleID are correlated over time. 
  #Note that soil NH4 data were natural log-transformed prior to analysis.  
  nh4 = lme4::lmer(log(NH4_wet_normalized) ~ Date * VegZone + (1|Sample.ID),data=data,REML=TRUE)
  anova(nh4)
  # run nh4 model again using lmerTest Satterthwaite approximation of degrees of freedom in lmerTest:
  options(contrasts = c("contr.sum", "contr.poly"))
  nh4 <- lmerTest::lmer(log(NH4_wet_normalized) ~ Date*VegZone + (1|Sample.ID), data = data,REML = TRUE)
  anova(nh4)
  
  # Pairwise comparisons by month*vegetation zone: 
  posthoc.nh4.emm <- emmeans(nh4,c("Date"),adjust="tukey",conf=.95)
  pairs(posthoc.nh4.emm,type="response")
  
  # Compare effect sizes: 
  summary(posthoc.nh4.emm,type="response")
  
## Check model assumptions:
  par(mfrow=c(2,2),mar=c(3.0,1.8,1.8,1.0),oma=c(2,1.5,1.5,1.5))
  res.nh4 <- residuals(nh4)
  qqnorm(res.nh4)
  qqline(res.nh4)
  hist(res.nh4,freq=FALSE)
  plot(fitted(nh4),residuals(nh4))
  shapiro.test(res.nh4)                                                 # test assumption of normality
  car::leveneTest(log(NH4_wet_normalized) ~ unique.trt, data = data)    # test assumption of homogeneity of variance
  
  # Log transformation improves the normality and homoscedasticity of the model residuals
  # Q-Q plot looks OK, and there's no obvious patterning in the residuals plot
  
  
## ============================================================================= ## 
##                             3. SOIL MOISTURE                                  ##
## ============================================================================= ## 
  
## Plot data
  
  # Plot soil moisture fraction, including possible main (Month, VegZone), random (SampleID), and interactive effects. 
  plot_grid(
    data %>% ggplot() + geom_boxplot(aes(x=Date,y=Soil.Moisture.Fraction))+theme_bw(),
    data %>% ggplot() + geom_boxplot(aes(x=VegZone,y=Soil.Moisture.Fraction))+theme_bw(),
    data %>% ggplot() + geom_boxplot(aes(x=Sample.ID,y=Soil.Moisture.Fraction,color=VegZone))+theme_bw()+scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3"))+theme(legend.position="none"),
    # interaction plot:
    ggplot(data, aes(x = Date, y = (Soil.Moisture.Fraction), color = VegZone)) + geom_boxplot(size=0.5) +
      stat_summary(fun = mean, geom = "line", aes(group = VegZone), size = 1.2) + theme_bw() +
      scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3")),
    ncol=2,nrow=2,align="h")
  
## Check - is the experimental design balanced?
  table(data[which(!is.na(data$Soil.Moisture.Fraction)),]$VegZone, data[which(!is.na(data$Soil.Moisture.Fraction)),]$Date)
  
## Analyze data across season and salt marsh zone:
  
  #Two-way ANOVA with SampleID (i.e. sub-plot) as a random effect, assuming that samples from the same sampleID are correlated over time. 
  soilmoisture = lme4::lmer(Soil.Moisture.Fraction ~ Date*VegZone + (1|Sample.ID),data=data,REML=TRUE)
  anova(soilmoisture)
  # run soilmoisture model again using lmerTest Satterthwaite approximation of degrees of freedom in lmerTest:
  options(contrasts = c("contr.sum", "contr.poly"))
  soilmoisture <- lmerTest::lmer(Soil.Moisture.Fraction ~ Date*VegZone + (1|Sample.ID), data = data,REML = TRUE)
  anova(soilmoisture)
  
## Check model assumptions:
  par(mfrow=c(2,2),mar=c(3.0,1.8,1.8,1.0),oma=c(2,1.5,1.5,1.5))
  res.SM <- residuals(soilmoisture)
  qqnorm(res.SM)
  qqline(res.SM)
  hist(res.SM,freq=FALSE)
  plot(fitted(soilmoisture),residuals(soilmoisture))
  shapiro.test(res.SM)                                                 # test assumption of normality
  car::leveneTest(Soil.Moisture.Fraction ~ unique.trt, data = data)    # test assumption of homogeneity of variance
  
  # Q-Q plot looks OK, and there's no obvious patterning in the residuals plot
  
## Estimate variation in soil moisture by vegetation zone:
  data %>% group_by(VegZone) %>%
           summarize(mean.SM = mean(Soil.Moisture.Fraction,na.rm=T),
                     sd.SM = sd(Soil.Moisture.Fraction,na.rm=T)) %>%
           mutate(cv.SM = (sd.SM/mean.SM)*100)
  
  
## ============================================================================= ## 
##                          4. BELOWGROUND BIOMASS                               ##
## ============================================================================= ## 
  
## Plot data
  
  # Plot belowground biomass, including possible main (Month, VegZone), random (SampleID), and interactive effects. 
  plot_grid(
    data %>% ggplot() + geom_boxplot(aes(x=Date,y=Belowground.Biomass..g...m2.))+theme_bw(),
    data %>% ggplot() + geom_boxplot(aes(x=VegZone,y=Belowground.Biomass..g...m2.))+theme_bw(),
    data %>% ggplot() + geom_boxplot(aes(x=Sample.ID,y=Belowground.Biomass..g...m2.,color=VegZone))+theme_bw()+scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3"))+theme(legend.position="none"),
    # interaction plot:
    ggplot(data, aes(x = Date, y = (Belowground.Biomass..g...m2.), color = VegZone)) + geom_boxplot(size=0.5) +
      stat_summary(fun = mean, geom = "line", aes(group = VegZone), size = 1.2) + theme_bw() +
      scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3")),
    ncol=2,nrow=2,align="h")
  
## Check - is the experimental design balanced?
  table(data[which(!is.na(data$Belowground.Biomass..g...m2.)),]$VegZone, data[which(!is.na(data$Belowground.Biomass..g...m2.)),]$Date)

## Analyze data across season and salt marsh zone:
  
  #Two-way ANOVA with SampleID (i.e. sub-plot) as a random effect, assuming that samples from the same sampleID are correlated over time. 
  biomass = lme4::lmer(Belowground.Biomass..g...m2. ~ Date*VegZone + (1|Sample.ID),data=data,REML=TRUE)
  anova(biomass)
  # run biomass model again using lmerTest Satterthwaite approximation of degrees of freedom in lmerTest:
  options(contrasts = c("contr.sum", "contr.poly"))
  biomass <- lmerTest::lmer(Belowground.Biomass..g...m2. ~ Date*VegZone + (1|Sample.ID), data = data,REML = TRUE)
  anova(biomass)
  
  # Pairwise comparisons by month*vegetation zone: 
  posthoc.biomass.emm <- emmeans(biomass,c("VegZone"),adjust="tukey",conf=.95)
  pairs(posthoc.biomass.emm,type="response")
  
  # Compare effect sizes:
  summary(posthoc.biomass.emm,type="response")
  
## Check model assumptions:
  par(mfrow=c(2,2),mar=c(3.0,1.8,1.8,1.0),oma=c(2,1.5,1.5,1.5))
  res.biomass <- residuals(biomass)
  qqnorm(res.biomass)
  qqline(res.biomass)
  hist(res.biomass,freq=FALSE)
  plot(fitted(biomass),residuals(biomass))
  shapiro.test(res.biomass)                                                  # test assumption of normality
  car::leveneTest(Belowground.Biomass..g...m2. ~ unique.trt, data = data)    # test assumption of homogeneity of variance
  
  # Q-Q plot looks OK, and there's no obvious patterning in the residuals plot
  
  
## ============================================================================= ## 
##                          FIGURE 3: N CYCLING RATES                            ##
## ============================================================================= ## 
  
  # Rename variables for figure labels:
  data$Month2 <- factor(data$Date,levels = c("May","June","July","August","October"),
                        labels = c("May","June","July","Aug","Oct"))
  data$VegZone2 <- factor(data$VegZone,levels=c("LM","HM","PH"),
                          labels=c("Low Marsh","High Marsh", "Phragmites"))
  
  # Create denitrification panels:
  fig3.potdenit <- ggplot(data) + 
    geom_boxplot(aes(x=Month2,y=DenitRate),alpha=.7)+
    facet_grid(. ~ VegZone2)+
    scale_y_log10()+
    labs(y=expression(atop(Potential~denitrification, (ng~N~"/"~hr~"/"~g~dry~soil))),
         x="Month")+
    theme_bw() +
    theme(legend.position="right",
          axis.title = element_text(size=11),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size=10,colour = "black"),
          #axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 6, b = 0, l = 0)),
          strip.text.x = element_text(size=12,margin = margin(.1, 0, .1, 0, "cm")))
  
  # Create N2O production panels:
  fig3.n2oprod <- ggplot(data) + 
    geom_boxplot(aes(x=Month2,y=(N2ORate+1)),alpha=.7)+
    facet_grid(. ~ VegZone2)+
    scale_y_log10()+
    labs(y=expression(atop(paste(N[2],"O")~production, (ng~N~"/"~hr~"/"~g~dry~soil))),
         x="Month")+
    theme_bw() +
    theme(legend.position="right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=11),
          axis.text = element_text(size=10,colour = "black"),
          #axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          #axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r = 6, b = 0, l = 0)),
          #strip.text.x = element_text(size=14,margin = margin(.1, 0, .1, 0, "cm")))
          strip.text.x = element_blank())
  
  # Create N2O yield panels:
  fig3.n2oyield <- ggplot(data) + 
    geom_boxplot(aes(x=Month2,y=(N2O.yield+0.01)),alpha=.7)+
    facet_grid(. ~ VegZone2)+
    scale_y_log10()+
    labs(y=expression(paste(N[2],"O")~yield),
         x="Month")+
    theme_bw() +
    theme(legend.position="right",
          axis.title = element_text(size=11),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size=10,colour = "black"),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
          #axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 6, b = 0, l = 0)),
          #strip.text.x = element_text(size=12,margin = margin(.1, 0, .1, 0, "cm")),
          strip.text.x = element_blank())
  
  fig3 <- fig3.potdenit / fig3.n2oprod / fig3.n2oyield
  
  jpeg("./output/figures/Figure3.jpg",width = 6.75,height=6.75,units = "in",res=300)
  print(fig3)
  dev.off()
  
  

## ============================================================================= ## 
##                          FIGURE 4: SOIL CONDITIONS                            ##
## ============================================================================= ##    
  
  fig4.moisture <- ggplot(data) + 
    geom_boxplot(aes(x=Month2,y=Soil.Moisture.Fraction),alpha=.7)+
    facet_grid(. ~ VegZone2)+
    coord_cartesian(ylim=c(0.5,1))+
    labs(y=expression(Soil~moisture),
         x="Month")+
    theme_bw() +
    theme(legend.position="right",
          axis.title = element_text(size=11),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size=10,colour = "black"),
          #axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
          strip.text.x = element_text(size=12,margin = margin(.1, 0, .1, 0, "cm")))
  
  fig4.so4 <- ggplot(data) + 
    geom_boxplot(aes(x=Month2,y=SO4_wet_normalized),alpha=.7)+
    facet_grid(. ~ VegZone2)+
    labs(x="Month",y=expression(atop(Sulfate~concentration, (mg~SO[4]^2^"-"~"/"~g~wet~soil))))+
    theme_bw() +
    theme(legend.position="right",
          axis.title = element_text(size=11),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size=10,colour = "black"),
          #axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
          strip.text.x = element_blank())
  
  fig4.nh4 <- ggplot(data) + 
    geom_boxplot(aes(x=Month2,y=NH4_wet_normalized),alpha=.7)+
    facet_grid(. ~ VegZone2)+
    scale_y_log10()+
    labs(x="Month",y=expression(atop(Ammonium~concentration, (mg~NH[4]^"+"~"/"~g~wet~soil))))+
    theme_bw() +
    theme(legend.position="right",
          axis.title = element_text(size=11),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size=10,colour = "black"),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
          #axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
          strip.text.x = element_blank())
  
  fig4 <- fig4.moisture / fig4.so4 / fig4.nh4
  
  jpeg("./output/figures/Figure4.jpg",width = 6.75,height=6.75,units = "in",res=300)
  print(fig4)
  dev.off()

