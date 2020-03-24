#### import libraries ####
# data wrangling

library(tidyverse)
library(reshape)
library(data.table)
library(xlsx)
library(gridExtra)
library(grid)
library(chron)
library(devtools)
library(SOfun)
library(usethis)
library(reprex)

# data visualization 

library(GGally)
library(RColorBrewer)
library(proj4)
library(leaflet)
library(leaflet.minicharts)
library(RColorBrewer)
library(mapview)
library(htmlwidgets)
library(corrplot)
library(mice)
library(VIM)
library(ggmosaic)
library(esquisse)

# data analysis - machine learning - statistics

library(Rcpp)
library(vegan)
library(cluster)
library(MuMIn) # R2c and R2m
library(nlme)
library(gamm4)
library(plotrix)
library(lme4)
library(dunn.test)  # Kruskal Wallis test
library(car)
library(psych)
library(psycho)
library(remotes)
library(mlr)
library(randomForest)
library(Metrics)
library(Hmisc)
library(xgboost)
library(checkmate)
library(ranger)



#### import data #### 

aday <- read.csv("20190516_24h.csv", header = TRUE)
aday_mis <- aday[0,]

for (i in 1:ncol(aday_mis)){
    if (is.factor(aday_mis[,i])){
        aday_mis[,i] <- as.numeric(as.character(aday_mis[,i]))
    }
    aday_mis[1,i] <- sum(is.na(aday[,i]))
}
#### processing data 4C's ####

aday$Date <- as.Date(aday$Date, format = "%y-%m-%d")

aday$Time <- ordered(aday$Time, levels = moveMe(levels(aday$Time), "0:22 after 23:12; 1:20 after 0:22; 2:35 after 1:20"))

#### Correcting dissolved gas concentrations ####

V.headspace = 6 # mL 
V.aq        = 6 # mL

R        = 1.98719   # gas constant in cal / K mol
T        = aday$T_w        # assumed river temperature 
pressure = 2 # atm 
pressure2 = 1

# Henry coefficients for N2O
A.n2o    = -180.95   
B.n2o    =  13205.8     
C.n2o    =  20.0399   
D.n2o    =  0.0238544


# HENRY's LAW
# H = 1/exp[{-a+b/(T+273)+C*ln(T+273)+d*(T+273)}]/R	
H.n2o    = 1/exp((A.n2o + B.n2o / (T + 273.15) + C.n2o * log(T + 273.15) + D.n2o * (T + 273.15)) / R)	

Cg.n2o   = aday$Dis_N2O_1 * 10^-6 * pressure  # vol atm N2O / vol sample; Cg


Ca.n2o   = 55.5 * (Cg.n2o/H.n2o) * 44 * 10^3  # mg N2O/L H2O
Cah.n2o  = (V.headspace/V.aq * Cg.n2o * (44/22.4) * (273.15/(T + 273.15))*10^3) # mg N2O/L H2O

n2o.aq   = (Ca.n2o + Cah.n2o) * 10^3 # ug N2O/L H2O
aday$Dis_N2O_cor <- n2o.aq


# Henry coefficients for ch4

A.ch4    = -365.183    
B.ch4    =  18106.7      
C.ch4    =  49.7554    
D.ch4    = -0.00028503 

# HENRY's LAW
# H = 1/exp[{-a+b/(T+273)+C*ln(T+273)+d*(T+273)}]/R	
H.ch4    = 1/exp((A.ch4 + B.ch4 / (T + 273.15) + C.ch4 * log(T + 273.15) + D.ch4 * (T + 273.15)) / R)	
Cg.ch4   = aday$Dis_CH4_1 * 10^-6  * pressure # vol atm ch4 / vol sample; Cg



Ca.ch4   = 55.5 * (Cg.ch4/H.ch4) * 16 * 10^3  # mg ch4/L H2O
Cah.ch4  = (V.headspace/V.aq * Cg.ch4 * (16/22.4) * (273.15/(T + 273.15))*10^3) # mg ch4/L H2O

ch4.aq   = (Ca.ch4 + Cah.ch4)  * 10^3 # ug ch4/L H2O
ch4.aq

aday$Dis_CH4_cor <- ch4.aq

# Henry coefficients for co2

A.co2    = -317.658    
B.co2    =  1737.12     
C.co2    =  43.0607    
D.co2    = -0.00219107 


# HENRY's LAW
# H = 1/exp[{-a+b/(T+273)+C*ln(T+273)+d*(T+273)}]/R 
H.co2    = 1/exp((A.co2 + B.co2 / (T + 273.15) + C.co2 * log(T + 273.15) + D.co2 * (T + 273.15)) / R)   
Cg.co2   = aday$Dis_CO2_1 * 10^-6 * pressure   # vol atm co2 / vol sample; Cg


Ca.co2   = 55.5 * (Cg.co2/H.co2) * 44 * 10^3  # mg co2/L H2O
Cah.co2  = (V.headspace/V.aq * Cg.co2 * (44/22.4) * (273.15/(T + 273.15))*10^3) # mg co2/L H2O

co2.aq   = (Ca.co2 + Cah.co2)  # mg co2/L H2O
co2.aq

aday$Dis_CO2_cor <- co2.aq


#### Boxplot of all countinuous variables regarding ponds #### 

# This  is for make a list of graphs

boxplot_aday <- function(data, column, column2){
    ggplot(data) +
        geom_boxplot(aes_string(x = column2 , y = column)) +
        xlab(column2) +
        ylab(column) 
}

boxplot_aday <- lapply(colnames(aday)[6:26], boxplot_aday, data = aday, column2 = "Pond")

# print the graphs

lapply(boxplot_aday, print)

# save the graphs into folder

for (i in 1:length(boxplot_aday)){
    tiff(filename = paste("Boxplot_24_",colnames(aday)[6:26][i],".tiff", sep =""),units = 'px',height = 1800,width = 1800,res = 300,pointsize = 12)
    print(boxplot_aday[[i]])
    dev.off()
}

# Put graphs with similar variables together

ggsave("Boxplot_24_1.tiff",marrangeGrob(boxplot_aday[1:4],nrow = 2, ncol= 2, top = ''),
       units = 'cm', height = 20, width = 20, dpi = 300)

ggsave("Boxplot_24_2.tiff",marrangeGrob(boxplot_aday[5:8],nrow = 2, ncol= 2, top = ''),
       units = 'cm', height = 20, width = 20, dpi = 300)

ggsave("Boxplot_24_3.tiff",marrangeGrob(boxplot_aday[9:11],nrow = 1, ncol= 3, top = ''),
       units = 'cm', height = 15, width = 30, dpi = 300)

ggsave("Boxplot_24_4.tiff",marrangeGrob(boxplot_aday[12:17],nrow = 2, ncol= 3, top = ''),
       units = 'cm', height = 15, width = 30, dpi = 300)

ggsave("Boxplot_24_5.tiff",marrangeGrob(boxplot_aday[18:21],nrow = 2, ncol= 2, top = ''),
       units = 'cm', height = 20, width = 20, dpi = 300)


#### Line chart of all countinuous variables regarding times #### 

# This  is for make a list of graphs

line_aday <- function(data, column){
    ggplot(data, aes_string(x = "Time", y = column, group = "Pond")) +
        geom_line(aes(color = Pond)) + 
        geom_point(aes(color = Pond)) 
}

line_pond <- lapply(colnames(aday)[6:26], line_aday, data = aday)

# print the graphs

lapply(line_pond, print)

# Put graphs with similar variables together

ggsave("Line_chart_1.tiff",marrangeGrob(line_pond[1:4],nrow = 2, ncol= 2, top = '')
       , units = 'cm', height = 20, width = 60,dpi = 300)

ggsave("Line_chart_2.tiff",marrangeGrob(line_pond[5:8],nrow = 2, ncol= 2, top = '')
       , units = 'cm', height = 20, width = 60,dpi = 300)

ggsave("Line_chart_3.tiff",marrangeGrob(line_pond[9:11],nrow = 2, ncol= 2, top = '')
       , units = 'cm', height = 20, width = 60,dpi = 300)

ggsave("Line_chart_4.tiff",marrangeGrob(line_pond[12:17],nrow = 3, ncol= 2, top = '')
       , units = 'cm', height = 40, width = 60,dpi = 300)

ggsave("Line_chart_5.tiff",marrangeGrob(line_pond[18:21],nrow = 2, ncol= 2, top = '')
       , units = 'cm', height = 20, width = 60,dpi = 300)

#### Line chart of DG1 and DG2 through times ####

line_DG_pond <- lapply(colnames(aday)[36:41], line_aday, data = aday)

for (i in 1:length(line_DG_pond)){
    ggsave(filename = paste("Line_chart_DG", i, ".tiff", sep =""), line_DG_pond[[i]], units = 'cm', height = 20, width = 40,dpi = 300)
}

ggsave("Line_chart_DG_all.tiff", marrangeGrob(line_DG_pond[1:6],nrow = 3, ncol= 2, top = '')
       , units = 'cm', height = 40, width = 60,dpi = 300)

#### Line chart of DG1 and DG2 regarding times and ponds ####

aday2 <- cbind(aday[,3:5],aday[,36:41])
aday2 <- aday2 %>% gather(key = 'Dissolved_gases', value = 'Concentration', - Pond, - Time, - Sample)

ggsave("Linechart_24h_CO2.tiff", aday2 %>% filter(str_detect(Dissolved_gases, "Dis_CO2_1|Dis_CO2_2")) %>% 
           ggplot(aes(x = Time, y = Concentration, group = as.factor(Dissolved_gases))) +
           geom_line(aes(color = as.factor(Dissolved_gases))) + 
           geom_point(aes(color = as.factor(Dissolved_gases))) +
           xlab("Time") +
           ylab("CO2 (ppm)")+
           facet_wrap(.~ Pond)
       , units = 'cm', height = 20, width = 60,dpi = 300)

ggsave("Linechart_24h_CH4.tiff", aday2 %>% filter(str_detect(Dissolved_gases, "Dis_CH4_1|Dis_CH4_2")) %>% 
           ggplot(aes(x = Time, y = Concentration, group = as.factor(Dissolved_gases))) +
           geom_line(aes(color = as.factor(Dissolved_gases))) + 
           geom_point(aes(color = as.factor(Dissolved_gases))) +
           xlab("Time") +
           ylab("CH4 (ppm)")+
           facet_wrap(.~ Pond)
       , units = 'cm', height = 20, width = 60,dpi = 300)

ggsave("Linechart_24h_N2O.tiff", aday2 %>% filter(str_detect(Dissolved_gases, "Dis_N2O_1|Dis_N2O_2")) %>% 
           ggplot(aes(x = Time, y = Concentration, group = as.factor(Dissolved_gases))) +
           geom_line(aes(color = as.factor(Dissolved_gases))) + 
           geom_point(aes(color = as.factor(Dissolved_gases))) +
           xlab("Time") +
           ylab("N2O (ppm)")+
           facet_wrap(.~ Pond, scales = "free")
       , units = 'cm', height = 20, width = 60,dpi = 300)


#### Line chart of DG_cor through times ####

line_DG_cor_pond <- lapply(colnames(aday)[47:49], line_aday, data = aday)

for (i in 1:length(line_DG_cor_pond)){
    ggsave(filename = paste("Line_chart_DG_cor", i, ".tiff", sep =""), line_DG_cor_pond[[i]], units = 'cm', height = 20, width = 40,dpi = 300)
}

ggsave("Line_chart_DG_all_cor.tiff", marrangeGrob(line_DG_cor_pond[1:3],nrow = 3, ncol= 1, top = ''), units = 'cm', height = 20, width = 40,dpi = 300)

#### Line chart of DG_cor regarding times and ponds ####

aday2 <- cbind(aday[,3:5],aday[,47:49])
aday2 <- aday2 %>% gather(key = 'Dissolved_gases', value = 'Concentration', - Pond, - Time, - Sample)

ggsave("Linechart_24h_CO2_cor.tiff", aday2 %>% filter(str_detect(Dissolved_gases, "Dis_CO2_cor")) %>% 
           ggplot(aes(x = Time, y = Concentration, group = as.factor(Dissolved_gases))) +
           geom_line(aes(color = as.factor(Dissolved_gases))) + 
           geom_point(aes(color = as.factor(Dissolved_gases))) +
           xlab("Time") +
           ylab("CO2 (ug/L)")+
           facet_wrap(.~ Pond)
       , units = 'cm', height = 20, width = 60,dpi = 300)

ggsave("Linechart_24h_CH4_cor.tiff", aday2 %>% filter(str_detect(Dissolved_gases, "Dis_CH4_cor")) %>% 
           ggplot(aes(x = Time, y = Concentration, group = as.factor(Dissolved_gases))) +
           geom_line(aes(color = as.factor(Dissolved_gases))) + 
           geom_point(aes(color = as.factor(Dissolved_gases))) +
           xlab("Time") +
           ylab("CH4 (ug/L)")+
           facet_wrap(.~ Pond)
       , units = 'cm', height = 20, width = 60,dpi = 300)

ggsave("Linechart_24h_N2O_cor.tiff", aday2 %>% filter(str_detect(Dissolved_gases, "Dis_N2O_cor")) %>% 
           ggplot(aes(x = Time, y = Concentration, group = as.factor(Dissolved_gases))) +
           geom_line(aes(color = as.factor(Dissolved_gases))) + 
           geom_point(aes(color = as.factor(Dissolved_gases))) +
           xlab("Time") +
           ylab("N2O (ug/L)")+
           facet_wrap(.~ Pond, scales = "free")
       , units = 'cm', height = 20, width = 60,dpi = 300)


#### Boxplot of differences between DG1 and DG2 regarding ponds #### 

ggsave("Boxplot_24h_CO2.tiff", aday2 %>% filter(str_detect(Dissolved_gases, "Dis_CO2_1|Dis_CO2_2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("CO2"))

ggsave("Boxplot_24h_NO2.tiff", aday2 %>% filter(str_detect(Dissolved_gases, "Dis_N2O_1|Dis_N2O_2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("NO2"))

ggsave("Boxplot_24h_CH4.tiff", aday2 %>% filter(str_detect(Dissolved_gases, "Dis_CH4_1|Dis_CH4_2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("CH4"))

# Plot them for different ponds

ggsave("Boxplot_pond_24h_CO2.tiff", aday2 %>% filter(str_detect(Dissolved_gases, "Dis_CO2_1|Dis_CO2_2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("CO2") +
           facet_wrap(.~Pond))

ggsave("Boxplot_pond_24h_NO2.tiff", aday2 %>% filter(str_detect(Dissolved_gases, "Dis_N2O_1|Dis_N2O_2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("NO2") +
           facet_wrap(.~Pond))

ggsave("Boxplot_pond_24h_CH4.tiff", aday2 %>% filter(str_detect(Dissolved_gases, "Dis_CH4_1|Dis_CH4_2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("CH4") +
           facet_wrap(.~Pond))
#### Boxplot of DG_cor regarding ponds #### 

ggsave("Boxplot_24h_CO2_cor.tiff", aday2 %>% filter(str_detect(Dissolved_gases, "Dis_CO2_cor")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("CO2"))

ggsave("Boxplot_24h_NO2_cor.tiff", aday2 %>% filter(str_detect(Dissolved_gases, "Dis_N2O_cor")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("NO2"))

ggsave("Boxplot_24h_CH4_cor.tiff", aday2 %>% filter(str_detect(Dissolved_gases, "Dis_CH4_cor")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("CH4"))

# Plot them for different ponds

ggsave("Boxplot_pond_24h_CO2_cor.tiff", aday2 %>% filter(str_detect(Dissolved_gases, "Dis_CO2_cor")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("CO2") +
           facet_wrap(.~Pond))

ggsave("Boxplot_pond_24h_NO2_cor.tiff", aday2 %>% filter(str_detect(Dissolved_gases, "Dis_N2O_cor")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("NO2") +
           facet_wrap(.~Pond))

ggsave("Boxplot_pond_24h_CH4_cor.tiff", aday2 %>% filter(str_detect(Dissolved_gases, "Dis_CH4_cor")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("CH4") +
           facet_wrap(.~Pond))

#### Line chart of fluxes through times ####

line_DG_pond <- lapply(colnames(aday)[30:35], line_aday, data = aday)

for (i in 1:length(line_DG_pond)){
    ggsave(filename = paste("Line_chart_Fluxes", i, ".tiff", sep =""), line_DG_pond[[i]], units = 'cm', height = 20, width = 40,dpi = 300)
}

ggsave("Line_chart_Fluxes_all.tiff", marrangeGrob(line_DG_pond[1:6],nrow = 3, ncol= 2, top = '')
       , units = 'cm', height = 40, width = 60,dpi = 300)

#### Line chart of Fluxes (mmol.m-2.d-1) regarding times and ponds ####

aday2 <- cbind(aday[,3:5],aday[,30:32])
colnames(aday2)[4:6] <- c("CO2","N2O","CH4")
aday2 <- aday2 %>% gather(key = 'Fluxes', value = 'Concentration', - Pond, - Time, - Sample)

ggsave("Linechart_24h_CO2_flux.tiff", aday2 %>% filter(str_detect(Fluxes, "CO2")) %>% 
           ggplot(aes(x = Time, y = Concentration, group = as.factor(Fluxes))) +
           geom_line(aes(color = as.factor(Fluxes))) + 
           geom_point(aes(color = as.factor(Fluxes))) +
           xlab("Time") +
           ylab("CO2 (mmol.m-2.d-1)")+
           facet_wrap(.~ Pond)
       , units = 'cm', height = 20, width = 60,dpi = 300)

ggsave("Linechart_24h_CH4_flux.tiff", aday2 %>% filter(str_detect(Fluxes, "CH4")) %>% 
           ggplot(aes(x = Time, y = Concentration, group = as.factor(Fluxes))) +
           geom_line(aes(color = as.factor(Fluxes))) + 
           geom_point(aes(color = as.factor(Fluxes))) +
           xlab("Time") +
           ylab("CH4 (mmol.m-2.d-1)")+
           facet_wrap(.~ Pond)
       , units = 'cm', height = 20, width = 60,dpi = 300)

ggsave("Linechart_24h_N2O_flux.tiff", aday2 %>% filter(str_detect(Fluxes, "N2O")) %>% 
           ggplot(aes(x = Time, y = Concentration, group = as.factor(Fluxes))) +
           geom_line(aes(color = as.factor(Fluxes))) + 
           geom_point(aes(color = as.factor(Fluxes))) +
           xlab("Time") +
           ylab("N2O (mmol.m-2.d-1)")+
           facet_wrap(.~ Pond, scales = "free")
       , units = 'cm', height = 20, width = 60,dpi = 300)


#### Line chart of Fluxes (g.m-2.d-1) regarding times and ponds ####

aday3 <- cbind(aday[,3:5],aday[,33:35])
colnames(aday3)[4:6] <- c("CO2","N2O","CH4")
aday3 <- aday3 %>% gather(key = 'Fluxes', value = 'Concentration', - Pond, - Time, - Sample)

ggsave("Linechart_24h_CO2_flux_g.tiff", aday3 %>% filter(str_detect(Fluxes, "CO2")) %>% 
           ggplot(aes(x = Time, y = Concentration, group = as.factor(Fluxes))) +
           geom_line(aes(color = as.factor(Fluxes))) + 
           geom_point(aes(color = as.factor(Fluxes))) +
           xlab("Time") +
           ylab("CO2 (g.m-2.d-1)")+
           facet_wrap(.~ Pond)
       , units = 'cm', height = 20, width = 60,dpi = 300)

ggsave("Linechart_24h_CH4_flux_g.tiff", aday3 %>% filter(str_detect(Fluxes, "CH4")) %>% 
           ggplot(aes(x = Time, y = Concentration, group = as.factor(Fluxes))) +
           geom_line(aes(color = as.factor(Fluxes))) + 
           geom_point(aes(color = as.factor(Fluxes))) +
           xlab("Time") +
           ylab("CH4 (g.m-2.d-1)")+
           facet_wrap(.~ Pond)
       , units = 'cm', height = 20, width = 60,dpi = 300)

ggsave("Linechart_24h_N2O_flux_g.tiff", aday3 %>% filter(str_detect(Fluxes, "N2O")) %>% 
           ggplot(aes(x = Time, y = Concentration, group = as.factor(Fluxes))) +
           geom_line(aes(color = as.factor(Fluxes))) + 
           geom_point(aes(color = as.factor(Fluxes))) +
           xlab("Time") +
           ylab("N2O (g.m-2.d-1)")+
           facet_wrap(.~ Pond, scales = "free")
       , units = 'cm', height = 20, width = 60,dpi = 300)


#### Boxplot of Fluxes (mmol.m-2.d-1) regarding ponds #### 

ggsave("Boxplot_24h_CO2_Flux.tiff", aday2 %>% filter(str_detect(Fluxes, "CO2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (mmol.m-2.d-1)") +
           ylab("CO2"))

ggsave("Boxplot_24h_NO2_Flux.tiff", aday2 %>% filter(str_detect(Fluxes, "N2O")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (mmol.m-2.d-1)") +
           ylab("NO2"))

ggsave("Boxplot_24h_CH4_Flux.tiff", aday2 %>% filter(str_detect(Fluxes, "CH4")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (mmol.m-2.d-1)") +
           ylab("CH4"))

# Plot them for different ponds

ggsave("Boxplot_pond_24h_CO2_Flux.tiff", aday2 %>% filter(str_detect(Fluxes, "CO2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (mmol.m-2.d-1)") +
           ylab("CO2") +
           facet_wrap(.~Pond))

ggsave("Boxplot_pond_24h_NO2_Flux.tiff", aday2 %>% filter(str_detect(Fluxes, "N2O")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (mmol.m-2.d-1)") +
           ylab("NO2") +
           facet_wrap(.~Pond))

ggsave("Boxplot_pond_24h_CH4_Flux.tiff", aday2 %>% filter(str_detect(Fluxes, "CH4")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (mmol.m-2.d-1)") +
           ylab("CH4") +
           facet_wrap(.~Pond))

#### Boxplot of Fluxes (g.m-2.d-1) regarding ponds #### 

ggsave("Boxplot_24h_CO2_Flux_g.tiff", aday3 %>% filter(str_detect(Fluxes, "CO2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (g.m-2.d-1)") +
           ylab("CO2"))

ggsave("Boxplot_24h_NO2_Flux_g.tiff", aday3 %>% filter(str_detect(Fluxes, "N2O")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (g.m-2.d-1)") +
           ylab("NO2"))

ggsave("Boxplot_24h_CH4_Flux_g.tiff", aday3 %>% filter(str_detect(Fluxes, "CH4")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (g.m-2.d-1)") +
           ylab("CH4"))

# Plot them for different ponds

ggsave("Boxplot_pond_24h_CO2_Flux_g.tiff", aday3 %>% filter(str_detect(Fluxes, "CO2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (g.m-2.d-1)") +
           ylab("CO2") +
           facet_wrap(.~Pond))

ggsave("Boxplot_pond_24h_NO2_Flux_g.tiff", aday3 %>% filter(str_detect(Fluxes, "N2O")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (g.m-2.d-1)") +
           ylab("NO2") +
           facet_wrap(.~Pond))

ggsave("Boxplot_pond_24h_CH4_Flux_g.tiff", aday3 %>% filter(str_detect(Fluxes, "CH4")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (g.m-2.d-1)") +
           ylab("CH4") +
           facet_wrap(.~Pond))

#### Correlation coefficients #### 


ggpairs(aday[6:26],
        lower = list(continuous = wrap("smooth", color = "deepskyblue")),
        upper = list(continuous = wrap("cor", size = 3, color = "tomato"))
) + theme(panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank())


#### Comparisons test (t-test, ANOVA, etc) #### 
# lack of data --> non parameteric analysis

