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

WSP <- read.csv("20190516_WSP.csv", header = TRUE)
WSP_mis <- WSP[0,]

for (i in 1:ncol(WSP_mis)){
    if (is.factor(WSP_mis[,i])){
        WSP_mis[,i] <- as.numeric(as.character(WSP_mis[,i]))
    }
    WSP_mis[1,i] <- sum(is.na(WSP[,i]))
}

#### processing data 4C's ####

WSP$Line <- as.factor(WSP$Line)
WSP$Position <- as.factor(WSP$Position)

#### Correcting dissolved gas concentrations ####

V.headspace = 6 # mL 
V.aq        = 6 # mL

R        = 1.98719   # gas constant in cal / K mol
T        = WSP$T_w        # assumed river temperature 
pressure = 2 # atm 

# Henry coefficients for N2O
A.n2o    = -180.95   
B.n2o    =  13205.8     
C.n2o    =  20.0399   
D.n2o    =  0.0238544

# HENRY's LAW
# H = 1/exp[{-a+b/(T+273)+C*ln(T+273)+d*(T+273)}]/R	
H.n2o    = 1/exp((A.n2o + B.n2o / (T + 273.15) + C.n2o * log(T + 273.15) + D.n2o * (T + 273.15)) / R)	

Cg.n2o   = WSP$Dis_N2O_1 * 10^-6 * pressure  # vol atm N2O / vol sample; Cg


Ca.n2o   = 55.5 * (Cg.n2o/H.n2o) * 44 * 10^3  # mg N2O/L H2O
Cah.n2o  = (V.headspace/V.aq * Cg.n2o * (44/22.4) * (273.15/(T + 273.15))*10^3) # mg N2O/L H2O

n2o.aq   = (Ca.n2o + Cah.n2o) * 10^3 # ug N2O/L H2O
WSP$Dis_N2O_cor <- n2o.aq

# Henry coefficients for ch4

A.ch4    = -365.183    
B.ch4    =  18106.7      
C.ch4    =  49.7554    
D.ch4    = -0.00028503 

# HENRY's LAW
# H = 1/exp[{-a+b/(T+273)+C*ln(T+273)+d*(T+273)}]/R	
H.ch4    = 1/exp((A.ch4 + B.ch4 / (T + 273.15) + C.ch4 * log(T + 273.15) + D.ch4 * (T + 273.15)) / R)	
Cg.ch4   = WSP$Dis_CH4_1 * 10^-6  * pressure # vol atm ch4 / vol sample; Cg

Ca.ch4   = 55.5 * (Cg.ch4/H.ch4) * 16 * 10^3  # mg ch4/L H2O
Cah.ch4  = (V.headspace/V.aq * Cg.ch4 * (16/22.4) * (273.15/(T + 273.15))*10^3) # mg ch4/L H2O

ch4.aq   = (Ca.ch4 + Cah.ch4)  * 10^3 # ug ch4/L H2O
ch4.aq

WSP$Dis_CH4_cor <- ch4.aq

# Henry coefficients for co2

A.co2    = -317.658    
B.co2    =  1737.12     
C.co2    =  43.0607    
D.co2    = -0.00219107 

# HENRY's LAW
# H = 1/exp[{-a+b/(T+273)+C*ln(T+273)+d*(T+273)}]/R 
H.co2    = 1/exp((A.co2 + B.co2 / (T + 273.15) + C.co2 * log(T + 273.15) + D.co2 * (T + 273.15)) / R)   
Cg.co2   = WSP$Dis_CO2_1 * 10^-6 * pressure   # vol atm co2 / vol sample; Cg


Ca.co2   = 55.5 * (Cg.co2/H.co2) * 44 * 10^3  # mg co2/L H2O
Cah.co2  = (V.headspace/V.aq * Cg.co2 * (44/22.4) * (273.15/(T + 273.15))*10^3) # mg co2/L H2O

co2.aq   = (Ca.co2 + Cah.co2)  # mg co2/L H2O
co2.aq

WSP$Dis_CO2_cor <- co2.aq

#### Boxplot of all countinuous variables regarding ponds #### 

# This  is for make a list of graphs

plot_WSP <- function(data, column, column2){
    ggplot(data) +
        geom_boxplot(aes_string(x = column2 , y = column)) +
        xlab(column2) +
        ylab(column) 
}

plot_pond <- lapply(colnames(WSP)[7:27], plot_WSP, data = WSP, column2 = "Pond")
 
# print the graphs

lapply(plot_pond, print)

# save the graphs into folder

for (i in 1:length(plot_pond)){
    tiff(filename = paste("Boxplot_pond_",colnames(WSP)[7:27][i],".tiff", sep =""),units = 'px',height = 1800,width = 1800,res = 300,pointsize = 12)
    print(plot_pond[[i]])
    dev.off()
}

# Put graphs with similar variables together

ggsave("Boxplot_pond_1.tiff",marrangeGrob(plot_pond[1:4],nrow = 2, ncol= 2, top = ''),
       units = 'cm', height = 20, width = 20, dpi = 300)

ggsave("Boxplot_pond_2.tiff",marrangeGrob(plot_pond[5:8],nrow = 2, ncol= 2, top = ''),
       units = 'cm', height = 20, width = 20, dpi = 300)

ggsave("Boxplot_pond_3.tiff",marrangeGrob(plot_pond[9:11],nrow = 1, ncol= 3, top = ''),
       units = 'cm', height = 15, width = 30, dpi = 300)

ggsave("Boxplot_pond_4.tiff",marrangeGrob(plot_pond[12:17],nrow = 2, ncol= 3, top = ''),
       units = 'cm', height = 15, width = 30, dpi = 300)

ggsave("Boxplot_pond_5.tiff",marrangeGrob(plot_pond[18:21],nrow = 2, ncol= 2, top = ''),
       units = 'cm', height = 20, width = 20, dpi = 300)


#### Boxplot of all countinuous variables regarding ponds and position #### 

# This  is for make a list of graph

plot_WSP2 <- function (data, column, column2) {
    ggplot(data) +
        geom_boxplot(aes_string(x = column2 , y = column)) +
        xlab(column2) +
        ylab(column) +
        facet_wrap(.~Position)
}

plot_pond_position <- lapply(colnames(WSP)[7:27], plot_WSP2, data = WSP, column2 = "Pond")

# print the plot

lapply(plot_pond_position, print)

# save the plot into folder

# for (i in 1:length(plot_pond_position)){
#     tiff(filename = paste("Boxplot_pond_position_",colnames(WSP)[7:27][i],".tiff", sep =""),units = 'px',height = 1800,width = 1800,res = 300,pointsize = 12)
#     print(plot_pond_position[[i]])
#     dev.off()
# }

# Put graphs with similar variables together

ggsave("Boxplot_pond_position1.tiff",marrangeGrob(plot_pond_position[1:4],nrow = 2, ncol= 2, top = '')
       , units = 'cm', height = 20, width = 20,dpi = 300)

ggsave("Boxplot_pond_position2.tiff",marrangeGrob(plot_pond_position[5:8],nrow = 2, ncol= 2, top = '')
       , units = 'cm', height = 20, width = 20,dpi = 300)

ggsave("Boxplot_pond_position3.tiff",marrangeGrob(plot_pond_position[9:11],nrow = 1, ncol= 3, top = '')
       , units = 'cm', height = 15, width = 30,dpi = 300)

ggsave("Boxplot_pond_position4.tiff",marrangeGrob(plot_pond_position[12:17],nrow = 2, ncol= 3, top = '')
       , units = 'cm', height = 15, width = 30,dpi = 300)

ggsave("Boxplot_pond_position5.tiff",marrangeGrob(plot_pond_position[18:21],nrow = 2, ncol= 2, top = '')
       , units = 'cm', height = 20, width = 20,dpi = 300)


#### Boxplot of DG1 and DG2 regarding ponds ####

# transform to dissolved gas dataframe

WSP2 <- cbind(WSP[,4:6],WSP[,34:39])
WSP2 <- WSP2 %>% gather(key = 'Dissolved_gases', value = 'Concentration', - Pond, - Line, - Position)

# Plot them for the whole system 

ggsave("Boxplot_CO2.tiff", WSP2 %>% filter(str_detect(Dissolved_gases, "Dis_CO2_1|Dis_CO2_2")) %>% 
    ggplot() +
    geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
    xlab("Dissolved_gas") +
    ylab("CO2"))

ggsave("Boxplot_NO2.tiff", WSP2 %>% filter(str_detect(Dissolved_gases, "Dis_N2O_1|Dis_N2O_2")) %>% 
    ggplot() +
    geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
    xlab("Dissolved_gas") +
    ylab("NO2"))

ggsave("Boxplot_CH4.tiff", WSP2 %>% filter(str_detect(Dissolved_gases, "Dis_CH4_1|Dis_CH4_2")) %>% 
    ggplot() +
    geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
    xlab("Dissolved_gas") +
    ylab("CH4"))

# Plot them for different ponds

ggsave("Boxplot_pond_CO2.tiff", WSP2 %>% filter(str_detect(Dissolved_gases, "Dis_CO2_1|Dis_CO2_2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("CO2") +
           facet_wrap(.~Pond))

ggsave("Boxplot_pond_NO2.tiff", WSP2 %>% filter(str_detect(Dissolved_gases, "Dis_N2O_1|Dis_N2O_2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("NO2") +
           facet_wrap(.~Pond))

ggsave("Boxplot_pond_CH4.tiff", WSP2 %>% filter(str_detect(Dissolved_gases, "Dis_CH4_1|Dis_CH4_2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("CH4") +
           facet_wrap(.~Pond))


# Plot them for different ponds and lines

ggsave("Boxplot_pond_line_CO2.tiff", WSP2 %>% filter(str_detect(Dissolved_gases, "Dis_CO2_1|Dis_CO2_2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("CO2") +
           facet_grid(Line~Pond))
ggsave("Boxplot_pond_line_NO2.tiff", WSP2 %>% filter(str_detect(Dissolved_gases, "Dis_N2O_1|Dis_N2O_2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("NO2") +
           facet_grid(Line~Pond, scales = "free"))
ggsave("Boxplot_pond_line_CH4.tiff", WSP2 %>% filter(str_detect(Dissolved_gases, "Dis_CH4_1|Dis_CH4_2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("CH4") +
           facet_grid(Line~Pond))



#### Boxplot of DG_cor regarding ponds ####

# transform to dissolved gas dataframe

WSP2 <- cbind(WSP[,4:6],WSP[,44:46])
WSP2 <- WSP2 %>% gather(key = 'Dissolved_gases', value = 'Concentration', - Pond, - Line, - Position)

# Plot them for the whole system 

ggsave("Boxplot_CO2_cor.tiff", WSP2 %>% filter(str_detect(Dissolved_gases, "Dis_CO2_cor")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("CO2"))

ggsave("Boxplot_NO2_cor.tiff", WSP2 %>% filter(str_detect(Dissolved_gases, "Dis_N2O_cor")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("NO2"))

ggsave("Boxplot_CH4_cor.tiff", WSP2 %>% filter(str_detect(Dissolved_gases, "Dis_CH4_cor")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("CH4"))

# Plot them for different ponds

ggsave("Boxplot_pond_CO2_cor.tiff", WSP2 %>% filter(str_detect(Dissolved_gases, "Dis_CO2_cor")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("CO2") +
           facet_wrap(.~Pond))

ggsave("Boxplot_pond_NO2_cor.tiff", WSP2 %>% filter(str_detect(Dissolved_gases, "Dis_N2O_cor")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("NO2") +
           facet_wrap(.~Pond))

ggsave("Boxplot_pond_CH4_cor.tiff", WSP2 %>% filter(str_detect(Dissolved_gases, "Dis_CH4_cor")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("CH4") +
           facet_wrap(.~Pond))


# Plot them for different ponds and lines

ggsave("Boxplot_pond_line_CO2_cor.tiff", WSP2 %>% filter(str_detect(Dissolved_gases, "Dis_CO2_cor")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("CO2") +
           facet_grid(Line~Pond))
ggsave("Boxplot_pond_line_NO2_cor.tiff", WSP2 %>% filter(str_detect(Dissolved_gases, "Dis_N2O_cor")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("NO2") +
           facet_grid(Line~Pond, scales = "free"))
ggsave("Boxplot_pond_line_CH4_cor.tiff", WSP2 %>% filter(str_detect(Dissolved_gases, "Dis_CH4_cor")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Dissolved_gases), y = Concentration)) +
           xlab("Dissolved_gas") +
           ylab("CH4") +
           facet_grid(Line~Pond))



#### Boxplot of fluxes (mmol.m-2.d-1) regarding ponds ####

# transform to fluxes dataframe

WSP3 <- cbind(WSP[,4:6],WSP[,31:33])
colnames(WSP3)[4:6] <- c("CO2","N2O","CH4")
WSP3 <- WSP3 %>% gather(key = 'Fluxes', value = 'Concentration', - Pond, - Line, - Position)

WSP4 <- cbind(WSP[,4:6],WSP[,34:36])
colnames(WSP4)[4:6] <- c("CO2","N2O","CH4")
WSP4 <- WSP4 %>% gather(key = 'Fluxes', value = 'Concentration', - Pond, - Line, - Position)


# Plot them for the whole system 

ggsave("Boxplot_CO2_flux.tiff", WSP3 %>% filter(str_detect(Fluxes, "CO2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (mmol.m-2.d-1)") +
           ylab("CO2"))

ggsave("Boxplot_NO2_flux.tiff", WSP3 %>% filter(str_detect(Fluxes, "N2O")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (mmol.m-2.d-1)") +
           ylab("NO2"))

ggsave("Boxplot_CH4_flux.tiff", WSP3 %>% filter(str_detect(Fluxes, "CH4")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (mmol.m-2.d-1)") +
           ylab("CH4"))

# Plot them for different ponds

ggsave("Boxplot_pond_CO2_flux.tiff", WSP3 %>% filter(str_detect(Fluxes, "CO2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (mmol.m-2.d-1)") +
           ylab("CO2") +
           facet_wrap(.~Pond, scales = "free"))

ggsave("Boxplot_pond_NO2_flux.tiff", WSP3 %>% filter(str_detect(Fluxes, "N2O")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (mmol.m-2.d-1)") +
           ylab("NO2") +
           facet_wrap(.~Pond, scales = "free"))

ggsave("Boxplot_pond_CH4_flux.tiff", WSP3 %>% filter(str_detect(Fluxes, "CH4")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (mmol.m-2.d-1)") +
           ylab("CH4") +
           facet_wrap(.~Pond, scales = "free"))


# Plot them for different ponds and lines

ggsave("Boxplot_pond_line_CO2_flux.tiff", WSP3 %>% filter(str_detect(Fluxes, "CO2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (mmol.m-2.d-1)") +
           ylab("CO2") +
           facet_grid(Line~Pond, scales = "free"))
ggsave("Boxplot_pond_line_NO2_flux.tiff", WSP3 %>% filter(str_detect(Fluxes, "N2O")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (mmol.m-2.d-1)") +
           ylab("NO2") +
           facet_grid(Line~Pond, scales = "free"))
ggsave("Boxplot_pond_line_CH4_flux.tiff", WSP3 %>% filter(str_detect(Fluxes, "CH4")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (mmol.m-2.d-1)") +
           ylab("CH4") +
           facet_grid(Line~Pond, scales = "free"))



#### Boxplot of fluxes (g.m-2.d-1) regarding ponds ####

# transform to fluxes dataframe


WSP4 <- cbind(WSP[,4:6],WSP[,34:36])
colnames(WSP4)[4:6] <- c("CO2","N2O","CH4")
WSP4 <- WSP4 %>% gather(key = 'Fluxes', value = 'Concentration', - Pond, - Line, - Position)


# Plot them for the whole system 

ggsave("Boxplot_CO2_flux_2.tiff", WSP4 %>% filter(str_detect(Fluxes, "CO2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (g.m-2.d-1)") +
           ylab("CO2"))

ggsave("Boxplot_NO2_flux_2.tiff", WSP4 %>% filter(str_detect(Fluxes, "N2O")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (g.m-2.d-1)") +
           ylab("NO2"))

ggsave("Boxplot_CH4_flux_2.tiff", WSP4 %>% filter(str_detect(Fluxes, "CH4")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (g.m-2.d-1)") +
           ylab("CH4"))

# Plot them for different ponds

ggsave("Boxplot_pond_CO2_flux_2.tiff", WSP4 %>% filter(str_detect(Fluxes, "CO2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (g.m-2.d-1)") +
           ylab("CO2") +
           facet_wrap(.~Pond, scales = "free"))

ggsave("Boxplot_pond_NO2_flux_2.tiff", WSP4 %>% filter(str_detect(Fluxes, "N2O")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (g.m-2.d-1)") +
           ylab("NO2") +
           facet_wrap(.~Pond, scales = "free"))

ggsave("Boxplot_pond_CH4_flux_2.tiff", WSP4 %>% filter(str_detect(Fluxes, "CH4")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (g.m-2.d-1)") +
           ylab("CH4") +
           facet_wrap(.~Pond, scales = "free"))


# Plot them for different ponds and lines

ggsave("Boxplot_pond_line_CO2_flux_2.tiff", WSP4 %>% filter(str_detect(Fluxes, "CO2")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (g.m-2.d-1)") +
           ylab("CO2") +
           facet_grid(Line~Pond, scales = "free"))
ggsave("Boxplot_pond_line_NO2_flux_2.tiff", WSP4 %>% filter(str_detect(Fluxes, "N2O")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (g.m-2.d-1)") +
           ylab("NO2") +
           facet_grid(Line~Pond, scales = "free"))
ggsave("Boxplot_pond_line_CH4_flux_2.tiff", WSP4 %>% filter(str_detect(Fluxes, "CH4")) %>% 
           ggplot() +
           geom_boxplot(aes(x = as.factor(Fluxes), y = Concentration)) +
           xlab("Flux (g.m-2.d-1)") +
           ylab("CH4") +
           facet_grid(Line~Pond, scales = "free"))

#### Correlation coefficients #### 

ggpairs(WSP[7:27],
        lower = list(continuous = wrap("smooth", color = "deepskyblue")),
        upper = list(continuous = wrap("cor", size = 3, color = "tomato"))
) + theme(panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank())

#### Comparisons test (t-test, ANOVA, etc) #### 
# lack of data --> non parameteric analysis
