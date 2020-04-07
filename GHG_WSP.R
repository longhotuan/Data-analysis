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
library(lubridate)

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
library(ggcorrplot)
library(ggbiplot)
library(viridis)
library(wesanderson)
library(rvg)
library(officer)


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
library(rstatix)
library(robustHD)
library(MASS)
library(cluster)
library(ClusterR)
library(clusterSim)
library(clue)
library(parallel)
library(parallelMap)



#### import data #### 

WSP <- read_csv("WSP.csv")
WSP_24h <- read_csv("WSP_24h.csv")


#### 4C's process ####

WSP_24h$Line <- 2
WSP_24h$Position <- 2
WSP$Sample <- rep(seq(1:6),3)
WSP_24h$Pond <- str_replace_all(WSP_24h$Pond, "2","")


WSP$Pond <- as.factor(WSP$Pond)
WSP$Line <- as.factor(WSP$Line)
WSP$Position <- as.factor(WSP$Position)
WSP$Sample <- as.factor(WSP$Sample)

WSP_24h$Pond <- as.factor(WSP_24h$Pond)
WSP_24h$Line <- as.factor(WSP_24h$Line)
WSP_24h$Position <- as.factor(WSP_24h$Position)
WSP_24h$Sample <- as.factor(WSP_24h$Sample)

WSP$Date <- dmy(WSP$Date)
WSP$Time <- hms(WSP$Time)
WSP_24h$Date <- dmy(WSP_24h$Date)
WSP_24h$Time <- hms(WSP_24h$Time)


# FP2 position 3 with mean value

WSP$`BOD5 (mg/L)`[which(is.na(WSP$`BOD5 (mg/L)`))] <- mean(WSP$`BOD5 (mg/L)`[which(WSP$Pond == "FP" & WSP$Line == 2)], na.rm = TRUE)
WSP$`COD (mg/L)`[which(is.na(WSP$`COD (mg/L)`))] <- mean(WSP$`COD (mg/L)`[which(WSP$Pond == "FP" & WSP$Line == 2)], na.rm = TRUE)
WSP$`TN (mg/L)`[which(is.na(WSP$`TN (mg/L)`))] <- mean(WSP$`TN (mg/L)`[which(WSP$Pond == "FP" & WSP$Line == 2)], na.rm = TRUE)
WSP$`TP (mg/L)`[which(is.na(WSP$`TP (mg/L)`))] <- mean(WSP$`TP (mg/L)`[which(WSP$Pond == "FP" & WSP$Line == 2)], na.rm = TRUE)

# 24h datasets

WSP_24h$`BOD5 (mg/L)`[which(is.na(WSP_24h$`BOD5 (mg/L)`) & WSP_24h$Pond == "FP")] <- mean(WSP_24h$`BOD5 (mg/L)`[which(WSP_24h$Pond == "FP")], 
                                                                                           na.rm = TRUE)
WSP_24h$`COD (mg/L)`[which(is.na(WSP_24h$`COD (mg/L)`) & WSP_24h$Pond == "FP")] <- mean(WSP_24h$`COD (mg/L)`[which(WSP_24h$Pond == "FP")], 
                                                                                         na.rm = TRUE)
WSP_24h$`TN (mg/L)`[which(is.na(WSP_24h$`TN (mg/L)`) & WSP_24h$Pond == "FP")] <- mean(WSP_24h$`TN (mg/L)`[which(WSP_24h$Pond == "FP")], 
                                                                                       na.rm = TRUE)
WSP_24h$`TP (mg/L)`[which(is.na(WSP_24h$`TP (mg/L)`) & WSP_24h$Pond == "FP")] <- mean(WSP_24h$`TP (mg/L)`[which(WSP_24h$Pond == "FP")], 
                                                                                       na.rm = TRUE)
WSP_24h$`BOD5 (mg/L)`[which(is.na(WSP_24h$`BOD5 (mg/L)`) & WSP_24h$Pond == "MP")] <- mean(WSP_24h$`BOD5 (mg/L)`[which(WSP_24h$Pond == "MP")], 
                                                                                           na.rm = TRUE)

#### WSP_Correcting dissolved gas concentrations  ####

V.headspace = 6 # mL 
V.aq        = 6 # mL

R        = 1.98719   # gas constant in cal / K mol
T        = WSP$`Water temperature (oC)`  # assumed river temperature 
pressure = 2 # atm 

# Henry coefficients for N2O
A.n2o    = -180.95   
B.n2o    =  13205.8     
C.n2o    =  20.0399   
D.n2o    =  0.0238544

# HENRY's LAW
# H = 1/exp[{-a+b/(T+273)+C*ln(T+273)+d*(T+273)}]/R	
H.n2o    = 1/exp((A.n2o + B.n2o / (T + 273.15) + C.n2o * log(T + 273.15) + D.n2o * (T + 273.15)) / R)	

Cg.n2o   = WSP$`N2O dissolved gas (mg/L)` * 10^-6 * pressure  # vol atm N2O / vol sample; Cg


Ca.n2o   = 55.5 * (Cg.n2o/H.n2o) * 44 * 10^3  # mg N2O/L H2O
Cah.n2o  = (V.headspace/V.aq * Cg.n2o * (44/22.4) * (273.15/(T + 273.15))*10^3) # mg N2O/L H2O

n2o.aq   = (Ca.n2o + Cah.n2o) * 10^3 # ug N2O/L H2O
WSP$`N2O Dissolved Gas (ug/L)` <- n2o.aq

# Henry coefficients for ch4

A.ch4    = -365.183    
B.ch4    =  18106.7      
C.ch4    =  49.7554    
D.ch4    = -0.00028503 

# HENRY's LAW
# H = 1/exp[{-a+b/(T+273)+C*ln(T+273)+d*(T+273)}]/R	
H.ch4    = 1/exp((A.ch4 + B.ch4 / (T + 273.15) + C.ch4 * log(T + 273.15) + D.ch4 * (T + 273.15)) / R)	
Cg.ch4   = WSP$`CH4 dissolved gas (mg/L)` * 10^-6  * pressure # vol atm ch4 / vol sample; Cg

Ca.ch4   = 55.5 * (Cg.ch4/H.ch4) * 16 * 10^3  # mg ch4/L H2O
Cah.ch4  = (V.headspace/V.aq * Cg.ch4 * (16/22.4) * (273.15/(T + 273.15))*10^3) # mg ch4/L H2O

ch4.aq   = (Ca.ch4 + Cah.ch4)  * 10^3 # ug ch4/L H2O
ch4.aq

WSP$`CH4 Dissolved Gas (ug/L)` <- ch4.aq

# Henry coefficients for co2

A.co2    = -317.658    
B.co2    =  1737.12     
C.co2    =  43.0607    
D.co2    = -0.00219107 

# HENRY's LAW
# H = 1/exp[{-a+b/(T+273)+C*ln(T+273)+d*(T+273)}]/R 
H.co2    = 1/exp((A.co2 + B.co2 / (T + 273.15) + C.co2 * log(T + 273.15) + D.co2 * (T + 273.15)) / R)   
Cg.co2   = WSP$`CO2 dissolved gas (mg/L)` * 10^-6 * pressure   # vol atm co2 / vol sample; Cg


Ca.co2   = 55.5 * (Cg.co2/H.co2) * 44 * 10^3  # mg co2/L H2O
Cah.co2  = (V.headspace/V.aq * Cg.co2 * (44/22.4) * (273.15/(T + 273.15))*10^3) # mg co2/L H2O

co2.aq   = (Ca.co2 + Cah.co2)*10^3  # ug co2/L H2O
co2.aq

WSP$`CO2 Dissolved Gas (ug/L)` <- co2.aq

#### WSP_24h_Correcting dissolved gas concentrations  ####

V.headspace = 6 # mL 
V.aq        = 6 # mL

R        = 1.98719   # gas constant in cal / K mol
T        = WSP_24h$`Water temperature (oC)`  # assumed river temperature 
pressure = 2 # atm 

# Henry coefficients for N2O
A.n2o    = -180.95   
B.n2o    =  13205.8     
C.n2o    =  20.0399   
D.n2o    =  0.0238544

# HENRY's LAW
# H = 1/exp[{-a+b/(T+273)+C*ln(T+273)+d*(T+273)}]/R	
H.n2o    = 1/exp((A.n2o + B.n2o / (T + 273.15) + C.n2o * log(T + 273.15) + D.n2o * (T + 273.15)) / R)	

Cg.n2o   = WSP_24h$`N2O dissolved gas (mg/L)` * 10^-6 * pressure  # vol atm N2O / vol sample; Cg


Ca.n2o   = 55.5 * (Cg.n2o/H.n2o) * 44 * 10^3  # mg N2O/L H2O
Cah.n2o  = (V.headspace/V.aq * Cg.n2o * (44/22.4) * (273.15/(T + 273.15))*10^3) # mg N2O/L H2O

n2o.aq   = (Ca.n2o + Cah.n2o) * 10^3 # ug N2O/L H2O
WSP_24h$`N2O Dissolved Gas (ug/L)` <- n2o.aq

# Henry coefficients for ch4

A.ch4    = -365.183    
B.ch4    =  18106.7      
C.ch4    =  49.7554    
D.ch4    = -0.00028503 

# HENRY's LAW
# H = 1/exp[{-a+b/(T+273)+C*ln(T+273)+d*(T+273)}]/R	
H.ch4    = 1/exp((A.ch4 + B.ch4 / (T + 273.15) + C.ch4 * log(T + 273.15) + D.ch4 * (T + 273.15)) / R)	
Cg.ch4   = WSP_24h$`CH4 dissolved gas (mg/L)` * 10^-6  * pressure # vol atm ch4 / vol sample; Cg

Ca.ch4   = 55.5 * (Cg.ch4/H.ch4) * 16 * 10^3  # mg ch4/L H2O
Cah.ch4  = (V.headspace/V.aq * Cg.ch4 * (16/22.4) * (273.15/(T + 273.15))*10^3) # mg ch4/L H2O

ch4.aq   = (Ca.ch4 + Cah.ch4)  * 10^3 # ug ch4/L H2O
ch4.aq

WSP_24h$`CH4 Dissolved Gas (ug/L)` <- ch4.aq

# Henry coefficients for co2

A.co2    = -317.658    
B.co2    =  1737.12     
C.co2    =  43.0607    
D.co2    = -0.00219107 

# HENRY's LAW
# H = 1/exp[{-a+b/(T+273)+C*ln(T+273)+d*(T+273)}]/R 
H.co2    = 1/exp((A.co2 + B.co2 / (T + 273.15) + C.co2 * log(T + 273.15) + D.co2 * (T + 273.15)) / R)   
Cg.co2   = WSP_24h$`CO2 dissolved gas (mg/L)` * 10^-6 * pressure   # vol atm co2 / vol sample; Cg


Ca.co2   = 55.5 * (Cg.co2/H.co2) * 44 * 10^3  # mg co2/L H2O
Cah.co2  = (V.headspace/V.aq * Cg.co2 * (44/22.4) * (273.15/(T + 273.15))*10^3) # mg co2/L H2O

co2.aq   = (Ca.co2 + Cah.co2)*10^3  # ug co2/L H2O
co2.aq

WSP_24h$`CO2 Dissolved Gas (ug/L)` <- co2.aq

#### Combine two dataset ####
WSP_both <- bind_rows(WSP, WSP_24h)
WSP_both$Time <- c(WSP$Time, WSP_24h$Time)

WSP_both$Pond <- as.factor(WSP_both$Pond)
WSP_both$Line <- as.factor(WSP_both$Line)
WSP_both$Position <- as.factor(WSP_both$Position)
WSP_both$Sample <- as.factor(WSP_both$Sample)


#### Correct_Boxplot_Fluxes per ponds ####
#** WSP no separation ####

WSP_fluxes <- WSP %>% select(Pond, "Flux_CO2 (mg.m-2.d-1)","Flux_CH4 (mg.m-2.d-1)","Flux_NO2 (mg.m-2.d-1)") %>%
    pivot_longer(cols = -Pond, names_to = "GHGs", values_to = "Fluxes")
WSP_fluxes$GHGs <- as.factor(WSP_fluxes$GHGs)
WSP_fluxes$GHGs <- relevel(WSP_fluxes$GHGs,"Flux_CO2 (mg.m-2.d-1)")
WSP_fluxes$GHGs <- factor(WSP_fluxes$GHGs,
                          labels = c(expression("CO"["2"]), expression("CH"["4"]), 
                                     expression("N"["2"]*"O")))

ggsave("Fluxes_WSP.tiff", WSP_fluxes %>% ggplot(aes(x = Pond, y = Fluxes, fill = Pond)) +
           geom_boxplot() +
           stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
           theme_bw()+
           ylab(bquote("Fluxes (mg."*m^-2*"."*d^-1*")")) +
           facet_wrap(.~ GHGs, scales = "free", labeller = label_parsed) +
           scale_fill_brewer(palette = "Paired", name = "Ponds",
                             labels = c("Aerated Ponds", "Facultative Ponds", "Maturation Ponds"))+
           theme(text=element_text(size=14),
                 strip.text.x = element_text(size=14),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(),
                 legend.position="right",
                 legend.title = element_text(size = 14),
                 legend.text = element_text(size = 12),
                 legend.spacing.x = unit(0.5, 'cm')), 
       units = 'cm', height = 15, width = 30, dpi = 300
)

#** WSP with separation ####

WSP_fluxes_sep <- WSP %>% select("Flux_CO2 (mg.m-2.d-1)","Flux_CH4 (mg.m-2.d-1)","Flux_NO2 (mg.m-2.d-1)")
WSP_fluxes_sep$Pond <- paste(WSP$Pond, WSP$Line, sep = " ")
WSP_fluxes_sep <- WSP_fluxes_sep %>%  pivot_longer(cols = -Pond, names_to = "GHGs", values_to = "Fluxes")
WSP_fluxes_sep$Pond <- as.factor(WSP_fluxes_sep$Pond)
WSP_fluxes_sep$Pond <- factor(WSP_fluxes_sep$Pond,
                                 labels = c("Anaerobic Pond 1", "Anaerobic Pond 2", "Facultative Pond 1", 
                                            "Facultative Pond 2", "Maturation Pond 1", "Maturation Pond 2"))
WSP_fluxes_sep$GHGs <- as.factor(WSP_fluxes_sep$GHGs)
WSP_fluxes_sep$GHGs <- relevel(WSP_fluxes_sep$GHGs,"Flux_CO2 (mg.m-2.d-1)")
WSP_fluxes_sep$GHGs <- factor(WSP_fluxes_sep$GHGs,
                          labels = c(expression("CO"["2"]), expression("CH"["4"]), 
                                     expression("N"["2"]*"O")))

ggsave("Fluxes_WSP_separation.tiff", WSP_fluxes_sep %>% ggplot(aes(x = Pond, y = Fluxes, fill = Pond)) +
           geom_boxplot() +
           stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
           theme_bw()+
           ylab(bquote("Fluxes (mg."*m^-2*"."*d^-1*")")) +
           facet_wrap(.~ GHGs, scales = "free", labeller = label_parsed) +
           scale_fill_brewer(palette = "Paired", name = "Ponds")+
           theme(text=element_text(size=14),
                 strip.text.x = element_text(size=14),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(),
                 legend.position="right",
                 legend.title = element_text(size = 14),
                 legend.text = element_text(size = 12),
                 legend.spacing.x = unit(0.5, 'cm')), 
       units = 'cm', height = 15, width = 30, dpi = 300
)

#** WSP_24h ####

WSP_24h_fluxes <- WSP_24h %>% select(Pond, "Flux_CO2 (mg.m-2.d-1)","Flux_CH4 (mg.m-2.d-1)","Flux_NO2 (mg.m-2.d-1)") %>%
    pivot_longer(cols = -Pond, names_to = "GHGs", values_to = "Fluxes")
WSP_24h_fluxes$GHGs <- as.factor(WSP_24h_fluxes$GHGs)
WSP_24h_fluxes$GHGs <- relevel(WSP_24h_fluxes$GHGs,"Flux_CO2 (mg.m-2.d-1)")
WSP_24h_fluxes$GHGs <- factor(WSP_24h_fluxes$GHGs,
                          labels = c(expression("CO"["2"]), expression("CH"["4"]), 
                                     expression("N"["2"]*"O")))

ggsave("Fluxes_WSP_24h.tiff", WSP_24h_fluxes %>% ggplot(aes(x = Pond, y = Fluxes, fill = Pond)) +
           geom_boxplot() +
           stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red") +
           theme_bw()+
           ylab(bquote("Fluxes (mg."*m^-2*"."*d^-1*")")) +
           facet_wrap(.~ GHGs, scales = "free", labeller = label_parsed) +
           scale_fill_brewer(palette = "Paired", name = "Ponds",
                             labels = c("Facultative Pond 2", "Maturation Pond 2"))+
           theme(text=element_text(size=14),
                 strip.text.x = element_text(size=14),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(),
                 legend.position="right",
                 legend.title = element_text(size = 14),
                 legend.text = element_text(size = 12),
                 legend.spacing.x = unit(0.5, 'cm')), 
       units = 'cm', height = 15, width = 30, dpi = 300
)

#** WSP_both ####

WSP_both_sep <- WSP_both %>% select("Flux_CO2 (mg.m-2.d-1)","Flux_CH4 (mg.m-2.d-1)","Flux_NO2 (mg.m-2.d-1)")
WSP_both_sep$Pond <- paste(WSP_both$Pond, WSP_both$Line, sep = " ")
WSP_both_sep <- WSP_both_sep %>%  pivot_longer(cols = -Pond, names_to = "GHGs", values_to = "Fluxes")
WSP_both_sep$Pond <- as.factor(WSP_both_sep$Pond)
WSP_both_sep$Pond <- factor(WSP_both_sep$Pond,
                              labels = c("Anaerobic Pond 1", "Anaerobic Pond 2", "Facultative Pond 1", 
                                         "Facultative Pond 2", "Maturation Pond 1", "Maturation Pond 2"))
WSP_both_sep$GHGs <- as.factor(WSP_both_sep$GHGs)
WSP_both_sep$GHGs <- relevel(WSP_both_sep$GHGs,"Flux_CO2 (mg.m-2.d-1)")
WSP_both_sep$GHGs <- factor(WSP_both_sep$GHGs,
                              labels = c(expression("CO"["2"]), expression("CH"["4"]), 
                                         expression("N"["2"]*"O")))

ggsave("Fluxes_WSP_both_separation.tiff", WSP_both_sep %>% ggplot(aes(x = Pond, y = Fluxes, fill = Pond)) +
           geom_boxplot() +
           stat_summary(fun=mean, geom="point", shape=20, size=5, color="black", fill="black") +
           theme_bw()+
           ylab(bquote("Fluxes (mg."*m^-2*"."*d^-1*")")) +
           facet_wrap(.~ GHGs, scales = "free", labeller = label_parsed) +
           scale_fill_brewer(palette = "Paired", name = "Ponds")+
           theme(text=element_text(size=14),
                 strip.text.x = element_text(size=14),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(),
                 legend.position="right",
                 legend.title = element_text(size = 14),
                 legend.text = element_text(size = 12),
                 legend.spacing.x = unit(0.5, 'cm')), 
       units = 'cm', height = 15, width = 30, dpi = 300
)

#### Old_line_chart_WSP_24h ####


WSP_24h_line <- WSP_24h %>% select(Time, Pond, "Flux_CO2 (mg.m-2.d-1)","Flux_CH4 (mg.m-2.d-1)","Flux_NO2 (mg.m-2.d-1)") %>%
    pivot_longer(cols=c(-Pond, -Time),  names_to = "GHGs", values_to = "Fluxes")
WSP_24h_line$Pond <- as.factor(WSP_24h_line$Pond)
WSP_24h_line$Pond <- factor(WSP_24h_line$Pond,
                            labels = c("Facultative Pond", "Maturation Pond"))
WSP_24h_line$GHGs <- as.factor(WSP_24h_line$GHGs)
WSP_24h_line$GHGs <- relevel(WSP_24h_line$GHGs,"Flux_CO2 (mg.m-2.d-1)")
WSP_24h_line$GHGs <- factor(WSP_24h_line$GHGs,
                              labels = c(expression("CO"["2"]), expression("CH"["4"]), 
                                         expression("N"["2"]*"O")))

ggsave("Line_24h_CO2_fluxes.tiff", WSP_24h_line %>% 
           ggplot(aes(x = Time, y = Fluxes, group = GHGs)) +
           geom_line(aes(color = GHGs)) + 
           geom_point(aes(color = GHGs)) +
           theme_bw()+
           ylab(bquote("Fluxes (mg."*m^-2*"."*d^-1*")")) +
           facet_wrap(Pond ~ GHGs, scales = "free", labeller = labeller(GHGs = label_parsed))+
           scale_color_brewer(palette = "Dark2",)+
           theme(text=element_text(size=14),
                 strip.text.x = element_text(size=14),
                 axis.text.x = element_text(size=12),
                 axis.title.x = element_blank(),
                 legend.position = "none")
       , units = 'cm', height = 20, width = 50,dpi = 300)

#### Correct_summarise flux in each pond ####
#** WSP with separation ####

WSP_total <- WSP %>% select("Flux_CO2 (mg.m-2.d-1)","Flux_CH4 (mg.m-2.d-1)","Flux_NO2 (mg.m-2.d-1)")
WSP_total$Pond <- paste(WSP$Pond, WSP$Line, sep = " ")
WSP_total <- WSP_total %>% group_by(Pond) %>% summarise_each(funs(mean, std.error)) # mg.m-2.d-1
WSP_total$Area <- c(30000, 30000, 130000, 130000, 74000, 56000) # m2
colnames(WSP_total)[2:7] <- c("Mean_CO2", "Mean_CH4", "Mean_N2O", "SEM_CO2", "SEM_CH4", "SEM_N2O")
WSP_total <- WSP_total %>% mutate_at(vars(`Mean_CO2`, `Mean_CH4`, `Mean_N2O`), funs("area" = .*365*Area/1000000)) # kg.year-1
WSP_total <- WSP_total %>% mutate_at(vars("SEM_CO2", "SEM_CH4", "SEM_N2O"), funs("area" = .*365*Area/1000000)) # kg.year-1
WSP_total <- WSP_total %>% mutate_at(vars("Mean_CO2_area", "Mean_CH4_area", "Mean_N2O_area"), funs("percent" = .*100/ sum(.))) # percentage of area

WSP_total_stacked <- WSP_total[,c(1,15:16)] # negative N2O flux --> cannot calculate its total emissions
WSP_total_stacked <- WSP_total_stacked %>% pivot_longer(cols = -Pond, names_to = "GHGs", values_to = "Total emissions")
WSP_total_stacked$Pond <- as.factor(WSP_total_stacked$Pond)
WSP_total_stacked$Pond <- factor(WSP_total_stacked$Pond,
                                 labels = c("Anaerobic Pond 1", "Anaerobic Pond 2", "Facultative Pond 1", 
                                            "Facultative Pond 2", "Maturation Pond 1", "Maturation Pond 2"))
WSP_total_stacked$GHGs <- as.factor(WSP_total_stacked$GHGs)
WSP_total_stacked$GHGs <- relevel(WSP_total_stacked$GHGs, "Mean_CO2_area_percent")
WSP_total_stacked$GHGs <- factor(WSP_total_stacked$GHGs, 
                                 labels = c(expression("CO"["2"]), expression("CH"["4"])))

ggsave("Per_pond_total_emissions.tiff", WSP_total_stacked %>% 
         ggplot() +
         geom_bar(aes(y=`Total emissions`, x=GHGs,fill = Pond), stat = 'identity')+
         theme_bw() +
         # xlab("Year") +
         ylab("Fraction of the total emission per year (%)") +
         # facet_grid(.~Bank) +
         scale_fill_brewer(palette = "Paired") +
         scale_x_discrete(labels =c(bquote("CO"[2]), bquote("CH"[4])))+
         theme(text=element_text(size=14),
               strip.text.x = element_text(size=14),
               axis.text.x = element_text(size = 14),
               axis.ticks.x = element_blank(),
               axis.title.x = element_blank(),
               legend.position="right",
               legend.title = element_blank(),
               legend.text = element_text(size = 12),
               legend.spacing.x = unit(0.5, 'cm')),
     units = 'cm', height = 20, width = 20, dpi = 300)
#** WSP_both with separation ####

WSP_both_total <- WSP_both %>% select("Flux_CO2 (mg.m-2.d-1)","Flux_CH4 (mg.m-2.d-1)","Flux_NO2 (mg.m-2.d-1)")
WSP_both_total$Pond <- paste(WSP_both$Pond, WSP_both$Line, sep = " ")
WSP_both_total <- WSP_both_total %>% group_by(Pond) %>% summarise_each(funs(mean, std.error)) # mg.m-2.d-1
WSP_both_total$Area <- c(30000, 30000, 130000, 130000, 74000, 56000) # m2
colnames(WSP_both_total)[2:7] <- c("Mean_CO2", "Mean_CH4", "Mean_N2O", "SEM_CO2", "SEM_CH4", "SEM_N2O")
WSP_both_total <- WSP_both_total %>% mutate_at(vars(`Mean_CO2`, `Mean_CH4`, `Mean_N2O`), funs("area" = .*365*Area/1000000)) # kg.year-1
WSP_both_total <- WSP_both_total %>% mutate_at(vars("SEM_CO2", "SEM_CH4", "SEM_N2O"), funs("area" = .*365*Area/1000000)) # kg.year-1
WSP_both_total <- WSP_both_total %>% mutate_at(vars("Mean_CO2_area", "Mean_CH4_area", "Mean_N2O_area"), funs("percent" = .*100/ sum(.))) # percentage of area

WSP_both_total_stacked <- WSP_both_total[,c(1,15:16)] # negative N2O flux --> cannot calculate its total emissions
WSP_both_total_stacked <- WSP_both_total_stacked %>% pivot_longer(cols = -Pond, names_to = "GHGs", values_to = "Total emissions")
WSP_both_total_stacked$Pond <- as.factor(WSP_both_total_stacked$Pond)
WSP_both_total_stacked$Pond <- factor(WSP_both_total_stacked$Pond,
                                 labels = c("Anaerobic Pond 1", "Anaerobic Pond 2", "Facultative Pond 1", 
                                            "Facultative Pond 2", "Maturation Pond 1", "Maturation Pond 2"))
WSP_both_total_stacked$GHGs <- as.factor(WSP_both_total_stacked$GHGs)
WSP_both_total_stacked$GHGs <- relevel(WSP_both_total_stacked$GHGs, "Mean_CO2_area_percent")
WSP_both_total_stacked$GHGs <- factor(WSP_both_total_stacked$GHGs, 
                                 labels = c(expression("CO"["2"]), expression("CH"["4"])))

ggsave("Per_pond_total_emissions_both.tiff", WSP_both_total_stacked %>% 
           ggplot() +
           geom_bar(aes(y=`Total emissions`, x=GHGs,fill = Pond), stat = 'identity')+
           theme_bw() +
           # xlab("Year") +
           ylab("Fraction of the total emission per year (%)") +
           # facet_grid(.~Bank) +
           scale_fill_brewer(palette = "Paired") +
           scale_x_discrete(labels =c(bquote("CO"[2]), bquote("CH"[4])))+
           theme(text=element_text(size=14),
                 strip.text.x = element_text(size=14),
                 axis.text.x = element_text(size = 14),
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(),
                 legend.position="right",
                 legend.title = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.spacing.x = unit(0.5, 'cm')),
       units = 'cm', height = 20, width = 20, dpi = 300)


#### Correlation coefficients #### 

#** WSP ####

variable_WSP <- WSP[,c(7:30, 35:37)] %>% select(-`Rainfall (mm)`)

corr_WSP <- cor(variable_WSP, use = 'pairwise')
p.mat <- cor.mtest(variable_WSP)$p
colnames(p.mat) <- colnames(corr_WSP)
row.names(p.mat) <- colnames(corr_WSP)

# GGally Not really nice
ggsave("Corr_coeff.jpeg", ggpairs(variable_WSP,
                                  lower = list(continuous = wrap("smooth", color = "deepskyblue")),
                                  upper = list(continuous = wrap("cor", size = 3, color = "tomato"))) +
           theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()),
units = 'cm', height = 50, width = 50, dpi = 300)


# corrplot nicer but cannot handle the categorical variables

tiff("corr_coeff_2.tiff",units = 'cm',height = 50,width = 50,res = 300, pointsize = 12)
corrplot(corr_WSP, p.mat = p.mat, method = "circle", type = "upper",
         sig.level = 0.05, insig = "blank", order = "alphabet")
dev.off()

# use ggcorrplot 

ggsave("Corr_coeff_3.jpeg",ggcorrplot(corr_WSP,
                                      hc.order = TRUE,
                                      type = "lower",
                                      lab = TRUE),
       units = 'cm', height = 20, width = 20, dpi = 300)
ggsave("Corr_coeff_4.jpeg",ggcorrplot(corr_WSP,
                                      method="circle",
                                      hc.order = TRUE,
                                      type = "lower",
                                      lab = FALSE),
       units = 'cm', height = 20, width = 20, dpi = 300)

#** WSP_both ####

variable_WSP_both <- WSP_both[,c(7:30, 35:37)] %>% select(-`Rainfall (mm)`)

corr_WSP_both <- cor(variable_WSP_both, use = 'pairwise')
p.mat <- cor.mtest(variable_WSP_both)$p
colnames(p.mat) <- colnames(corr_WSP_both)
row.names(p.mat) <- colnames(corr_WSP_both)

# GGally Not really nice
ggsave("Corr_coeff_both.jpeg", ggpairs(variable_WSP_both,
                                  lower = list(continuous = wrap("smooth", color = "deepskyblue")),
                                  upper = list(continuous = wrap("cor", size = 3, color = "tomato"))) +
           theme(panel.grid.minor = element_blank(),
                 panel.grid.major = element_blank()),
       units = 'cm', height = 50, width = 50, dpi = 300)


# corrplot nicer but cannot handle the categorical variables

tiff("corr_coeff_2_both.tiff",units = 'cm',height = 50,width = 50,res = 300, pointsize = 12)
corrplot(corr_WSP_both, p.mat = p.mat, method = "circle", type = "upper",
         sig.level = 0.05, insig = "blank", order = "alphabet")
dev.off()

# use ggcorrplot 

ggsave("Corr_coeff_3_both.jpeg",ggcorrplot(corr_WSP_both,
                                      hc.order = TRUE,
                                      type = "lower",
                                      lab = TRUE),
       units = 'cm', height = 20, width = 20, dpi = 300)
ggsave("Corr_coeff_4_both.jpeg",ggcorrplot(corr_WSP_both,
                                      method="circle",
                                      hc.order = TRUE,
                                      type = "lower",
                                      lab = FALSE),
       units = 'cm', height = 20, width = 20, dpi = 300)


#### Friedmann test ####
#* CH4 ####
# were there differences or not? counting for heterogeneity samples unlike KW 
# cannot for date as it's not an unreplicated complete block design
#** Pond + Line variability 

WSP_fried <- WSP[,c( 4, 5,7:30, 35:37)]
WSP_fried_1 <- aggregate(data = WSP_fried, .~Line+Pond, mean)
WSP_fried_1 %>% friedman_test(`Flux_CH4 (mg.m-2.d-1)`~Pond|Line) #insignificant
WSP_fried_1 %>% friedman_test(`Flux_CH4 (mg.m-2.d-1)`~Line|Pond) #insignificant

#** Pond + Position variability 

WSP_fried <- WSP[,c( 4, 6,7:30, 35:37)]
WSP_fried_2 <- aggregate(data = WSP_fried, .~Position+Pond, mean)
WSP_fried_2 %>% friedman_test(`Flux_CH4 (mg.m-2.d-1)`~Pond|Position) # significant
WSP_fried_2 %>% friedman_test(`Flux_CH4 (mg.m-2.d-1)`~Position|Pond) # significant

#** Line + Position variability 

WSP_fried <- WSP[,c( 5, 6,7:30, 35:37)]
WSP_fried_3 <- aggregate(data = WSP_fried, .~Position+Line, mean)
WSP_fried_3 %>% friedman_test(`Flux_CH4 (mg.m-2.d-1)`~Line|Position) # insignificant
WSP_fried_3 %>% friedman_test(`Flux_CH4 (mg.m-2.d-1)`~Position|Line) # insignificant

#** Date + Position variability 

WSP_fried <- WSP[,c( 2, 6,7:30, 35:37)]
WSP_fried_4 <- aggregate(data = WSP_fried, .~Position+Date, mean)
WSP_fried_4 %>% friedman_test(`Flux_CH4 (mg.m-2.d-1)`~Date|Position) #insignificant
WSP_fried_4 %>% friedman_test(`Flux_CH4 (mg.m-2.d-1)`~Position|Date) #insignificant


#* CO2 ####
#** Pond + Line variability 

WSP_fried <- WSP[,c( 4, 5,7:30, 35:37)]
WSP_fried_1 <- aggregate(data = WSP_fried, .~Line+Pond, mean)
WSP_fried_1 %>% friedman_test(`Flux_CO2 (mg.m-2.d-1)`~Pond|Line) #insignificant
WSP_fried_1 %>% friedman_test(`Flux_CO2 (mg.m-2.d-1)`~Line|Pond) #insignificant

#** Pond + Position variability 

WSP_fried <- WSP[,c( 4, 6,7:30, 35:37)]
WSP_fried_2 <- aggregate(data = WSP_fried, .~Position+Pond, mean)
WSP_fried_2 %>% friedman_test(`Flux_CO2 (mg.m-2.d-1)`~Pond|Position) # significant
WSP_fried_2 %>% friedman_test(`Flux_CO2 (mg.m-2.d-1)`~Position|Pond) #insignificant

#** Line + Position variability 

WSP_fried <- WSP[,c( 5, 6,7:30, 35:37)]
WSP_fried_3 <- aggregate(data = WSP_fried, .~Position+Line, mean)
WSP_fried_3 %>% friedman_test(`Flux_CO2 (mg.m-2.d-1)`~Line|Position) # insignificant
WSP_fried_3 %>% friedman_test(`Flux_CO2 (mg.m-2.d-1)`~Position|Line) # insignificant

#** Date + Position variability 

WSP_fried <- WSP[,c( 2, 6,7:30, 35:37)]
WSP_fried_4 <- aggregate(data = WSP_fried, .~Position+Date, mean)
WSP_fried_4 %>% friedman_test(`Flux_CO2 (mg.m-2.d-1)`~Date|Position) #insignificant
WSP_fried_4 %>% friedman_test(`Flux_CO2 (mg.m-2.d-1)`~Position|Date) #insignificant


#* N2O ####
#** Pond + Line variability 

WSP_fried <- WSP[,c( 4, 5,7:30, 35:37)]
WSP_fried_1 <- aggregate(data = WSP_fried, .~Line+Pond, mean)
WSP_fried_1 %>% friedman_test(`Flux_NO2 (mg.m-2.d-1)`~Pond|Line) #insignificant
WSP_fried_1 %>% friedman_test(`Flux_NO2 (mg.m-2.d-1)`~Line|Pond) #insignificant

#** Pond + Position variability

WSP_fried <- WSP[,c( 4, 6,7:30, 35:37)]
WSP_fried_2 <- aggregate(data = WSP_fried, .~Position+Pond, mean)
WSP_fried_2 %>% friedman_test(`Flux_NO2 (mg.m-2.d-1)`~Pond|Position) # significant
WSP_fried_2 %>% friedman_test(`Flux_NO2 (mg.m-2.d-1)`~Position|Pond) #insignificant

#** Line + Position variability

WSP_fried <- WSP[,c( 5, 6,7:30, 35:37)]
WSP_fried_3 <- aggregate(data = WSP_fried, .~Position+Line, mean)
WSP_fried_3 %>% friedman_test(`Flux_NO2 (mg.m-2.d-1)`~Line|Position) # insignificant
WSP_fried_3 %>% friedman_test(`Flux_NO2 (mg.m-2.d-1)`~Position|Line) # insignificant

#** Date + Position variability 

WSP_fried <- WSP[,c( 2, 6,7:30, 35:37)]
WSP_fried_4 <- aggregate(data = WSP_fried, .~Position+Date, mean)
WSP_fried_4 %>% friedman_test(`Flux_NO2 (mg.m-2.d-1)`~Date|Position) #insignificant
WSP_fried_4 %>% friedman_test(`Flux_NO2 (mg.m-2.d-1)`~Position|Date) #insignificant

#### Wilcoxon signed-rank test ####
# identify which groups were different
#* CH4 ####
WSP %>% wilcox_test(`Flux_CH4 (mg.m-2.d-1)`~ Pond, paired = TRUE, p.adjust.method = "bonferroni")
WSP %>% wilcox_test(`Flux_CH4 (mg.m-2.d-1)`~ Line, paired = TRUE, p.adjust.method = "bonferroni")
WSP %>% wilcox_test(`Flux_CH4 (mg.m-2.d-1)`~ Position, paired = TRUE, p.adjust.method = "bonferroni")

#* CO2 ####
WSP %>% wilcox_test(`Flux_CO2 (mg.m-2.d-1)`~ Pond, paired = TRUE, p.adjust.method = "bonferroni")
WSP %>% wilcox_test(`Flux_CO2 (mg.m-2.d-1)`~ Line, paired = TRUE, p.adjust.method = "bonferroni")
WSP %>% wilcox_test(`Flux_CO2 (mg.m-2.d-1)`~ Position, paired = TRUE, p.adjust.method = "bonferroni")

#* N2O ####

WSP %>% wilcox_test(`Flux_NO2 (mg.m-2.d-1)`~ Pond, paired = TRUE, p.adjust.method = "bonferroni")
WSP %>% wilcox_test(`Flux_NO2 (mg.m-2.d-1)`~ Line, paired = TRUE, p.adjust.method = "bonferroni")
WSP %>% wilcox_test(`Flux_NO2 (mg.m-2.d-1)`~ Position, paired = TRUE, p.adjust.method = "bonferroni")

#### Mixed model_WSP ####
# try again with one random effect at a time
WSP_sta <- WSP

WSP_sta$log_N2O <- log(abs(WSP_sta$`Flux_NO2 (mg.m-2.d-1)`))
WSP_sta$sta_N2O <- standardize(WSP_sta$log_N2O) 
WSP_sta$log_CH4 <- log(WSP_sta$`Flux_CH4 (mg.m-2.d-1)`)
WSP_sta$sta_CH4 <- standardize(WSP_sta$log_CH4) 
WSP_sta$log_CO2 <- log(WSP_sta$`Flux_CO2 (mg.m-2.d-1)`)
WSP_sta$sta_CO2 <- standardize(WSP_sta$log_CO2) 
# Pond and Line-no significance ####
#* N2O #####
# diagnostic outliers 

# using box plot 

boxplot.stats(WSP_sta$sta_N2O)$out # no outliers

ggsave("Cleveland_N2O.tiff", ggplot(WSP_sta) +
           aes(x = sta_N2O, y = No) +
           geom_point(size = 3L, colour = "#0c4c8a") +
           xlab(bquote("Standardized Dissolved "*N[2]*"O")) +
           ylab("Order of the data")+
           theme_bw()+
           theme(axis.title.y = element_text(size = 14),
                 axis.title.x = element_text(size = 14),
                 text = element_text(size = 14)),
       units = 'cm', height = 20, width = 20, dpi = 300)

set.seed(1000)

WSP_lmm_N2O <- lmer(sta_N2O~1+(1|Pond/Line), data = WSP_sta)
r.squaredGLMM(WSP_lmm_N2O)
summary(WSP_lmm_N2O)
vcov(WSP_lmm_N2O)
var_N2O <- as.data.frame(VarCorr(WSP_lmm_N2O))
ICC_Line_N2O <- var_N2O$vcov[1]/sum(var_N2O$vcov)
ICC_WSP_N2O <- (var_N2O$vcov[1] + var_N2O$vcov[2])/sum(var_N2O$vcov)

#* CO2 #####
# diagnostic outliers 

# using box plot 

boxplot.stats(WSP_sta$sta_CO2)$out # no outliers

ggsave("Cleveland_CO2.tiff", ggplot(WSP_sta) +
           aes(x = sta_CO2, y = No) +
           geom_point(size = 3L, colour = "#0c4c8a") +
           xlab(bquote("Standardized Dissolved "*N[2]*"O")) +
           ylab("Order of the data")+
           theme_bw()+
           theme(axis.title.y = element_text(size = 14),
                 axis.title.x = element_text(size = 14),
                 text = element_text(size = 14)),
       units = 'cm', height = 20, width = 20, dpi = 300)

set.seed(1000)

WSP_lmm_CO2 <- lmer(sta_CO2~1+(1|Pond/Line), data = WSP_sta)
r.squaredGLMM(WSP_lmm_CO2)
summary(WSP_lmm_CO2)
vcov(WSP_lmm_CO2)
var_CO2 <- as.data.frame(VarCorr(WSP_lmm_CO2))
ICC_Line_CO2 <- var_CO2$vcov[1]/sum(var_CO2$vcov)
ICC_WSP_CO2 <- (var_CO2$vcov[1] + var_CO2$vcov[2])/sum(var_CO2$vcov)

#* CH4 #####
# diagnostic outliers 

# using box plot 

boxplot.stats(WSP_sta$sta_CH4)$out # no outliers

WSP

ggsave("Cleveland_CH4.tiff", ggplot(WSP_sta[WSP_sta$sta_CH4 < 2.46,]) +
           aes(x = sta_CH4, y = No) +
           geom_point(size = 3L, colour = "#0c4c8a") +
           xlab(bquote("Standardized Dissolved "*N[2]*"O")) +
           ylab("Order of the data")+
           theme_bw()+
           theme(axis.title.y = element_text(size = 14),
                 axis.title.x = element_text(size = 14),
                 text = element_text(size = 14)),
       units = 'cm', height = 20, width = 20, dpi = 300)

set.seed(1000)

WSP_lmm_CH4 <- lmer(sta_CH4~1+(1|Pond/Line), data = WSP_sta[WSP_sta$sta_CH4 < 2.46,])
r.squaredGLMM(WSP_lmm_CH4)
summary(WSP_lmm_CH4)
vcov(WSP_lmm_CH4)
var_CH4 <- as.data.frame(VarCorr(WSP_lmm_CH4))
ICC_Line_CH4 <- var_CH4$vcov[1]/sum(var_CH4$vcov)
ICC_WSP_CH4 <- (var_CH4$vcov[1] + var_CH4$vcov[2])/sum(var_CH4$vcov)

# Pond and Position-no significance ####
#* N2O #####

WSP_lmm_N2O <- lmer(sta_N2O~1+(1|Pond/Position), data = WSP_sta)
r.squaredGLMM(WSP_lmm_N2O)
summary(WSP_lmm_N2O)
vcov(WSP_lmm_N2O)
var_N2O <- as.data.frame(VarCorr(WSP_lmm_N2O))
ICC_Position_N2O <- var_N2O$vcov[1]/sum(var_N2O$vcov)
ICC_WSP2_N2O <- (var_N2O$vcov[1] + var_N2O$vcov[2])/sum(var_N2O$vcov)

#* CO2 #####

set.seed(1000)

WSP_lmm_CO2 <- lmer(sta_CO2~1+(1|Pond/Position), data = WSP_sta)
r.squaredGLMM(WSP_lmm_CO2)
summary(WSP_lmm_CO2)
vcov(WSP_lmm_CO2)
var_CO2 <- as.data.frame(VarCorr(WSP_lmm_CO2))
ICC_Position_CO2 <- var_CO2$vcov[1]/sum(var_CO2$vcov)
ICC_WSP2_CO2 <- (var_CO2$vcov[1] + var_CO2$vcov[2])/sum(var_CO2$vcov)

#* CH4 #####

set.seed(1000)

WSP_lmm_CH4 <- lmer(sta_CH4~1+(1|Pond/Position), data = WSP_sta[WSP_sta$sta_CH4 < 2.46,])
r.squaredGLMM(WSP_lmm_CH4)
summary(WSP_lmm_CH4)
vcov(WSP_lmm_CH4)
var_CH4 <- as.data.frame(VarCorr(WSP_lmm_CH4))
ICC_Position_CH4 <- var_CH4$vcov[1]/sum(var_CH4$vcov)
ICC_WSP2_CH4 <- (var_CH4$vcov[1] + var_CH4$vcov[2])/sum(var_CH4$vcov)

#### PCA_WSP_both ####
# use prcomp in the package MASS

WSP_both_PCA <- WSP_both[,c(7:8, 10, 12, 15:27)]
WSP_labels <- str_split_fixed(colnames(WSP_both_PCA), " \\(", n = 2)[,1]
WSP_labels <- str_replace_all(WSP_labels, " ", "_")
WSP_labels <- str_replace_all(WSP_labels, "-", "")
WSP_labels <- str_replace_all(WSP_labels, "\\+", "")
colnames(WSP_both_PCA) <- WSP_labels

WSP_cor <- cor(WSP_both_PCA, use = "na.or.complete")

WSP_both_PCA_2 <- prcomp(WSP_both_PCA, scale. = TRUE)

WSP_both_PCA_2
summary(WSP_both_PCA_2)

# Eigenvectors (= matrix U)
WSP_both_PCA_2$rotation

# Component matrix (= matrix A)
WSP_both_PCA_2$rotation %*% diag(WSP_both_PCA_2$sdev)

#' Component scores (= matrix C)
WSP_both_PCA_2$x
head(WSP_both_PCA_2$x)

# Scree plot
ggscreeplot(WSP_both_PCA_2) + geom_col()

# Biplot
WSP_CO2 <- WSP_both$`Flux_CO2`
WSP_CH4 <- WSP_both$`Flux_CH4`
WSP_Pond <- as.factor(WSP_both$Pond)

#** CO2 ####

ggbiplot(WSP_both_PCA_2, obs.scale = 1, var.scale = 1, alpha=0, varname.size = 4, labels.size= 6) +
    geom_point(aes(color = WSP_CO2 , shape = WSP_Pond, size = WSP_CO2)) +
    scale_color_gradient(low = "blue", high = "red") + 
    theme_classic() +
    scale_x_continuous(limits=c(-5,5)) +
    scale_y_continuous(limits=c(-5,5)) +
    theme(text=element_text(size=14),
          legend.position="right",
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.spacing.x = unit(0.5, 'cm'))

#** CH4 ####

ggbiplot(WSP_both_PCA_2, obs.scale = 1, var.scale = 1, alpha=0, varname.size = 4, labels.size= 6) +
    geom_point(aes(color = WSP_CH4 , shape = WSP_Pond, size = WSP_CH4)) +
    scale_color_gradient(low = "blue", high = "red") + 
    theme_classic() +
    scale_x_continuous(limits=c(-5,5)) +
    scale_y_continuous(limits=c(-5,5)) +
    theme(text=element_text(size=14),
          legend.position="right",
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.spacing.x = unit(0.5, 'cm'))

#### correct_K means_stats clusering (no hyperparameter tuning) ####
#** correct_using K means to classify PC1 and PC2 ####

WSP_both_kmean <- WSP_both_PCA_2$x[,1:2]

# Using the elbow method to find the optimal number of clusters

set.seed(1)
wcss <- vector()
for (i in 1:10) wcss[i] <- sum(kmeans(WSP_both_kmean, i)$withinss)
plot(1:10,
     wcss,
     type = 'b',
     main = paste('The Elbow Method'),
     xlab = 'Number of clusters',
     ylab = 'WCSS')

# Fitting K-Means to the dataset

set.seed(10)
WSP_both_kmeans <- kmeans(x = WSP_both_kmean, centers = 4)
WSP_both_y_kmeans <- WSP_both_kmeans$cluster

# Visualising the clusters
clusplot(WSP_both_kmean,
         WSP_both_y_kmeans,
         lines = 0,
         shade = TRUE,
         color = TRUE,
         labels = 2,
         plotchar = FALSE,
         span = TRUE,
         main = paste('Clusters of customers'),
         xlab = 'PC1',
         ylab = 'PC2')

#** CO2 ####

WSP_both_PCA_kmeans_CO2 <- ggbiplot(WSP_both_PCA_2, obs.scale = 1, var.scale = 1, alpha=0, varname.size =4, labels.size= 6, circle = TRUE) +
    geom_point(aes(color =  as.factor(WSP_both_kmeans$cluster) , shape = WSP_Pond, size = WSP_CO2)) +
    scale_color_brewer(palette = "Dark2") + 
    theme_classic() +
    ylab("PC2 (19.3%)")+
    xlab("PC1 (32.0%)") + 
    scale_x_continuous(limits=c(-5,5)) +
    scale_y_continuous(limits=c(-5,5)) +
    theme(text=element_text(size=14),
          legend.position="right",
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.spacing.x = unit(0.5, 'cm'))
WSP_both_PCA_kmeans_CO2
ggsave("PCA_Kmeans_CO2.jpg", WSP_both_PCA_kmeans,
       units = 'cm', height = 20, width = 20, dpi = 300)

WSP_both_PCA_kmeans_editable_CO2 <- dml(ggobj = WSP_both_PCA_kmeans_CO2)
WSP_both_PCA_kmeans_doc <- read_pptx()
WSP_both_PCA_kmeans_doc <- add_slide(WSP_both_PCA_kmeans_doc)
WSP_both_PCA_kmeans_doc <- ph_with(x = WSP_both_PCA_kmeans_doc, WSP_both_PCA_kmeans_editable_CO2,
               location = ph_location_type(type = "body") )
print(WSP_both_PCA_kmeans_doc, target = "WSP_both_PCA_kmeans_CO2.pptx")

#** CH4 ####

WSP_both_PCA_kmeans_CH4 <- ggbiplot(WSP_both_PCA_2, obs.scale = 1, var.scale = 1, alpha=0, varname.size = 4, labels.size= 6) +
    geom_point(aes(color =  as.factor(WSP_both_kmeans$cluster) , shape = WSP_Pond, size = WSP_CH4)) +
    scale_color_brewer(palette = "Dark2") + 
    theme_classic() +
    scale_x_continuous(limits=c(-5,5)) +
    scale_y_continuous(limits=c(-5,5)) +
    ylab("PC2 (19.3%)")+
    xlab("PC1 (32.0%)") + 
    theme(text=element_text(size=14),
          legend.position="right",
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.spacing.x = unit(0.5, 'cm'))
WSP_both_PCA_kmeans_CH4
ggsave("PCA_Kmeans_CH4.jpg", WSP_both_PCA_kmeans_CH4,
       units = 'cm', height = 20, width = 20, dpi = 300)

WSP_both_PCA_kmeans_editable_CH4 <- dml(ggobj = WSP_both_PCA_kmeans_CH4)
WSP_both_PCA_kmeans_doc <- read_pptx()
WSP_both_PCA_kmeans_doc <- add_slide(WSP_both_PCA_kmeans_doc)
WSP_both_PCA_kmeans_doc <- ph_with(x = WSP_both_PCA_kmeans_doc, WSP_both_PCA_kmeans_editable_CH4,
                                   location = ph_location_type(type = "body") )
print(WSP_both_PCA_kmeans_doc, target = "WSP_both_PCA_kmeans_CH4.pptx")


#** wrong_using K means to classify CO2 and CH4 ####

WSP_both_kmean <- WSP_both %>% select(`Flux_CO2 (mg.m-2.d-1)`,`Flux_CH4 (mg.m-2.d-1)`)

# Using the elbow method to find the optimal number of clusters

set.seed(1)
wcss <- vector()
for (i in 1:10) wcss[i] <- sum(kmeans(WSP_both_kmean, i)$withinss)
plot(1:10,
     wcss,
     type = 'b',
     main = paste('The Elbow Method'),
     xlab = 'Number of clusters',
     ylab = 'WCSS')

# Fitting K-Means to the dataset

set.seed(10)
WSP_both_kmeans <- kmeans(x = WSP_both_kmean, centers = 4, iter.max = 1000, nstart = 10)
WSP_both_y_kmeans <- WSP_both_kmeans$cluster

# Visualising the clusters
clusplot(WSP_both_kmean,
         WSP_both_y_kmeans,
         lines = 0,
         shade = TRUE,
         color = TRUE,
         labels = 2,
         plotchar = FALSE,
         span = TRUE,
         main = paste('Clusters of customers'),
         xlab = 'PC1',
         ylab = 'PC2')

labels_WSP_kmeans <- paste(WSP_both$Pond)

#** CO2 ####
WSP_both_PCA_kmeans_CO2_v2 <- ggbiplot(WSP_both_PCA_2, obs.scale = 1, var.scale = 1, alpha=0, varname.size =0, labels.size= 6) +
    geom_point(aes(color =  as.factor(WSP_both_kmeans$cluster) , shape = WSP_Pond, size = WSP_CO2)) +
    scale_color_brewer(palette = "Dark2") + 
    theme_classic() +
    ylab("PC2 (19.2%)")+
    xlab("PC1 (31.8%)") + 
    scale_x_continuous(limits=c(-5,5)) +
    scale_y_continuous(limits=c(-5,5)) +
    theme(text=element_text(size=14),
          legend.position="right",
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.spacing.x = unit(0.5, 'cm'))

ggsave("PCA_Kmeans_CO2_v2.jpeg", WSP_both_PCA_kmeans_CO2_v2,
       units = 'cm', height = 20, width = 20, dpi = 300)

WSP_both_PCA_kmeans_editable_CO2_v2 <- dml(ggobj = WSP_both_PCA_kmeans_CO2_v2)
WSP_both_PCA_kmeans_doc <- read_pptx()
WSP_both_PCA_kmeans_doc <- add_slide(WSP_both_PCA_kmeans_doc)
WSP_both_PCA_kmeans_doc <- ph_with(x = WSP_both_PCA_kmeans_doc, WSP_both_PCA_kmeans_editable_CO2_v2,
                                   location = ph_location_type(type = "body") )
print(WSP_both_PCA_kmeans_doc, target = "WSP_both_PCA_kmeans_CO2_v2.pptx")

#** CH4 ####

WSP_both_PCA_kmeans_CH4_v2 <- ggbiplot(WSP_both_PCA_2, obs.scale = 1, var.scale = 1, alpha=0, varname.size = 4, labels.size= 6) +
    geom_point(aes(color =  as.factor(WSP_both_kmeans$cluster) , shape = WSP_Pond, size = WSP_CH4)) +
    scale_color_brewer(palette = "Dark2") + 
    theme_classic() +
    scale_x_continuous(limits=c(-5,5)) +
    scale_y_continuous(limits=c(-5,5)) +
    ylab("PC2 (19.2%)")+
    xlab("PC1 (31.8%)") + 
    theme(text=element_text(size=14),
          legend.position="right",
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.spacing.x = unit(0.5, 'cm'))

ggsave("PCA_Kmeans_CH4_v2.jpeg", WSP_both_PCA_kmeans_CH4,
       units = 'cm', height = 20, width = 20, dpi = 300)

WSP_both_PCA_kmeans_editable_CH4_v2 <- dml(ggobj = WSP_both_PCA_kmeans_CH4_v2)
WSP_both_PCA_kmeans_doc <- read_pptx()
WSP_both_PCA_kmeans_doc <- add_slide(WSP_both_PCA_kmeans_doc)
WSP_both_PCA_kmeans_doc <- ph_with(x = WSP_both_PCA_kmeans_doc, WSP_both_PCA_kmeans_editable_CH4_v2,
                                   location = ph_location_type(type = "body"))
print(WSP_both_PCA_kmeans_doc, target = "WSP_both_PCA_kmeans_CH4_v2.pptx")

#### correct_RF (hyperparameter tuning) ####
#** CO2 ####

WSP_both_RF <- WSP_both[,c(7, 8, 10, 12, 15:28)]
WSP_labels <- str_split_fixed(colnames(WSP_both_RF), " \\(", n = 2)[,1]
WSP_labels <- str_replace_all(WSP_labels, " ", "_")
WSP_labels <- str_replace_all(WSP_labels, "-", "")
WSP_labels <- str_replace_all(WSP_labels, "\\+", "")
colnames(WSP_both_RF) <- WSP_labels

trainTask_CO2_2 <- makeRegrTask(data = WSP_both_RF,target = "Flux_CO2")

# create mlr learner
set.seed(1234)
lrn_CO2_2 <- makeLearner("regr.ranger")
# lrn_CO2_2$par.vals <- list(ntree = 100L, importance = TRUE)
cv_CO2_2 <- makeResampleDesc(method = "LOO")
# set parallel backend
parallelStartSocket(cpus = detectCores()-1)
res_CO2_2 <- resample(learner = lrn_CO2_2, task = trainTask_CO2_2, resampling = cv_CO2_2)
# Tuning hyperparameters

# Parameter Tuning: Mainly, there are three parameters in the random forest algorithm which you should look at (for tuning):
# ntree - The number of trees to grow. Larger the tree, it will be more computationally expensive to build models.
# mtry - It refers to how many variables we should select at a node split. 
# Also as mentioned above, the default value is p/3 for regression and sqrt(p) for classification. 
# We should always try to avoid using smaller values of mtry to avoid overfitting.
# nodesize - It refers to how many observations we want in the terminal nodes.
# This parameter is directly related to tree depth. Higher the number, lower the tree depth. 
# With lower tree depth, the tree might even fail to recognize useful signals from the data.

# To know which hyperparameter can be tuned using getParamSet
getParamSet(lrn_CO2_2)

# Only tune three abovementioned hyperparameters

params_CO2_2 <- makeParamSet(makeIntegerParam("mtry", lower = 2, upper = 10),
                             makeIntegerParam("min.node.size", lower = 2, upper = 25))
tc_CO2_2 <- makeTuneControlMBO(budget = 100)
tr_CO2_2 <- tuneParams(learner = lrn_CO2_2, task = trainTask_CO2_2, resampling = cv_CO2_2,
                       par.set = params_CO2_2, control = tc_CO2_2, show.info = TRUE)

parallelStop()
tr_CO2_2
# Tune result:
# Op. pars: mtry=6; min.node.size=2
# mse.test.mean=34718792.6322091
# Apply the optimal RF

# using ranger package with permutation feature importance
set.seed(1234)
regressor_CO2_2_opt_ranger <- ranger(formula = Flux_CO2~., data = WSP_both_RF, num.trees =  1000, mtry = 6, 
                                     importance = 'permutation', min.node.size = 2) 
imp_CO2_2_opt_per <- importance_pvalues(x = regressor_CO2_2_opt_ranger, method = 'janitza', num.permutations = 100)

featureImportance_CO2_2_opt_per <- data.frame(Feature=row.names(imp_CO2_2_opt_per), Importance=imp_CO2_2_opt_per[,1], 
                                              pvalue=imp_CO2_2_opt_per[,2])
# remove the variables with pvalue >0.05
featureImportance_CO2_2_opt_per <- featureImportance_CO2_2_opt_per %>% filter(pvalue <=0.05)

labels_CO2_2_per <- c("DO", "Air ~ temperature", "pH", "Total ~ Dissolved ~ Solids", "PO[4]^3^-{}", 
                      "BOD[5]", "Solar ~ radiation","Water ~ temperature", "Chlorophyll*~alpha", "Wind")
labels_CO2_2_per <- rev(labels_CO2_2_per)
labels_CO2_2_per_parse <- parse(text = labels_CO2_2_per)

ggsave("RF_CO2_2_opt_per.tiff", ggplot(featureImportance_CO2_2_opt_per[1:10,], aes(x=reorder(Feature, Importance), y=Importance)) +
           geom_bar(stat="identity", fill="tomato") +
           coord_flip() + 
           theme_bw(base_size=20) +
           scale_x_discrete(labels = labels_CO2_2_per_parse) +
           scale_y_continuous(labels = seq(from = 0.1, to = 0.4, length.out = 4))+
           labs(x = NULL, y = "Scale Importance") + 
           ggtitle(bquote("C"*O[2]*"")) +
           theme(plot.title=element_text(size=18)),
       units = 'cm', height = 20, width = 20, dpi = 300)

#** CH4 ####

WSP_both_RF <- WSP_both[,c(7, 8, 10, 12, 15:27, 29)]
WSP_both_RF <- WSP_both_RF[-5,]
WSP_labels <- str_split_fixed(colnames(WSP_both_RF), " \\(", n = 2)[,1]
WSP_labels <- str_replace_all(WSP_labels, " ", "_")
WSP_labels <- str_replace_all(WSP_labels, "-", "")
WSP_labels <- str_replace_all(WSP_labels, "\\+", "")
colnames(WSP_both_RF) <- WSP_labels

trainTask_CH4_2 <- makeRegrTask(data = WSP_both_RF,target = "Flux_CH4")

# create mlr learner
set.seed(1234)
lrn_CH4_2 <- makeLearner("regr.ranger")
# lrn_CH4_2$par.vals <- list(ntree = 100L, importance = TRUE)
cv_CH4_2 <- makeResampleDesc(method = "LOO")
# set parallel backend
parallelStartSocket(cpus = detectCores()-1)
res_CH4_2 <- resample(learner = lrn_CH4_2, task = trainTask_CH4_2, resampling = cv_CH4_2)
# Tuning hyperparameters

# To know which hyperparameter can be tuned using getParamSet
getParamSet(lrn_CH4_2)

# Only tune three abovementioned hyperparameters

params_CH4_2 <- makeParamSet(makeIntegerParam("mtry", lower = 2, upper = 10),
                             makeIntegerParam("min.node.size", lower = 2, upper = 25))
tc_CH4_2 <- makeTuneControlMBO(budget = 100)
tr_CH4_2 <- tuneParams(learner = lrn_CH4_2, task = trainTask_CH4_2, resampling = cv_CH4_2,
                       par.set = params_CH4_2, control = tc_CH4_2, show.info = TRUE)

parallelStop()
tr_CH4_2
# Tune result:
# Op. pars: mtry=10; min.node.size=5
# mse.test.mean=16952445.4909653
# Apply the optimal RF

# using ranger package with permutation feature importance
set.seed(1234)
regressor_CH4_2_opt_ranger <- ranger(formula = Flux_CH4~., data = WSP_both_RF, num.trees =  1000, mtry = 9, 
                                     importance = 'permutation', min.node.size = 4) 
imp_CH4_2_opt_per <- importance_pvalues(x = regressor_CH4_2_opt_ranger, method = 'janitza', num.permutations = 100)

featureImportance_CH4_2_opt_per <- data.frame(Feature=row.names(imp_CH4_2_opt_per), Importance=imp_CH4_2_opt_per[,1], 
                                              pvalue=imp_CH4_2_opt_per[,2])
# remove the variables with pvalue >0.05
featureImportance_CH4_2_opt_per <- featureImportance_CH4_2_opt_per %>% filter(pvalue <=0.05)

labels_CH4_2_per <- c("BOD[5]", "DO", "Total ~ Dissolved ~ Solids", "Air ~ temperature", "pH", "PO[4]^3^-{}", "COD","Water ~ temperature", "TN", 
                      "TP")
labels_CH4_2_per <- rev(labels_CH4_2_per)
labels_CH4_2_per_parse <- parse(text = labels_CH4_2_per)

ggsave("RF_CH4_2_opt_per.tiff", ggplot(featureImportance_CH4_2_opt_per[1:10,], aes(x=reorder(Feature, Importance), y=Importance)) +
           geom_bar(stat="identity", fill="tomato") +
           coord_flip() + 
           theme_bw(base_size=20) +
           scale_x_discrete(labels = labels_CH4_2_per_parse) +
           scale_y_continuous(labels = seq(from = 0.1, to = 0.5, length.out = 5))+
           labs(x = NULL, y = "Scale Importance") +
           ggtitle(bquote("C"*H[4]*"")) +
           theme(plot.title=element_text(size=18)),
       units = 'cm', height = 20, width = 20, dpi = 300)
