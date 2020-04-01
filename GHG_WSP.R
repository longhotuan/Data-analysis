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
#** WSP no separation between lines ####

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

#** WSP with separation between lines ####

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

WSP_total <- WSP[,c(28:30)] 
WSP_total$Pond <- paste(WSP$Pond, WSP$Line, sep = " ")
WSP_total <- WSP_total %>% group_by(Pond) %>% summarise_each(funs(mean, std.error)) # mg.m-2.d-1
WSP_total$Area <- c(30000, 30000, 130000, 130000, 74000, 56000) # m2
colnames(WSP_total)[2:7] <- c("Mean_CO2", "Mean_CH4", "Mean_N2O", "SEM_CO2", "SEM_CH4", "SEM_N2O")
WSP_total <- WSP_total %>% mutate_at(vars(`Mean_CO2`, `Mean_CH4`, `Mean_N2O`), funs("area" = .*365*Area/1000000)) # kg.year-1
WSP_total <- WSP_total %>% mutate_at(vars("SEM_CO2", "SEM_CH4", "SEM_N2O"), funs("area" = .*365*Area/1000000)) # kg.year-1
WSP_total <- WSP_total %>% mutate_at(vars("Mean_CO2_area", "Mean_CH4_area", "Mean_N2O_area"), funs("percent" = .*100/ sum(.))) # percentage of area

# Separate lines 1 and 2 

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


#### Correlation coefficients #### 
variable_river <- cbind(river[,6:13], river[,22:33]) %>% select(-Rain)

corr_river <- cor(variable_river, use = 'pairwise')
p.mat <- cor.mtest(variable_river)$p
colnames(p.mat) <- colnames(corr_river)
row.names(p.mat) <- colnames(corr_river)

# GGally Not really nice
ggsave("Corr_coeff.tiff", ggpairs(variable_river,
                                  lower = list(continuous = wrap("smooth", color = "deepskyblue")),
                                  upper = list(continuous = wrap("cor", size = 3, color = "tomato"))
) + theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()),
units = 'cm', height = 50, width = 50, dpi = 300)


# corrplot nicer but cannot handle the categorical variables

tiff("corr_coeff_2.tiff",units = 'cm',height = 50,width = 50,res = 300, pointsize = 12)
corrplot(corr_river, p.mat = p.mat, method = "circle", type = "upper",
         sig.level = 0.05, insig = "blank", order = "alphabet")
dev.off()

# Using mosaic plot to represent the relationship among two or more categorical variables 
# in this case, only for river and hydromorphological data

ggsave("Mosaic_river_LB.tiff", ggplot(river)+
           geom_mosaic(aes(x= product(River), fill = "Lelf Bank"))+
           labs(x ="", y = ""),
       units = 'cm', height = 30, width = 40, dpi = 300)
ggsave("Mosaic_river_RB.tiff", ggplot(river)+
           geom_mosaic(aes(x= product(River), fill = "Right Bank"))+
           labs(x ="", y = ""),
       units = 'cm', height = 30, width = 40, dpi = 300)
ggsave("Mosaic_river_FV.tiff",ggplot(river)+
           geom_mosaic(aes(x= product(River), fill = "Flow variability"))+
           labs(x ="", y = ""),
       units = 'cm', height = 30, width = 40, dpi = 300)
ggsave("Mosaic_river_shading.tiff",ggplot(river)+
           geom_mosaic(aes(x= product(River), fill = Shading))+
           labs(x ="", y = ""),
       units = 'cm', height = 30, width = 40, dpi = 300)
ggsave("Mosaic_river_erosion.tiff",ggplot(river)+
           geom_mosaic(aes(x= product(River), fill = Erosion))+
           labs(x ="", y = ""),
       units = 'cm', height = 30, width = 40, dpi = 300)
ggsave("Mosaic_river_pool.tiff",ggplot(river)+
           geom_mosaic(aes(x= product(River), fill = "Pool Class"))+
           labs(x ="", y = ""),
       units = 'cm', height = 30, width = 40, dpi = 300)


#### Spatio-temporal variability ####
#*** Friedmann test ####



river_fried <- river %>% select(c(2, 5, 36:38))
river_fried_1 <- aggregate(data = river_fried, .~Date + River, mean)



#### Wrong_Permutation testing  ####
# lack of data --> non parameteric analysis (not sure about the distribution of the data)
# --> using Permanova for multivariate comparison for testing the simultaneous response of one or more variables to one or more
# factors in an ANOVA experimental design on the basis of any distance measure, using permutation methods
# to accommodate random effects, hierarchical models, mixed models, quantitative covariates,
# repeated measures, unbalanced and/or asymmetrical designs, and, most recently, heterogeneous dispersions among groups.
# or Fried.mann for univariate comparison repeated measures with block effects to avoid dependent samples.
# Test the multivariate homogeneity of groups dispersions

mod <- betadisper(daisy(GHGes, metric = "euclidean", stand = TRUE), group = river$River) # using betadisper is a multivariate analogue of Levene's test for homogeneity of variances.
permutest(mod)
anova(mod)
plot(mod, hull=FALSE, ellipse=TRUE)
boxplot(mod)

# p value > 0.05 --> homogeneity of multivariate dispersions.
# PERMANOVA (likeANOVA) is very robust to heterogeneity for balanced designs but not unbalanced designs.
# Fortunately, in this case, it is homoegenous

# using Permanova anyway

pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'euclidean', p.adjust.m ='bonferroni'){
    library(vegan)

    co = combn(unique(as.character(factors)),2)
    pairs = c()
    F.Model =c()
    R2 = c()
    p.value = c()

    for(elem in 1:ncol(co)){
        if(sim.function == 'daisy'){
            library(cluster)
            x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
        } else {
            x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)
        }

        ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
        pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
        F.Model =c(F.Model,ad$aov.tab[1,4]);
        R2 = c(R2,ad$aov.tab[1,5]);
        p.value = c(p.value,ad$aov.tab[1,6])
    }

    p.adjusted = p.adjust(p.value,method=p.adjust.m)
    sig = c(rep('',length(p.adjusted)))
    sig[p.adjusted <= 0.05] <-'.'
    sig[p.adjusted <= 0.01] <-'*'
    sig[p.adjusted <= 0.001] <-'**'
    sig[p.adjusted <= 0.0001] <-'***'

    pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
    print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
    return(pairw.res)

}

GHGes <- river[,46:48]
GHGes_dis_matrix <- daisy(GHGes, metric = "euclidean", stand = TRUE)

set.seed(2805)

permanova_river_phys <- adonis(GHGes_dis_matrix ~ T_w + DO + pH + EC+ Sal + Turb + Chlr,
                               data = river, permutations = 999,
                               method = "euclidean", strata = river$River)
summary(permanova_river_phys)
permanova_river_phys
permanova_river_phys$aov.tab[,6]

pairwise.adonis(river[,6:13], river$River) #  physical


# more about the


