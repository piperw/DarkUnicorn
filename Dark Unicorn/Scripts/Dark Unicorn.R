#########################################
#  Dark Unicorn Analysis                # 
#  20 September 2018                    #
#  Created by Piper Wallingford         #  
#  Most recent updates: Oct 2021        #  
#########################################


################################### Surveys ####################################

# clear workspace
rm(list=ls())

# load libraries
library("car")
library("lme4")
library("lmerTest")
library("ggeffects")
library("bbmle")
library("tidyverse")
library("blme")
library("ggiraphExtra")
library("ggiraph")
library("ggpubr")
library("maps")
library("mapdata")
library("raster")
library("plotrix")
library("marmap")
library("gridBase")
library("ggplot2")


# read in data
Data <- read.csv("Data/SurveyData.csv") 


################################ Data Summaries  ###########################

# calculate densities
Data$Mdens = Data$M.lugubris/Data$Area
Data$Adens = Data$A.spirata/Data$Area
Data$Ndens = Data$N.emarginata/Data$Area
Data$Other.dens = (Data$Ceres + Data$Conus + Data$Max + Data$Olive)/Data$Area
Data$Native.dens = (Data$A.spirata + Data$N.emarginata)/Data$Area


# Split by Mex  presence
Data$M.PA <- NULL
Data$M.PA = ifelse(Data$M.lugubris!=0, 1, 0)
Data$N.PA <- NULL
Data$N.PA = ifelse(Data$Native!=0, 1, 0)


# Set Factors
Data <- Data %>% 
  mutate(Season=as.factor(Season), Site=as.factor(Site), Transect=as.factor(Transect), TH=as.numeric(TH))


# Average across transects, seasons, TH
Site.Data <- Data %>% 
  group_by(Site) %>%
  summarize(mean.as = mean(Adens),
            mean.ne = mean(Ndens),
            mean.other = mean(Other.dens),
            mean.native = mean(Adens + Ndens + Other.dens),
            mean.mex = mean(Mdens),
            se.mex = sd(Mdens)/sqrt(n()),
            se.nat = sd(Adens + Ndens + Other.dens)/sqrt(n()))

# Average across transects
Trnsct.Data <- Data %>% 
  group_by(Season, Site, TH) %>%
  summarize(mean.native= mean(Native.dens),
            se.nat = sd(Native.dens)/sqrt(n()),
            mean.mex = mean(Mdens),
            se.mex = sd(Mdens)/sqrt(n())) %>%
  mutate(N.PA=ifelse(mean.native>0,1,0)) %>%
  mutate(M.PA=ifelse(mean.mex>0,1,0))


# Site Data
Site.Means <- Site.Data[,-c(5,7,8)]
Site.Means <- gather(Site.Means, key, value, -Site)
colnames(Site.Means) <- c("Site","Species","Mean")
Site.Means$Species <- replace(Site.Means$Species, Site.Means$Species =="mean.mex", "Mexacanthina")
Site.Means$Species <- replace(Site.Means$Species, Site.Means$Species == "mean.as", "Acanthinucella")
Site.Means$Species <- replace(Site.Means$Species, Site.Means$Species == "mean.ne", "Nucella")
Site.Means$Species <- replace(Site.Means$Species, Site.Means$Species == "mean.other", "Other")
Site.Means$Site <- factor(Site.Means$Site, levels = c("LC","CC","SH","HP","TI","TS","DP","CR","SC","SR"), ordered = TRUE)
Site.Means$Status <- ifelse(Site.Means$Species == "Mexacanthina", "Mexacanthina", "Native")
Site.Means$Species <- factor(Site.Means$Species, levels = c("Acanthinucella","Nucella","Other","Mexacanthina"), ordered = TRUE)

Site.SE <- Site.Data[,-c(2:6)]
Site.SE <- gather(Site.SE, key, value, -Site)
colnames(Site.SE) <- c("Site","Species","SE")
Site.SE$Species <- replace(Site.SE$Species, Site.SE$Species =="se.nat", "Native")
Site.SE$Species <- replace(Site.SE$Species, Site.SE$Species =="se.mex", "Mexacanthina")
Site.SE$Site <- factor(Site.SE$Site, levels = c("LC","CC","SH","HP","TI","TS","DP","CR","SC","SR"), ordered = TRUE)

Site.Nat <- Site.Data[,-c(2:4,7,8)]
Site.Nat <- gather(Site.Nat, key, value, -Site)
colnames(Site.Nat) <- c("Site","Species","Mean")
Site.Nat$Species <- replace(Site.Nat$Species, Site.Nat$Species =="mean.native", "Native")
Site.Nat$Species <- replace(Site.Nat$Species, Site.Nat$Species =="mean.mex", "Mexacanthina")
Site.Nat$Site <- factor(Site.Nat$Site, levels = c("LC","CC","SH","HP","TI","TS","DP","CR","SC","SR"), ordered = TRUE)

Site.errors <- merge(Site.Nat, Site.SE, by = c("Site","Species"))


# Tide Height Data
TH.Site <- Data
TH.Site$MPA <- NULL
TH.Site$MPA <- ifelse(TH.Site$Site == "LC" | TH.Site$Site == "CC" | TH.Site$Site == "SH" | TH.Site$Site == "DP", "Mexacanthina Absent", "Mexacanthina Present")
TH.Site <- TH.Site %>%
  group_by(MPA, TH) %>%
  summarise(native = mean(Native.dens),
            mex = mean(Mdens),
            se.nat = sd(Native.dens)/sqrt(n()),
            se.mex = sd(Mdens)/sqrt(n()))
TH.Site$MPA <-factor(TH.Site$MPA, levels = c("Mexacanthina Absent", "Mexacanthina Present"))

TH.nat <- TH.Site[,c(1:3,5)]
TH.nat$Species <- "Native Whelks"
colnames(TH.nat)[3] <- "Density"
colnames(TH.nat)[4] <- "se"

TH.mex <- TH.Site[,c(1:2,4,6)]
TH.mex$Species <- "Mexacanthina"
colnames(TH.mex)[3] <- "Density"
colnames(TH.mex)[4] <- "se"

TH.MPA <- rbind(TH.nat, TH.mex)


# mapping datasets
data(us.cities)
sites <- read.csv("Data/Sites.csv")


# Clean up environment
rm(Site.Data, Site.SE, Site.Nat, TH.nat, TH.mex, TH.Site)


######################## Explorative Plots ###################################

variables <- Data[,c(4,25,29:31)]
pairs(variables)

ggplot(data = Data, aes(x=as.factor(TH), y=Native.dens, fill=factor(M.PA))) +
  geom_boxplot() 
ggplot(data = Data, aes(x=as.factor(TH), y=Mdens)) +
  geom_boxplot() 

ggplot(Data, aes(x=Native,color=as.factor(M.PA)))+
  geom_histogram() 

rm(variables)

################################ Site Plots #####################################

# Create plot
#png("Plots/Figure1a.png", width = 2000, height = 2400, res = 300)
#par(mar = c(5,5,4,2))

# prepare the plot
plot(sites$Long, sites$Lat, cex = 0.25, xlab = "Longitude", ylab = "Latitude", xlim= c(-119, -117), ylim = c(32.5,34.25), type = 'n', cex.axis = 1.25, cex.lab = 1.5)

#create main map
rect(-119.5, 32.25, -116.5, 34.5, density = NULL, angle = 45,
     col = "slategray1", border = NULL)

map.plot <- map("worldHires", xlim=c(-119, -117), ylim = c(32.5,34.25), col = 'darkseagreen3', fill = TRUE, border = FALSE, add = TRUE)

# add counties + mexico
map("county", region = c('california,san diego', 'california,orange', 'california,los angeles', 'california,ventura', 'california,riverside','california,san bernardino', 'california,kern', "california,santa barbara"), fill = TRUE, col = 'darkseagreen3', add = TRUE, border = TRUE)

map("worldHires","mexico", col="gray60", fill=TRUE, add=TRUE)  #add the adjacent parts of mexico 

# add city markers
map.cities(us.cities, country="CA", minpop = 1000000, pch = 17, label = FALSE, cex = 3)


# Add sites
points(x = jitter(sites$Long[sites$Species == "N"]), y = sites$Lat[sites$Species == "N"], pch=21, col = "black", bg="aquamarine4", cex= 1.5, lwd = 0.75)
points(x = jitter(sites$Long[sites$Species == "MN"]), y = sites$Lat[sites$Species == "MN"], pch=21, col = "black", bg="mediumpurple3", cex= 1.5, lwd = 0.75)
points(sites$Long[sites$Species == "None"], sites$Lat[sites$Species == "None"], pch=21, col="black", bg = "white", cex=1.5, lwd = 0.75)

# add labels
text(x = sites$Long[sites$Labels == "DR"] +0.07, y = sites$Lat[sites$Labels == "DR"]-0.01, labels = sites$Site[sites$Labels == "DR"], cex=1)
text(x = sites$Long[sites$Labels == "DL"] -0.05, y = sites$Lat[sites$Labels == "DL"]-0.03, labels = sites$Site[sites$Labels == "DL"], cex=1)
text(x = sites$Long[sites$Labels == "UR"] +0.05, y = sites$Lat[sites$Labels == "UR"]+0.03, labels = sites$Site[sites$Labels == "UR"], cex=1)
text(x = sites$Long[sites$Labels == "R"] +0.08, y = sites$Lat[sites$Labels == "R"], labels = sites$Site[sites$Labels == "R"], cex = 1)
text(x = sites$Long[sites$Labels == "L"]-0.08, y = sites$Lat[sites$Labels == "L"], labels = sites$Site[sites$Labels == "L"], cex=1)
text(x = sites$Long[sites$Labels == "DDL"]-0.01, y = sites$Lat[sites$Labels == "DDL"]-0.045, labels = sites$Site[sites$Labels == "DDL"], cex=1)
text(x = sites$Long[sites$Labels == "UUR"] +0.025, y = sites$Lat[sites$Labels == "UUR"]+0.04, labels = sites$Site[sites$Labels == "UUR"], cex=1)
text(x = sites$Long[sites$Labels == "UL"] -0.06, y = sites$Lat[sites$Labels == "UL"]-0.02, labels = sites$Site[sites$Labels == "UL"], cex=1)

#add legend
legend("bottomleft", expression(paste(italic("Mexacanthina"), paste(" and native whelks")), "Native whelks only", "No whelks"), pch = 21, pt.bg = c("mediumpurple3","aquamarine4","white"), bty = "n", cex = 1.25)

#dev.off()


#site densities
#png("Plots/Figure1b.png", width = 2000, height = 2400, res = 300)
ggplot(data=Site.Means, aes(y = Mean, x = Status, fill = Species)) + 
  geom_bar(stat="identity") + 
  facet_grid(Site ~ ., switch = "y") +
  theme(strip.background = element_blank(),
        strip.text.y = element_text(size = 14),
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0,"mm"),
        axis.text.y=element_blank(),
        axis.text.x=element_text(size=16),
        axis.title=element_text(size=18),
        axis.ticks.y=element_blank()) + 
  coord_flip() +
  scale_y_continuous(expand = c(0, 0), limits = c(0,5), breaks = scales::pretty_breaks(n = 5)) +
  labs(y = expression(paste('Whelks (per m'^{2},')')), x = "Site") +
  theme(axis.line.x = element_line(color="gray", size = 0.25),
        legend.position = c(1.15,1), 
        legend.justification = c("right", "top"), 
        legend.box.just = "right",
        legend.title.align=0.35,
        legend.title=element_text(size=16),
        legend.text=element_text(size=14, face = "italic")) +
  guides(fill=guide_legend(ncol=3)) + 
  geom_errorbar(data=Site.errors ,aes(x=Species, ymin=Mean-SE, ymax=Mean+SE),
                width=.05, position = position_dodge(0.9)) +
  scale_fill_manual(breaks=c("Mexacanthina","Acanthinucella","Nucella","Other",
                             "Native"),
                    values = c("mediumpurple3","aquamarine4","peachpuff","lemonchiffon","white"))
#dev.off()


################################ TH ###########################################
# Mex and All Nat by M Presence

#png("Plots/THdistributions.PA.png", width = 2400, height = 1600, res = 300)
ggplot(data=TH.MPA, aes(x = TH, y = Density, fill = Species, color = Species)) +
  facet_grid(~MPA) +
  geom_point(aes(y = Density, color = Species)) +
  geom_errorbar(aes(ymin=Density-se, ymax=Density+se, color = Species), width = 0.05) +
  geom_ribbon(aes(ymin = 0, ymax = Density, color = Species), alpha = 0.6) +
  scale_color_manual(values = c("mediumpurple3","aquamarine4")) + 
  scale_fill_manual(values = c("mediumpurple3","aquamarine4")) +
  xlab("Tidal Elevation (m)") +
  ylab(expression(paste('Density (per m'^{2},')'))) +
  scale_y_continuous(expand = c(0.0,0.05)) +
  scale_x_continuous(expand = c(0.001,0.001)) +
  theme(strip.background = element_blank(),
        strip.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 12),
        panel.border = element_rect(colour = "gray", fill=NA, size=1),
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(1,"cm"),
        axis.text.y=element_text(size = 14),
        axis.text.x=element_text(size=14),
        axis.title=element_text(size=16)) + 
  theme(legend.position = c(0.05, 0.95), 
        legend.key = element_rect(colour = 'white', fill = 'white', size = 0.5),
        legend.justification = c("left", "top"), 
        legend.box.just = "left",
        legend.title.align=0.35,
        legend.title=element_text(size=12),
        legend.text=element_text(size=12, face = "italic"))
#dev.off()


############################# Hurdle Models ###################################

# 1. Density models
binomial1 <- glm(N.PA ~ mean.mex * TH, data = Trnsct.Data, family = binomial(link = logit))
gamma1 <- glm(mean.native ~ mean.mex * TH, data = Trnsct.Data[Trnsct.Data$N.PA == 1,], family = Gamma(link = log), control = list(maxit = 50))

#binomial coefs
bin_coef1 <- plogis(is.numeric(coef(binomial1)[[1]]))
bin_confit1 <- plogis(confint(binomial1, method = "boot"))
Anova(binomial1, type = 3, singular.ok = TRUE)

#gamma coefs
gamma_coef1 <- exp(coef(gamma1)[[1]])
gamma_confit1 <- exp(confint(gamma1, method = "boot"))  # method = "boot"
Anova(gamma1, type = 3, singular.ok = TRUE)

plot(mean.mex ~ TH, data = Trnsct.Data)
plot(mean.native ~ TH, data = Trnsct.Data)
plot(mean.native ~ mean.mex, data = Trnsct.Data)


# 2. Presence/Absence models
binomial2 <- glm(N.PA ~ M.PA * TH, data = Trnsct.Data, family = binomial(link = logit))
gamma2 <- glm(mean.native ~ M.PA * TH, data = Trnsct.Data[Trnsct.Data$N.PA == 1,], family = Gamma(link = log))

# binomial coefs
bin_coef2 <- plogis(is.numeric(coef(binomial2)[[1]]))
bin_confit2 <- plogis(confint(binomial2, method = "boot"))
Anova(binomial2, type = 3, singular.ok = TRUE)

# gamma coefs
gamma_coef2 <- exp(coef(gamma2)[[1]])
gamma_confit2 <- exp(confint(gamma2, method = "boot"))
Anova(gamma2, type = 3, singular.ok = TRUE)


############################# Hurdle Model Plots ####################################
## density
plot(ggpredict(binomial1, terms = c("TH", "mean.mex [0,0.03278689,14.7765]"), interactive = TRUE, se = TRUE, facet = TRUE))
plot(ggpredict(gamma1, c("TH","mean.mex [0,14.7765]"), interactive = TRUE, se = TRUE, facet = TRUE))

## presence/absence
plot(ggpredict(binomial2, c("TH","M.PA"), interactive = TRUE, se = TRUE, facet = TRUE))
plot(ggpredict(gamma2, c("TH","M.PA"), interactive = TRUE, se = TRUE, facet = TRUE))

# Summarize distributions
summary(Trnsct.Data$mean.mex)
summary(Trnsct.Data$mean.mex[Trnsct.Data$mean.mex != 0])
summary(Trnsct.Data$mean.native)
summary(Trnsct.Data$mean.native[Trnsct.Data$mean.native != 0])


# Presence/absence (PA)
#png("Plots/Figure3c.png", width = 2400, height = 1600, res = 300)
a <- plot(ggpredict(binomial2, c("TH","M.PA"), interactive = TRUE, se = TRUE, facet = TRUE)) + 
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE) +
  scale_color_manual(name = expression(paste(italic("Mexacanthina "), "Presence")), labels = c("Absent", "Present"), values = c("aquamarine4","mediumpurple4")) +
  scale_fill_manual(values = c("aquamarine4","mediumpurple4")) +
  labs(x = "Tidal Elevation (m)", y = "Probability of Native Whelk Presence") +
  scale_x_continuous(expand = c(0.005, 0.005)) + 
  theme(strip.background = element_blank(),
        strip.text.y = element_text(size = 12),
        panel.border = element_rect(colour = "gray", fill=NA, size=1),
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0,"line"),
        axis.text.y=element_text(size = 14),
        axis.text.x=element_text(size=14),
        axis.title=element_text(size=15)) +
  theme(plot.margin = unit(c(0, 0.25, 0.5, 1), "cm")) +
  theme(axis.line.x = element_line(color="gray", size = 0.25),
        legend.position = "top", 
        legend.justification = "top", 
        legend.box.just = "right",
        legend.title=element_text(size=14),
        legend.key = element_rect(fill = NA),
        legend.text=element_text(size=12)) 
#dev.off()


# Density (PA)
#png("Plots/Figure3d.png", width = 2400, height = 1600, res = 300)
c <- plot(ggpredict(gamma2, terms = c("TH", "M.PA"), interactive = TRUE, se = TRUE, facet = TRUE)) + 
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE) +
  scale_color_manual(name = expression(paste(italic("Mexacanthina "), "Presence")), labels = c("Absent", "Present"), values = c("aquamarine4","mediumpurple4")) +
  scale_fill_manual(values = c("aquamarine4","mediumpurple4")) +
  labs(x = "Tidal Elevation (m)", y = "Density of Native Whelks") +
  scale_x_continuous(expand = c(0.005, 0.005)) + 
  scale_y_continuous(limits = c(0,3), breaks = c(0.5,1.0,1.5,2.0,2.5,3.0)) +
  theme(strip.background = element_blank(),
        strip.text.y = element_text(size = 12),
        panel.border = element_rect(colour = "gray", fill=NA, size=1),
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0,"line"),
        axis.text.y=element_text(size = 14),
        axis.text.x=element_text(size=14),
        axis.title=element_text(size=15),
        axis.title.y=element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(plot.margin = unit(c(0, 0.25, 0.5, 1), "cm")) +
  theme(axis.line.x = element_line(color="gray", size = 0.25),
        legend.position = "top", 
        legend.justification = "top", 
        legend.box.just = "right",
        legend.title=element_text(size=14),
        legend.key = element_rect(fill = NA),
        legend.text=element_text(size=12)) 
#dev.off()


# Presence/absence (density)
#png("Plots/Figure3a.png", width = 2400, height = 1600, res = 300)
b <- plot(ggpredict(binomial1, terms = c("TH", "mean.mex [0,2.55]"), interactive = TRUE, se = TRUE, facet = TRUE)) + 
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE) +
  scale_color_manual(name = expression(paste(italic("Mexacanthina "), "Density")), labels = c("None", "Mean (2.55)"), values = c("aquamarine4","mediumpurple4")) +
  scale_fill_manual(name = expression(paste(italic("Mexacanthina "), "Density")), labels = c("None", "Mean (2.55)"), values = c("aquamarine4","mediumpurple4")) +
  labs(x = "Tidal Elevation (m)", y = "Probability of Native Whelk Presence") +
  scale_x_continuous(expand = c(0.005, 0.005)) + 
  theme(strip.background = element_blank(),
        strip.text.y = element_text(size = 12),
        panel.border = element_rect(colour = "gray", fill=NA, size=1),
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0,"line"),
        axis.text.y=element_text(size = 14),
        axis.text.x=element_text(size=14),
        axis.title=element_text(size=15)) +
  theme(plot.margin = unit(c(0.1, 0.25, 0.5, 1), "cm")) +
  theme(axis.line.x = element_line(color="gray", size = 0.25),
        legend.position = "top", 
        legend.justification = "top", 
        legend.box.just = "right",
        legend.title=element_text(size=14),
        legend.key = element_rect(fill = NA),
        legend.text=element_text(size=12)) 
#dev.off()


# Density (density)
#png("Plots/Figure3b.png", width = 2400, height = 1600, res = 300)
d <- plot(ggpredict(gamma1, terms = c("TH", "mean.mex[0,2.54952]"), interactive = TRUE, se = TRUE, facet = TRUE)) + 
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE) +
  scale_color_manual(name = expression(paste(italic("Mexacanthina "), "Density")), labels = c("None", "Mean (2.55)"), values = c("aquamarine4","mediumpurple4")) +
  scale_fill_manual(name = expression(paste(italic("Mexacanthina "), "Density")), labels = c("None", "Mean (2.55)"), values = c("aquamarine4","mediumpurple4")) +
  labs(x = "Tidal Elevation (m)", y = "Density of Native Whelks") +
  scale_x_continuous(expand = c(0.005, 0.005)) + 
  scale_y_continuous(limits = c(0,3), breaks = c(0.5,1.0,1.5,2.0,2.5,3.0)) +
  theme(strip.background = element_blank(),
        strip.text.y = element_text(size = 12),
        panel.border = element_rect(colour = "gray", fill=NA, size=1),
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0,"line"),
        axis.text.y=element_text(size = 14),
        axis.text.x=element_text(size=14),
        axis.title=element_text(size=15),
        axis.title.y=element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(plot.margin = unit(c(0, 0.25, 0.5, 1), "cm")) +
  theme(axis.line.x = element_line(color="gray", size = 0.25),
        legend.position = "top", 
        legend.justification = "top", 
        legend.box.just = "right",
        legend.title=element_text(size=14),
        legend.key = element_rect(fill = NA),
        legend.text=element_text(size=12)) 
#dev.off()

ggarrange(a,c,b,d)

#tiff("Plots/Figure3.tif", width = 300, height = 250, units = "mm", res = 300)
#ggarrange(a, b, c, d, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
#dev.off()


#### Supplement 
# Presence/absence (density with all densities)
#png("Plots/Supp.PA.png", width = 2400, height = 1600, res = 300)
plot(ggpredict(binomial1, terms = c("TH", "mean.mex [0, 0.03278689, 2.549517, 14.77654]"), interactive = TRUE, se = TRUE, facet = TRUE)) + 
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE) +
  scale_color_manual(name = expression(paste(italic("Mexacanthina "), "Density")), labels = c("None", "Minimum (0.03)", "Mean (2.55)", "Maximum (14.78)"), values = c("aquamarine4", "mediumorchid1", "mediumpurple4", "#290924")) +
  scale_fill_manual(name = expression(paste(italic("Mexacanthina "), "Density")), labels = c("None", "Minimum (0.03)", "Mean (2.55)", "Maximum (14.78)"), values = c("aquamarine4", "mediumorchid1", "mediumpurple4", "#290924")) +
  labs(x = "Tidal Elevation (m)", y = "Probability of Native Whelk Presence") +
  scale_x_continuous(expand = c(0.005, 0.005)) + 
  theme(strip.background = element_blank(),
        strip.text.y = element_text(size = 12),
        panel.border = element_rect(colour = "gray", fill=NA, size=1),
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0,"line"),
        axis.text.y=element_text(size = 14),
        axis.text.x=element_text(size=14),
        axis.title=element_text(size=16)) +
  theme(axis.line.x = element_line(color="gray", size = 0.25),
        legend.position = c(0.95,0.95), 
        legend.justification = c("right", "top"), 
        legend.box.just = "right",
        legend.title.align=0.35,
        legend.title=element_text(size=14),
        legend.key = element_rect(fill = NA),
        legend.text=element_text(size=12)) 
#dev.off()

# Density (density with all densities)
#png("Plots/Supp.Density.png", width = 2400, height = 1600, res = 300)
plot(ggpredict(gamma1, terms = c("TH", "mean.mex [-0.1, 0.03278689, 2.549517, 4.77654]"), interactive = TRUE, se = TRUE, facet = TRUE)) + 
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE) +
  scale_color_manual(name = expression(paste(italic("Mexacanthina "), "Density")), labels = c("None", "Minimum (0.03)", "Mean (2.55)", "Maximum (14.78)"), values = c("aquamarine4", "mediumorchid1", "mediumpurple4", "#290924")) +
  scale_fill_manual(name = expression(paste(italic("Mexacanthina "), "Density")), labels = c("None", "Minimum (0.03)", "Mean (2.55)", "Maximum (14.78)"), values = c("aquamarine4", "mediumorchid1", "mediumpurple4", "#290924")) +
  labs(x = "Tidal Elevation (m)", y = "Native Whelk Density") +
  scale_x_continuous(expand=c(0,0), limits=c(0,1.75)) +
  scale_y_continuous(expand=c(0,0), limits=c(-10,10)) +
  coord_cartesian(xlim=c(0,1.75), ylim=c(0,2)) +
  scale_x_continuous(expand = c(0.005, 0.005)) + 
  theme(strip.background = element_blank(),
        strip.text.y = element_text(size = 12),
        panel.border = element_rect(colour = "gray", fill=NA, size=1),
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0,"line"),
        axis.text.y=element_text(size = 14),
        axis.text.x=element_text(size=14),
        axis.title=element_text(size=16)) +
  theme(axis.line.x = element_line(color="gray", size = 0.25),
        legend.position = c(0.95,0.95), 
        legend.justification = c("right", "top"), 
        legend.box.just = "right",
        legend.title.align=0.35,
        legend.title=element_text(size=14),
        legend.key = element_rect(fill = NA),
        legend.text=element_text(size=12)) 
#dev.off()





################################ Hunger Games ###############################

# clear workspace
rm(list=ls())

# load libraries
library("multcomp")
library("agricolae")

## read in Data
Whelks <- read.csv("Data/WhelkData_Final.csv")
Mussels <- read.csv("Data/Muss_Biomass.csv")
Barnacles <- read.csv("Data/Barn_Biomass.csv")


## Convert prey to biomass

Barn.mod <- lm(Mass ~ Balive, data = Barnacles)
qqnorm(resid(Barn.mod))
qqline(resid(Barn.mod))
summary(Barn.mod)

plot(Barnacles$Balive, Barnacles$Mass)
# barn.eq <- -0.224389 + 0.050620(alive) + 0.1871

#png("Plots/BarnMass.png", width = 2400, height = 1600, res = 300)
ggplot(data=Barnacles, aes(y = Mass, x = Balive)) + 
  geom_point(shape = 1) +
  geom_smooth(method='lm') +
  labs(x = "Number of Individuals", y = "Prey Biomass (g)") +
  theme(strip.background = element_blank(),
        strip.text.y = element_text(size = 12),
        panel.border = element_rect(colour = "gray", fill=NA, size=1),
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0,"line"),
        axis.text.y=element_text(size = 14),
        axis.text.x=element_text(size=14),
        axis.title=element_text(size=16))
#dev.off()


Muss.mod <- lm(Mass ~ Mlength + Mdepth + Mwidth, data = Mussels)
qqnorm(resid(Muss.mod))
qqline(resid(Muss.mod))
summary(Muss.mod)

Muss.mod2 <- lm(Mass ~ Mlength + Mwidth, data = Mussels)
qqnorm(resid(Muss.mod2))
qqline(resid(Muss.mod2))
summary(Muss.mod2)

plot(Mussels$Mlength, Mussels$Mass)
plot(Mussels$Mdepth, Mussels$Mass)
plot(Mussels$Mwidth, Mussels$Mass)


#png("Plots/MussMass.png", width = 2400, height = 1600, res = 300)
ggplot(data=Mussels, aes(y = Mass, x = Mlength)) + 
  geom_point(shape = 1) +
  geom_smooth(method='lm') +
  labs(x = "Length (mm)", y = "Prey Biomass (g)") +
  theme(strip.background = element_blank(),
        strip.text.y = element_text(size = 12),
        panel.border = element_rect(colour = "gray", fill=NA, size=1),
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0,"line"),
        axis.text.y=element_text(size = 14),
        axis.text.x=element_text(size=14),
        axis.title=element_text(size=16))
#dev.off()


# muss.eq <- -2.433694 + 0.123821(length) + 0.002142(depth) + 0.028417(width) + 0.1533
# muss.eq2 <- -2.43395 + 0.01849(length) + 0.02970(width) + 0.1515


## Experiment Prey

# Account for barnacle mortality

Barn.data <- read.csv("Data/BarnPrey.csv")
Barn.data$Ded <- Barn.data$Ialive - Barn.data$Falive

Control <- Barn.data[Barn.data$Species == "C",]
Barn.data <- Barn.data[Barn.data$Species != "C",]

# barn.eq <- -0.224389 + 0.050620(alive) + 0.1871

BarnMass <- Barn.data
BarnMass$Mass <- (Barn.data$Ded*0.050620) - 0.224389 + 0.1871
BarnMass <- BarnMass[,-c(8:10)]

# Mussels 

Muss.data <- read.csv("Data/MussPrey.csv")
Muss.data <- Muss.data[Muss.data$Species != "C",]
Muss.data <- Muss.data[Muss.data$Survived == "N",]

#muss.eq2 <- -2.43395 + 0.12418(length) + 0.02970(width) + 0.1515

MussMass <- Muss.data
MussMass$Mass <- (0.12418*MussMass$Mlength) + (0.02970*MussMass$Mwidth) - 2.43395 + 0.1515
MussMass <- MussMass[,-c(8:11)]


## Analyze Prey Consumed 

PreyMass <- rbind(MussMass,BarnMass)
PreyMass <- PreyMass %>%
  group_by(Tile, Pred, Prey, Density) %>%
  summarize(mass = sum(Mass))
PreyMass$mass <- ifelse(PreyMass$Density != "Low", PreyMass$mass/2, PreyMass$mass)


PredMass <- Whelks %>%
  group_by(Tile,Pred,Prey,Density) %>%
  summarize(PredMass = sum(iMass))
PredMass$PredMass <- ifelse(PredMass$Density != "Low", PredMass$PredMass/2, PredMass$PredMass)


BioMass <- merge(PreyMass,PredMass, by = c("Tile","Pred","Prey","Density"))
BioMass$Biomass <- BioMass$mass/BioMass$PredMass
BioMass <- BioMass %>% 
  mutate(Pred=as.factor(Pred), Prey=as.factor(Prey))



prey.ded <- aov(log(Biomass) ~ Pred * Prey, data = BioMass)
qqnorm(resid(prey.ded))
qqline(resid(prey.ded))
summary(prey.ded)
HSD.test(prey.ded, "Prey", console = TRUE)


Summary <- BioMass %>% 
  group_by(Pred, Prey, Density) %>%
summarize(mean.ded = mean(Biomass),
se.ded = sd(Biomass)/sqrt(n()))

log.Summary <- BioMass %>% 
  group_by(Pred, Prey, Density) %>%
  summarize(mean.ded = mean(log(Biomass)),
            se.ded = sd(log(Biomass))/sqrt(n()))

log.Summary$Prey <- factor(log.Summary$Prey, levels = c("Mu","MB","B"), ordered = TRUE)
log.Summary$Density <- factor(log.Summary$Density, levels = c("Low","High","Mixed"), ordered = TRUE)
log.Summary$Pred <- factor(log.Summary$Pred, levels = c("A","M","AA","MM","MA"), ordered = TRUE)


#### Plot

labels <- c(Low = "Individual", High = "With Conspecific", Mixed = "With Heterospecific")
labels.prey <- c(Mu = "Mussels", MB = "Mixed Prey", B = "Barnacles")

#png("Plots/Figure4a.png", width = 3200, height = 2200, res = 300)
ggplot(data=log.Summary, aes(y = mean.ded, x = Prey, fill = Pred)) + 
  geom_bar(position = "dodge", stat="identity") + 
  scale_x_discrete(labels = labels.prey) +
  facet_grid(. ~ Density, labeller=labeller(Density = labels)) +   # , switch = "y"
  scale_y_continuous(expand = c(0,0), limits = c(0,3), breaks = scales::pretty_breaks(n = 5)) +
  labs(x = "Prey Treatment", y = "Prey Biomass Consumed (g per Whelk Biomass in g)") +
  scale_fill_manual(name = "Whelk Species", labels = c("Acanthinucella",  "Mexacanthina"," "," ","Both Species"), values = c("aquamarine4","mediumpurple3","aquamarine4","mediumpurple3","steelblue4")) +
  geom_errorbar(data=log.Summary ,aes(x=Prey, ymin=mean.ded-se.ded, ymax=mean.ded + se.ded),width=.05, position = position_dodge(0.9)) + 
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 14),
        panel.border = element_rect(colour = "gray", fill=NA, size=1),
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0,"line"),
        axis.text.y=element_text(size = 14),
        axis.text.x=element_text(size=14),
        axis.title=element_text(size=16)) +
  theme(axis.line.x = element_line(color="gray", size = 0.25),
        legend.position = c(0.75,0.97), 
        legend.justification = c("left", "top"), 
        legend.box.just = "left",
        legend.title.align=0.35,
        legend.title=element_text(size=14),
        legend.text=element_text(size=12, face = "italic")) +
  guides(fill=guide_legend(ncol=1)) 
#dev.off()



## Analyze Whelks 

# Average across treatments
Whelks.Species <- Whelks %>% 
  group_by(Tile, Prey, Pred, Density, Species) %>%
  summarize(iMass = mean(iMass),
            fMass = mean(fMass),
            X.Mass = mean(X.Mass))


Whelks.Species$Pred <- ifelse(Whelks.Species$Density == "Mixed" & Whelks.Species$Species == "A", "AM", Whelks.Species$Pred)

Whelks.Species$Pred <- as.factor(Whelks.Species$Pred)
Whelks.Species$Prey <- as.factor(Whelks.Species$Prey)

mass.mod <- aov(X.Mass ~ Pred * Prey, data = Whelks.Species)
qqnorm(resid(mass.mod))
qqline(resid(mass.mod))
summary(mass.mod)
HSD.test(mass.mod, "Pred", console = TRUE)

PredSummary <- Whelks.Species %>% 
  group_by(Pred, Prey, Density) %>%
  summarize(mean.mass = mean(X.Mass),
            se.mass = sd(X.Mass)/sqrt(n()))

PredSummary$Pred <- factor(PredSummary$Pred, levels = c("A","M","AA","MM","AM","MA"), ordered = TRUE)
PredSummary$Density <- factor(PredSummary$Density, levels = c("Low","High","Mixed"), ordered = TRUE)


#png("Plots/Figure4b.png", width = 3200, height = 2200, res = 300)
ggplot(data=PredSummary, aes(y = mean.mass, x = Prey, fill = Pred)) + 
  geom_bar(position = "dodge", stat="identity") + 
  scale_x_discrete(labels = labels.prey) +
  facet_grid(. ~ Density, labeller=labeller(Density = labels)) +   # , switch = "y"
  scale_y_continuous(expand = c(0,0), limits = c(0,30), breaks = scales::pretty_breaks(n = 5)) +
  labs(x = "Prey Treatment", y = "Whelk Growth (% Change in Mass)") +
  scale_fill_manual(name = "Whelk Species", labels = c("Acanthinucella","Mexacanthina", "Acanthinucella", "Mexacanthina", "Acanthinucella", "Mexacanthina"), values = c("aquamarine4","mediumpurple3","aquamarine4","mediumpurple3","aquamarine4","mediumpurple3")) +
  geom_errorbar(data=PredSummary ,aes(x=Prey, ymin=mean.mass-se.mass, ymax=mean.mass+se.mass),width=.05, position = position_dodge(0.9)) + 
  #geom_text(aes(x = Density, y = mean.ded, label = Letters))
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 14),
        panel.border = element_rect(colour = "gray", fill=NA, size=1),
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0,"line"),
        axis.text.y=element_text(size = 14),
        axis.text.x=element_text(size=14),
        axis.title=element_text(size=16)) +
  theme(legend.position = "none")
#dev.off()


############################### Thermotolerance ##############################

# clear workspace
rm(list=ls())

# load libraries
library("survminer")
library("mgcv")
library("car")
library("MASS")
library("multcomp")


## Read in Data
TT.Data <- read.csv("Data/TT.csv")

TT.Data$Species <- factor(TT.Data$Species, levels = c("Nucella", "Acanthinucella", "Mexacanthina"))
TT.Data$Treatment <- as.numeric(TT.Data$Treatment)
TT.Data$Survival <- as.numeric(TT.Data$Survival)

TT.Data <- TT.Data[,c(2,4,5)]

TT.summary <- TT.Data %>%
  group_by(Species,Treatment) %>%
  summarize(prop.alive = mean(Survival),
            se.alive = sd(Survival)/sqrt(5))

mex.model <- glm(Survival ~ Treatment, data = TT.Data[TT.Data$Species == "Mexacanthina",], family = 'binomial')
summary(mex.model)
mex.dose <- dose.p(mex.model,p=0.5)
mex.mean <- matrix(mex.dose)
mex.se <- sd(TT.Data$Survival[TT.Data$Species == "Mexacanthina"])/sqrt(5)

ac.model <- glm(TT.Data$Survival[TT.Data$Species == "Acanthinucella"] ~ TT.Data$Treatment[TT.Data$Species == "Acanthinucella"], family = 'binomial')
ac.dose <- dose.p(ac.model,p=0.5)
ac.mean <- matrix(ac.dose)
ac.se <- sd(TT.Data$Survival[TT.Data$Species == "Acanthinucella"])/sqrt(5)

nuc.model <- glm(TT.Data$Survival[TT.Data$Species == "Nucella"] ~ TT.Data$Treatment[TT.Data$Species == "Nucella"], family = 'binomial')
nuc.dose <- dose.p(nuc.model,p=0.5)
nuc.mean <- matrix(nuc.dose)
nuc.se <- sd(TT.Data$Survival[TT.Data$Species == "Nucella"])/sqrt(5)

all.species <- glm(Survival ~ Species + Treatment, data = TT.Data, family = binomial)
summary(all.species)
Anova(all.species, type = 2)
summary(glht(all.species, mcp(Species="Tukey")))


# Results
Results.Matrix <- matrix(c( # first means
mex.mean, ac.mean, nuc.mean,
  # now SE
mex.se, ac.se, nuc.se), 3,2)

# label columns and rows
colnames(Results.Matrix)<-c("LT50","SE")
Results.Matrix



## Plots

#png("Plots/TT.all.png", width = 1600, height = 1200, res = 300)
ggplot(TT.Data, aes(x=Treatment, y=Survival, color = Species)) + 
  geom_jitter(width = 0.0, height = 0.075) +
  geom_smooth(method="glm", method.args=list(family="binomial"), se = FALSE, size = 1) +
  scale_color_manual(values = c("peachpuff","aquamarine4","mediumpurple3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "gray")) +
  theme(strip.background = element_blank(),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=14)) +
  labs(y = 'Proportion Survival', x = "Temperature") +
  scale_x_continuous(limits = c(18,42), breaks=seq(18,42,4)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray") +
  ylim(0,1) +
  theme(legend.position = c(0.05, 0.05), 
        legend.key = element_rect(colour = 'white', fill = 'white', size = 0.5),
                legend.justification = c("left", "bottom"), 
                legend.box.just = "left",
                legend.title.align=0.35,
                legend.title=element_text(size=12),
                legend.text=element_text(size=10, face = "italic"))
#dev.off()        


