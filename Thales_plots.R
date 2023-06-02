
# Packages
library(ggplot2)
library(egg)
library(tidyr)

# Plot Coverage
coverage<-read.csv("Thales_coverage.csv")
experiment<-unlist(lapply(strsplit(colnames(coverage)[-c(21:22)], "t"), "[[", 2))
experiment<-as.numeric(experiment)

# For main text lets just plot original, hetero and hetero with IT @ n = 25
cover_plot<-data.frame(coverage=unlist(c(coverage[8,-c(21,22)], coverage[11,-c(21,22)], coverage[17,-c(21,22)])), method=sort(rep(c("Het. LM", "Het. LM (incl. IT)", "Hom. LM"), 20)), ratio=experiment)
c(8,11,17)

a<-ggplot(cover_plot, aes(x=log(ratio), y=coverage, color=method)) + theme_bw() +
	ylim(0, 1) + geom_hline(yintercept=0.95) + geom_vline(xintercept=0, color="red") + geom_point(size=3) +
	xlab("Log Ratio Nutrients") + ylab("Coverage") + labs(subtitle="A.") +
	theme(legend.position=c(0.2,0.135), legend.title=element_blank(), text=element_text(size=15), legend.background=element_blank()) +
	scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

# SE
se<-read.csv("Thales_SES.csv")

# For main text lets just plot original, hetero and hetero with IT @ n = 25
se_plot<-data.frame(se=unlist(c(se[2,-c(21,22)], se[17,-c(21,22)])), method=sort(rep(c("Delta Method", "Simulated"), 20)), ratio=experiment)

b<-ggplot(se_plot, aes(x=log(ratio), y=se, color=method)) + theme_bw() +
	ylim(1, 14) + geom_vline(xintercept=0, color="red") + geom_point(size=3) +
	xlab("Log Ratio Nutrients") + ylab("Standard Error") + labs(subtitle="B.") +
	theme(legend.position=c(0.175,0.93), legend.title=element_blank(), text=element_text(size=15), legend.background=element_blank()) +
	scale_color_manual(values=c("#999999", "black"))


######################################################################
######################### Drosophila Data ############################
######################################################################

# Load functions from Header file
source("Thales_Header.R")

# Code for the first part of the analysis is taken from the supplement for: 

########## Nutrigonometry IV: Thales' theorem #################
## Scientific Reports
## Juliano Morimoto
#* Please cite the original paper if using the approach or the functions

## Packages - some redundancy removed by AMS
library(dplyr)
library(tidyr)
library(sp)
library(patchwork)
library(stringr)
library(purrr)

# Added by AMS
library(sandwich)
library(plyr)

### Exercise 1: Drosophila melanogaster data from Lee et al.

leedt <- read.csv("41598_2023_34722_MOESM2_ESM.csv", 
                  header = TRUE, 
                  strip.white = TRUE, 
                  stringsAsFactors = FALSE,
                  na.string = 'NA') %>%
  na.omit()

## Fixing the class of the varialbe in the data ##
leedt$carb_eaten <- as.numeric(leedt$carb_eaten)
leedt$protein_eaten <- as.numeric(leedt$protein_eaten)
leedt$dailyeggs<- as.numeric(leedt$dailyeggs)
leedt$lifespan <- as.numeric(leedt$lifespan)
leedt$lifetimeegg <- as.numeric(leedt$lifetimeegg)

##
leedt <- na.omit(leedt)

### Lee intake target ###
leedt_intakeTarget <- data.frame(y = (148/5)*4, 
                                 x = (148/5)*1)
                                 
# Added by AMS - create a psuedo IT to test how relocating the IT induces variance in the angle
leedt_intakeTarget_FALSE <- data.frame(y = (148/5)*1, 
                                 x = (148/5)*1)
                                 
# AMS added the psudo angle and total food intake and vector from the origin        
compromise_angleDros <- data.frame(ratio = leedt$Ratio,
                                   beta = calculate_Thales(leedt[c("carb_eaten", "protein_eaten")], leedt_intakeTarget),
                                   false_angle = calculate_Thales(leedt[c("carb_eaten", "protein_eaten")], leedt_intakeTarget_FALSE),
                                   total_intake = leedt$carb_eaten + leedt$protein_eaten,
                                   vector = sqrt(leedt$carb_eaten^2 + leedt$protein_eaten^2))


### Testing statistically - original model
compromise_angleDros$ratio <- as.factor(compromise_angleDros$ratio)
model_angle_Dros<- lm(beta - 90 ~ 0 + ratio, data = compromise_angleDros)
summary(model_angle_Dros)
anova(model_angle_Dros)
confint(model_angle_Dros)

# This matches the in text from the original paper - so far so good

# Now lets calculate the variance in the angle, the food intake and the 'pseudo angle'
my_summary<-ddply(compromise_angleDros, .(ratio), summarise, sd_intake=sd(total_intake), sd_angle=sd(beta), sd_false=sd(false_angle))

# Order for plotting
plot_order<-c("(0:1)", "(1:16)", "(1:8)", "(1:4)", "(1:2)", "(1:1)", "(2:1)")
my_summary<-my_summary[match(plot_order, my_summary$ratio),]

# Make data in long format for plotting
plot_sds<-data.frame(ratio = rep(my_summary$ratio, 3), sds = c(my_summary[,2], my_summary[,3], my_summary[,4]), statistic = c(rep("Total Intake", 7), rep("Angle, IT = (1:4)", 7), rep("Angle, False IT = (1:1)", 7)))
plot_sds$ratio<-factor(plot_sds$ratio, levels=unique(plot_sds$ratio))

c<-ggplot(plot_sds, aes(x=ratio, y=sds, color=statistic, group=statistic)) +
	geom_point(size=3) +
	geom_path() +
	theme_bw() + 
	xlab("PC Ratio") + ylab("Standard Deviation") + labs(subtitle="C.") +
	theme(legend.position=c(0.24,0.92), legend.title=element_blank(), text=element_text(size=15), legend.background=element_blank()) +
	scale_color_manual(values=c("#E69F00", "#009E73", "#000000"))


pdf("Figure.1.pdf", height=5, width=15)
ggarrange(a,b,c, ncol=3, nrow=1)
dev.off()

# Do welch test for food group on the intake and vector from the origin
library(onewaytests)
closest_distance_test<-welch.test(vector ~ ratio, data=compromise_angleDros)
equal_distance_test<-welch.test(total_intake ~ ratio, data=compromise_angleDros)

# Supplementary figure for other methods and n
# Convert to long format

long_cover<-gather(coverage, ratio, cov, diet0.13:diet8)
long_cover$ratio<-as.numeric(unlist(lapply(strsplit(long_cover$ratio, "t"), "[[", 2)))

a<-ggplot(long_cover[which(long_cover$n == 5),], aes(x=log(ratio), y=cov, color=method)) + theme_bw() +
	ylim(0, 1) + geom_hline(yintercept=0.95) + geom_vline(xintercept=0, color="red") + geom_point(size=3) +
	xlab("Log Ratio Nutrients") + ylab("Coverage") + labs(subtitle="A. (n = 5)") +
	theme(legend.position=c(0.23,0.2), legend.title=element_blank(), text=element_text(size=15), legend.background=element_blank())

b<-ggplot(long_cover[which(long_cover$n == 25),], aes(x=log(ratio), y=cov, color=method)) + theme_bw() +
	ylim(0, 1) + geom_hline(yintercept=0.95) + geom_vline(xintercept=0, color="red") + geom_point(size=3) +
	xlab("Log Ratio Nutrients") + ylab("Coverage") + labs(subtitle="B. (n = 25)") +
	theme(legend.position="none", legend.title=element_blank(), text=element_text(size=15), legend.background=element_blank())

c<-ggplot(long_cover[which(long_cover$n == 50),], aes(x=log(ratio), y=cov, color=method)) + theme_bw() +
	ylim(0, 1) + geom_hline(yintercept=0.95) + geom_vline(xintercept=0, color="red") + geom_point(size=3) +
	xlab("Log Ratio Nutrients") + ylab("Coverage") + labs(subtitle="C. (n = 50)") +
	theme(legend.position="none", legend.title=element_blank(), text=element_text(size=15), legend.background=element_blank())

pdf("SuppFig.1.pdf", height=5, width=15)
ggarrange(a,b,c, ncol=3, nrow=1)
dev.off()




