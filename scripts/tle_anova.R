setwd("~/Mote_nutrient_experiment/data")

library(ggplot2)
library(dunn.test)
library(ggpubr)
library(dplyr)
library(agricolae)
library(multcompView)

tle <- read.csv("tle.csv")
tle$Weeks <- factor(tle$Weeks, levels = c("3", "6")) 
tle$trt_weeks <- factor(tle$trt_weeks, levels = c("NT3", "NT6", "LA3" , "LA6", "HA3", "HA6", "LN3", "LN6", "HN3", "HN6", "LP3", "LP6", "HP3", "HP6", "LC3", "LC6", "HC3", "HC6")) 

#tle is not normally distributed
ggdensity(tle$tle)
shapiro.test(tle$tle) #fails shapiro test, low p value
#passes fligner-killeen test of homogeneity of variances for non-normally distributed data, but fails levene test
model <- lm(tle ~ trt_weeks, data = tle)

plot(model, which = 2) #q-q plot
plot(model, which = 3) #variance

#log transform
tle <- mutate(tle, log_tle = log(tle))
model <- lm(log_tle ~ trt_weeks, data = tle)
ggdensity(tle$log_tle) #looks better
shapiro.test(tle$log_tle) #passes

#modeling relationship of tle to treatment
model <- lm(log_tle ~ trt_weeks, data = tle)
summary(model)
anova(model)
#omnibus test shows main effect is significant
 

tle3 <- tle[ which(tle$Weeks=='3'),]
tle6 <- tle[ which(tle$Weeks=='6'),]

#treatment_weeks
#3
aov.model <- aov(log_tle ~ trt_weeks, data = tle3)
summary(aov.model)
tukey_stats <- TukeyHSD(aov.model)
tukey_groups3 <- HSD.test(aov.model, "trt_weeks", group=TRUE)
tukey_groups3 <- tukey_groups3$groups
tukey_stats <- as.data.frame(tukey_stats$trt_weeks)
write.table(tukey_stats, file = "Tukey_tle3.txt", sep = "\t")
#6 
aov.model <- aov(log_tle ~ trt_weeks, data = tle6)
summary(aov.model)
tukey_stats <- TukeyHSD(aov.model)
tukey_groups6 <- HSD.test(aov.model, "trt_weeks", group=TRUE) #made this into pairwise matrix manually
tukey_groups6 <- tukey_groups6$groups
tukey_stats <- as.data.frame(tukey_stats$trt_weeks)
write.table(tukey_stats, file = "Tukey_tle6.txt", sep = "\t") 

tle3$trt_weeks <- factor(tle3$trt_weeks, levels = c("NT3","LA3", "HA3", "LN3", "HN3", "LP3", "HP3", "LC3", "HC3")) 
tle6$trt_weeks <- factor(tle6$trt_weeks, levels = c("NT6","LA6", "HA6", "LN6", "HN6", "LP6", "HP6", "LC6", "HC6")) 

#manually turned into triangular matrix in excel
tle_stats <- read.csv("tle_3_pairwise.csv", header = TRUE, row.names = 1)
myletters<-multcompLetters(tle_stats,compare="<=",threshold=0.05,Letters=letters)
myletters

colors = c('#848482', '#FCEDA2', '#FBD92E', '#E26060', '#BE0032', '#74A4D1', '#106EC8', '#70BD64','#1A9906')
colors_full = c('#848482', '#848482', '#FCEDA2', '#FCEDA2', '#FBD92E', '#FBD92E', '#E26060', '#E26060', '#BE0032', '#BE0032', '#74A4D1', '#74A4D1', '#106EC8', '#106EC8', '#70BD64', '#70BD64', '#1A9906', '#1A9906')

boxplot3 <- ggplot(tle3, aes(x=trt_weeks, y=tle, fill=trt_weeks)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = colors) +
  labs(x="Treatment Weeks", y = "Linear Extension (mm)") + 
  ylim(0,15) +
  stat_summary(geom = 'text', label = c("a","ab", "ab", "abc", "abc", "abcd", "bcd", "cd", "d"), fun = min, vjust = 2, size = 3.5) +
  theme(legend.position = "none")

ggsave(filename="~/Mote_nutrient_experiment/plots/tle_boxplot3.png", plot=boxplot3, device="png",  width = 5, height = 5, dpi=500)
ggsave(filename="~/Mote_nutrient_experiment/plots/tle_boxplot3.svg", plot=boxplot3, device="svg", dpi=500)

#I actually didn't make the 6 week stats into a triangle because only one comparison was significant - HC6 vs NT6
# tle_stats <- read.csv("tle_6_pairwise.csv", header = TRUE, row.names = 1)
# myletters<-multcompLetters(tle_stats,compare="<=",threshold=0.05,Letters=letters)
# myletters

boxplot6 <- ggplot(tle6, aes(x=trt_weeks, y=tle, fill=trt_weeks)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = colors) +
  ylim(0,15) +
  labs(x="Treatment Weeks", y = "Linear Extension (mm)") + 
  stat_summary(geom = 'text', label = c("a","ab","ab", "ab", "ab", "ab", "ab", "ab", "b"), fun = min, vjust = 2, size = 3.5) +
  theme(legend.position = "none")

ggsave(filename="~/Mote_nutrient_experiment/plots/tle_boxplot6.png", plot=boxplot6, device="png",  width = 5, height = 5, dpi=500)
ggsave(filename="~/Mote_nutrient_experiment/plots/tle_boxplot6.svg", plot=boxplot6, device="svg", dpi=500)

boxplot <- ggplot(tle, aes(x=trt_weeks, y=tle, fill=trt_weeks)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = colors_full) +
  labs(x="Treatment Weeks", y = "Linear Extension (mm)") + 
  theme(legend.position = "none")
ggsave(filename="~/Mote_nutrient_experiment/plots/tle_boxplot_full.png", plot=boxplot, device="png",  width = 6, height = 5, dpi=500)


#modeling relationship of abundance of ASV1 to treatment
model <- lm(rickabun ~ trt_weeks, data = tle)
summary(model)
aov.model <- aov(rickabun ~ trt_weeks, data = tle)
summary(aov.model)
TukeyHSD(aov.model)
#no significant relationship

#simple linear regression of tle to abundance of ASV1
#subset of only tle at 6 weeks
ggplot(tle6, aes(log_tle, rickabun)) +
  geom_point() +
  stat_smooth(method = lm) +
  xlab("log(Total Linear Extension (mm))") + ylab("Relative Abundance ASV1")

model <- lm(rickabun ~ log_tle, data = tle6)
summary(model)
