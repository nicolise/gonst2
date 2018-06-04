#03-20-17
#Nicole E Soltis
#Analysis for Botrytis x GONST2 in Arabidopsis project
#----------------------------------------------------------
rm(list=ls())
setwd("~/Projects/gonst2/data")
ModDat <- read.csv("AllResults.csv")
PlantKey <- read.csv("PlantKey.csv")
ModDat <- ModDat[,c(1:10, 153)]
#remove rows with NA as isolate
names(ModDat)
names(PlantKey)
#merge in Plant Genotype by "Flat_Col_Row"
PlantKey <- PlantKey[,c("Genotype", "Tray", "Flat_Col_Row")]
ModDat <- merge(PlantKey,ModDat,by="Flat_Col_Row")
ModDat <- ModDat[!(is.na(ModDat$Isolate)),]
#add a variable for scaled lesion size
ModDat$Lesion.Size.cm <- ModDat$Lesion.Size / (ModDat$PixelsPerCM^2)
plot(ModDat$Lesion.Size.cm ~ ModDat$Isolate)
plot(ModDat$Lesion.Size.cm ~ ModDat$Genotype)

#extract mean +-SE for figure 5
SummDat <- ddply(ModDat, c("Genotype", "Time"), summarise,
                   mLS   = mean(Lesion.Size.cm),
                   sdLS = sd(Lesion.Size.cm),
                   nLS = length(Lesion.Size.cm),
                 seLS = sdLS/(nLS^0.5))

SummDat.2 <- ddply(ModDat, c("Genotype"), summarise,
                 mLS   = mean(Lesion.Size.cm),
                 sdLS = sd(Lesion.Size.cm),
                 nLS = length(Lesion.Size.cm),
                 seLS = sdLS/(nLS^0.5))

#linear model
names(ModDat)
ModDat <- ModDat[ModDat$Isolate!="Control",]
ModDat$Iso.Geno <- paste(ModDat$Isolate, ModDat$Genotype, sep="_")
ModDat.96 <- ModDat[ModDat$Time=="96",]
ModDat.72 <- ModDat[ModDat$Time=="72",]

#how many observations per genotype combination?
ModDat.72.G1 <- ModDat.72[ModDat.72$Tray=="G1",]
table(ModDat.72.G1$Iso.Geno)
ModDat.72.G2 <- ModDat.72[ModDat.72$Tray=="G2",]
table(ModDat.72.G2$Iso.Geno)
table(ModDat.72$Iso.Geno)
table(ModDat.96$Iso.Geno)

#fixed effects models
mymod.72 <- lm(Lesion.Size.cm ~ Isolate * Genotype + Tray + Tray/Flat + Tray:Genotype + Tray:Isolate, data=ModDat.72)
anova(mymod.72)

mymod.72.s <- lm(Lesion.Size.cm ~ Isolate * Genotype + Tray + Tray/Flat, data=ModDat.72)
anova(mymod.72.s)

mymod.96 <- lm(Lesion.Size.cm ~ Isolate * Genotype + Tray + Tray/Flat, data=ModDat.96)
anova(mymod.96)

mymod.96.s <- lm(Lesion.Size.cm ~ Isolate * Genotype + Tray + Tray/Flat, data=ModDat.96)
anova(mymod.96.s)

mymod <- lm(Lesion.Size.cm ~ Isolate * Genotype + Tray + Tray/Flat + Time, data=ModDat)
anova(mymod)

unique(ModDat$Isolate)
ModDat.72 <- ModDat[ModDat$Time=="72",]

write.csv(ModDat, "AllResults_NESgonst2.csv")

#mixed effects models (+ rand)
library("lme4"); library(car); library(lmerTest)
mymod.72.s.mm <- lmer(Lesion.Size.cm ~ Isolate * Genotype + (1|Tray) + (1|Tray:Flat), data=ModDat.72)
anova(mymod.72.s.mm)
Anova(mymod.72.s.mm, type=2)
ranova(mymod.72.s.mm)

mymod.96.s.mm <- lmer(Lesion.Size.cm ~ Isolate * Genotype + (1|Tray) + (1|Tray:Flat), data=ModDat.96)
anova(mymod.96.s.mm)
Anova(mymod.96.s.mm, type=2)
ranova(mymod.96.s.mm)

mymod.s.mm <- lmer(Lesion.Size.cm ~ Isolate * Genotype + (1|Tray) + (1|Tray:Flat) + Time, data=ModDat)
anova(mymod.s.mm)
Anova(mymod.s.mm, type=2)
ranova(mymod.s.mm)


