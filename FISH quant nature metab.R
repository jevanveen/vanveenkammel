library(tidyverse)
library(ggpubr)
library(cowplot)
library(scales)
library(pzfx)
library(Hmisc)
library(svglite)



## Import cellprofiler outputs ####


cytoplasm.rprm.sirna <- read.csv("~/Dropbox/JEV/CellProfiler/121619_rprm_siRNACytoplasm.csv")
cytoplasm.rprm.sirna <- read.csv("~/Desktop/121619_rprm_siRNACytoplasm.csv")
cytoplasm.rprm.sirna$Metadata_Mouse <- as.factor(cytoplasm.rprm.sirna$Metadata_Mouse)

nuc.rprm.sirna <- read.csv("~/Dropbox/JEV/CellProfiler/121619_rprm_siRNANuclei.csv")
nuc.rprm.sirna <- read.csv("~/Desktop/121619_rprm_siRNANuclei.csv")
nuc.rprm.sirna$Metadata_Mouse <- as.factor(nuc.rprm.sirna$Metadata_Mouse)


nuccyto.rprm.sirna <- bind_cols(cytoplasm.rprm.sirna, nuc.rprm.sirna)
levels(nuccyto.rprm.sirna$Metadata_Mouse)
nuccyto.rprm.sirna$KO <- ifelse(nuccyto.rprm.sirna$Metadata_Mouse == 1947 | nuccyto.rprm.sirna$Metadata_Mouse == 1950 | 
                                         nuccyto.rprm.sirna$Metadata_Mouse == 1951, "KO", "WT")
nuccyto.rprm.sirna$KO <- as.factor(nuccyto.rprm.sirna$KO)
nuccyto.rprm.sirna$KO <- ordered(nuccyto.rprm.sirna$KO, c("WT", "KO"))

ggplot(data = nuccyto.rprm.sirna, aes(x = KO, y = Intensity_MedianIntensity_proteinx, color = Metadata_Mouse)) +
  geom_point(position = "jitter")


#nkx kos in males quant


cytoplasm.rprm.nkxera.males <- read.csv("~/Dropbox/JEV/CellProfiler/rprmnkx1esr1_malesCytoplasm.csv")
cytoplasm.rprm.nkxera.males$Metadata_Mouse <- as.factor(cytoplasm.rprm.nkxera.males$Metadata_Mouse)

nuc.rprm.nkxera.males <- read.csv("~/Dropbox/JEV/CellProfiler/rprmnkx1esr1_malesNuclei.csv")
nuc.rprm.nkxera.males$Metadata_Mouse <- as.factor(nuc.rprm.nkxera.males$Metadata_Mouse)

levels(nuccyto.rprm.nkxera.males$Metadata_Mouse)
nuccyto.rprm.nkxera.males <- bind_cols(cytoplasm.rprm.nkxera.males, nuc.rprm.nkxera.males)
nuccyto.rprm.nkxera.males$KO <- ifelse(nuccyto.rprm.nkxera.males$Metadata_Mouse == 1923 | nuccyto.rprm.nkxera.males$Metadata_Mouse == 1926 | 
                                         nuccyto.rprm.nkxera.males$Metadata_Mouse == 1929, "KO", "WT")
nuccyto.rprm.nkxera.males$KO <- as.factor(nuccyto.rprm.nkxera.males$KO)
nuccyto.rprm.nkxera.males$KO <- ordered(nuccyto.rprm.nkxera.males$KO, c("WT", "KO"))

ggplot(data = nuccyto.rprm.nkxera.males, aes(x = Metadata_Mouse, y = Intensity_MedianIntensity_era)) +
  geom_point(position = "jitter")

#now tac1

cytoplasm.tac1.nkxera.males <- read.csv("~/Dropbox/JEV/CellProfiler/tac1nkx1esr1_malesCytoplasm.csv")
cytoplasm.tac1.nkxera.males$Metadata_Mouse <- as.factor(cytoplasm.tac1.nkxera.males$Metadata_Mouse)

nuc.tac1.nkxera.males <- read.csv("~/Dropbox/JEV/CellProfiler/tac1nkx1esr1_malesNuclei.csv")
nuc.tac1.nkxera.males$Metadata_Mouse <- as.factor(nuc.tac1.nkxera.males$Metadata_Mouse)

nuccyto.tac1.nkxera.males <- bind_cols(cytoplasm.tac1.nkxera.males, nuc.tac1.nkxera.males)
levels(nuccyto.tac1.nkxera.males$Metadata_Mouse)
nuccyto.tac1.nkxera.males$KO <- ifelse(nuccyto.tac1.nkxera.males$Metadata_Mouse == 1932 | nuccyto.tac1.nkxera.males$Metadata_Mouse == 1935 | 
                                         nuccyto.tac1.nkxera.males$Metadata_Mouse == 1937, "KO", "WT")
nuccyto.tac1.nkxera.males$KO <- as.factor(nuccyto.tac1.nkxera.males$KO)
nuccyto.tac1.nkxera.males$KO <- ordered(nuccyto.tac1.nkxera.males$KO, c("WT", "KO"))

ggplot(data = nuccyto.tac1.nkxera.males, aes(x = Metadata_Mouse, y = Intensity_MedianIntensity_era)) +
  geom_point(position = "jitter")



#rpr/tac1 coloc - DAPI stain was bad, have to identify each color independently

cyto.tacrprm <- read.csv("~/Dropbox/JEV/CellProfiler/rprmtac1Cytoplasm.csv")
cyto.tacrprm$Metadata_Mouse <- as.factor(cyto.tacrprm$Metadata_Mouse)
cyto.tacrprm$Metadata_Region <- as.factor(cyto.tacrprm$Metadata_Region)

cells.tacrprm <- read.csv("~/Dropbox/JEV/CellProfiler/rprmtac1Cells.csv")
cells.tacrprm$Metadata_Mouse <- as.factor(cells.tacrprm$Metadata_Mouse)
cells.tacrprm$Metadata_Region <- as.factor(cells.tacrprm$Metadata_Region)

levels(cyto.tacrprm$Metadata_Mouse)

#rprm/era coloc

cytoplasm.rprm.era <- read.csv("~/Dropbox/JEV/CellProfiler/rprmeracolocCytoplasm.csv")
cytoplasm.rprm.era$Metadata_Mouse <- as.factor(cytoplasm.rprm.era$Metadata_Mouse)
cytoplasm.rprm.era$Metadata_Region <- as.factor(cytoplasm.rprm.era$Metadata_Region)

nuc.rprm.era <- read.csv("~/Dropbox/JEV/CellProfiler/rprmeracolocNuclei.csv")
nuc.rprm.era$Metadata_Mouse <- as.factor(nuc.rprm.era$Metadata_Mouse)
nuc.rprm.era$Metadata_Region <- as.factor(nuc.rprm.era$Metadata_Region)

nuccyto.rprm.era <- bind_cols(cytoplasm.rprm.era, nuc.rprm.era)

nuccyto.rprm.era$MF <- ifelse(nuccyto.rprm.era$Metadata_Mouse == 182 | nuccyto.rprm.era$Metadata_Mouse == 183 | 
                                nuccyto.rprm.era$Metadata_Mouse == 184 | nuccyto.rprm.era$Metadata_Mouse == 185, "Female", "Male")

males <- filter(nuccyto.rprm.era, MF == "Male")
females <- filter(nuccyto.rprm.era, MF == "Female")
levels(nuccyto.rprm.era$Metadata_Mouse)



#tac1/era coloc

cytoplasm.tac1.era <- read.csv("~/Dropbox/JEV/CellProfiler/tac1eracolocCytoplasm.csv")
cytoplasm.tac1.era$Metadata_Mouse <- as.factor(cytoplasm.tac1.era$Metadata_Mouse)
cytoplasm.tac1.era$Metadata_Region <- as.factor(cytoplasm.tac1.era$Metadata_Region)

nuc.tac1.era <- read.csv("~/Dropbox/JEV/CellProfiler/tac1eracolocNuclei.csv")
nuc.tac1.era$Metadata_Mouse <- as.factor(nuc.tac1.era$Metadata_Mouse)
nuc.tac1.era$Metadata_Region <- as.factor(nuc.tac1.era$Metadata_Region)

nuccyto.tac1.era <- bind_cols(cytoplasm.tac1.era, nuc.tac1.era)



nuccyto.tac1.era$MF <- ifelse(nuccyto.tac1.era$Metadata_Mouse == 186 | nuccyto.tac1.era$Metadata_Mouse == 187 | 
                                nuccyto.tac1.era$Metadata_Mouse == 188 | nuccyto.tac1.era$Metadata_Mouse == 189 |
                                nuccyto.tac1.era$Metadata_Mouse == 1810, "Female", "Male")
males <- filter(nuccyto.tac1.era, MF == "Male")
females <- filter(nuccyto.tac1.era, MF == "Female")

levels(males$Metadata_Mouse)
levels(females$Metadata_Mouse)

#sst/era coloc

cytoplasm.sst.era <- read.csv("~/Dropbox/JEV/CellProfiler/ssteracolocCytoplasm.csv")
cytoplasm.sst.era$Metadata_Mouse <- as.factor(cytoplasm.sst.era$Metadata_Mouse)
cytoplasm.sst.era$Metadata_Region <- as.factor(cytoplasm.sst.era$Metadata_Region)

nuc.sst.era <- read.csv("~/Dropbox/JEV/CellProfiler/ssteracolocNuclei.csv")
nuc.sst.era$Metadata_Mouse <- as.factor(nuc.sst.era$Metadata_Mouse)
nuc.sst.era$Metadata_Region <- as.factor(nuc.sst.era$Metadata_Region)

nuccyto.sst.era <- bind_cols(cytoplasm.sst.era, nuc.sst.era)

nuccyto.sst.era.erapos <- filter(nuccyto.sst.era, Intensity_MedianIntensity_era >= (mean(Intensity_MedianIntensity_era) + 0* sd(Intensity_MedianIntensity_era)))

#pdyn/era coloc

cytoplasm.pdyn.era <- read.csv("~/Dropbox/JEV/CellProfiler/pdyneracolocCytoplasm.csv")
cytoplasm.pdyn.era$Metadata_Mouse <- as.factor(cytoplasm.pdyn.era$Metadata_Mouse)
cytoplasm.pdyn.era$Metadata_Region <- as.factor(cytoplasm.pdyn.era$Metadata_Region)

nuc.pdyn.era <- read.csv("~/Dropbox/JEV/CellProfiler/pdyneracolocNuclei.csv")
nuc.pdyn.era$Metadata_Mouse <- as.factor(nuc.pdyn.era$Metadata_Mouse)
nuc.pdyn.era$Metadata_Region <- as.factor(nuc.pdyn.era$Metadata_Region)

nuccyto.pdyn.era <- bind_cols(cytoplasm.pdyn.era, nuc.pdyn.era)


nuccyto.pdyn.era$MF <- ifelse(nuccyto.pdyn.era$Metadata_Mouse == 1816 | nuccyto.pdyn.era$Metadata_Mouse == 1817 | 
                                nuccyto.pdyn.era$Metadata_Mouse == 1818 | nuccyto.pdyn.era$Metadata_Mouse == 1819 |
                                nuccyto.pdyn.era$Metadata_Mouse == 1820, "Female", "Male")

males <- filter(nuccyto.pdyn.era, MF == "Male")
females <- filter(nuccyto.pdyn.era, MF == "Female")

levels(nuccyto.pdyn.era$Metadata_Mouse)



nuccyto.pdyn.era.erapos <- filter(nuccyto.pdyn.era, Intensity_MedianIntensity_era >= mean(filter(nuccyto.pdyn.era, MF == "Female")$Intensity_MedianIntensity_era))

ggplot(data = nuccyto.pdyn.era, aes(x = Metadata_Mouse, y = Intensity_MedianIntensity_era)) +
  geom_violin() +
  theme_classic()

#rprm measurements in sf1:erko


cytoplasm.rprm.sf1era <- read.csv("~/Dropbox/JEV/CellProfiler/rprmsf1esr1v3Cytoplasm.csv")
cytoplasm.rprm.sf1era$Metadata_Mouse <- as.factor(cytoplasm.rprm.sf1era$Metadata_Mouse)
cytoplasm.rprm.sf1era$Metadata_Region <- as.factor(cytoplasm.rprm.sf1era$Metadata_Region)

nuc.rprm.sf1era <- read.csv("~/Dropbox/JEV/CellProfiler/rprmsf1esr1v3Nuclei.csv")
nuc.rprm.sf1era$Metadata_Mouse <- as.factor(nuc.rprm.sf1era$Metadata_Mouse)
nuc.rprm.sf1era$Metadata_Region <- as.factor(nuc.rprm.sf1era$Metadata_Region)


nuccyto.rprm.sf1era <- bind_cols(cytoplasm.rprm.sf1era, nuc.rprm.sf1era)
nuccyto.rprm.sf1era$KO <- ifelse(nuccyto.rprm.sf1era$Metadata_Mouse == 1863 | nuccyto.rprm.sf1era$Metadata_Mouse == 1865 | 
                                   nuccyto.rprm.sf1era$Metadata_Mouse == 1871 | nuccyto.rprm.sf1era$Metadata_Mouse == 1872, "KO", "WT")
nuccyto.rprm.sf1era$KO <- as.factor(nuccyto.rprm.sf1era$KO)
nuccyto.rprm.sf1era$KO <- ordered(nuccyto.rprm.sf1era$KO, c("WT", "KO"))

nuccyto.rprm.sf1era.vmh <- filter(nuccyto.rprm.sf1era, Metadata_Region == "vmh")
nuccyto.rprm.sf1era.arc <- filter(nuccyto.rprm.sf1era, Metadata_Region == "arc")



#rprm measurements in nkx:erko

cytoplasm.rprm.nkxera <- read.csv("~/Dropbox/JEV/CellProfiler/rprmnkxesr1v2Cytoplasm.csv")
cytoplasm.rprm.nkxera$Metadata_Mouse <- as.factor(cytoplasm.rprm.nkxera$Metadata_Mouse)
cytoplasm.rprm.nkxera$Metadata_Region <- as.factor(cytoplasm.rprm.nkxera$Metadata_Region)

nuc.rprm.nkxera <- read.csv("~/Dropbox/JEV/CellProfiler/rprmnkxesr1v2Nuclei.csv")
nuc.rprm.nkxera$Metadata_Mouse <- as.factor(nuc.rprm.nkxera$Metadata_Mouse)
nuc.rprm.nkxera$Metadata_Region <- as.factor(nuc.rprm.nkxera$Metadata_Region)


nuccyto.rprm.nkxera <- bind_cols(cytoplasm.rprm.nkxera, nuc.rprm.nkxera)
nuccyto.rprm.nkxera$KO <- ifelse(nuccyto.rprm.nkxera$Metadata_Mouse == 1858 | nuccyto.rprm.nkxera$Metadata_Mouse == 1860 | 
                                   nuccyto.rprm.nkxera$Metadata_Mouse == 1862, "KO", "WT")
nuccyto.rprm.nkxera$KO <- as.factor(nuccyto.rprm.nkxera$KO)
nuccyto.rprm.nkxera$KO <- ordered(nuccyto.rprm.nkxera$KO, c("WT", "KO"))

nuccyto.rprm.nkxera.vmh <- filter(nuccyto.rprm.nkxera, Metadata_Region == "vmh")
nuccyto.rprm.nkxera.arc <- filter(nuccyto.rprm.nkxera, Metadata_Region == "arc")



#tac1 measurements in nkx:esr1f animals

cytoplasm.tac1.nkxera <- read.csv("~/Dropbox/JEV/CellProfiler/tac1nkxesr1Cytoplasm.csv")
cytoplasm.tac1.nkxera$Metadata_Mouse <- as.factor(cytoplasm.tac1.nkxera$Metadata_Mouse)
cytoplasm.tac1.nkxera$Metadata_Region <- as.factor(cytoplasm.tac1.nkxera$Metadata_Region)

nuc.tac1.nkxera <- read.csv("~/Dropbox/JEV/CellProfiler/tac1nkxesr1Nuclei.csv")
nuc.tac1.nkxera$Metadata_Mouse <- as.factor(nuc.tac1.nkxera$Metadata_Mouse)
nuc.tac1.nkxera$Metadata_Region <- as.factor(nuc.tac1.nkxera$Metadata_Region)


nuccyto.tac1.nkxera <- bind_cols(cytoplasm.tac1.nkxera, nuc.tac1.nkxera)
nuccyto.tac1.nkxera$KO <- ifelse(nuccyto.tac1.nkxera$Metadata_Mouse == 1858 | nuccyto.tac1.nkxera$Metadata_Mouse == 1860 | 
                                   nuccyto.tac1.nkxera$Metadata_Mouse == 1862 | nuccyto.tac1.nkxera$Metadata_Mouse == 1875 | nuccyto.tac1.nkxera$Metadata_Mouse == 1876 | nuccyto.tac1.nkxera$Metadata_Mouse == 1877, "KO", "WT")
nuccyto.tac1.nkxera$KO <- as.factor(nuccyto.tac1.nkxera$KO)
nuccyto.tac1.nkxera$KO <- ordered(nuccyto.tac1.nkxera$KO, c("WT", "KO"))

nuccyto.tac1.nkxera.vmh <- filter(nuccyto.tac1.nkxera, Metadata_Region == "vmh")
nuccyto.tac1.nkxera.arc <- filter(nuccyto.tac1.nkxera, Metadata_Region == "arc")


#tac1 measurements in sf1:esr1f animals

cytoplasm.tac1.sf1era <- read.csv("~/Dropbox/JEV/CellProfiler/tac1sf1esr1Cytoplasm.csv")
cytoplasm.tac1.sf1era$Metadata_Mouse <- as.factor(cytoplasm.tac1.sf1era$Metadata_Mouse)
cytoplasm.tac1.sf1era$Metadata_Region <- as.factor(cytoplasm.tac1.sf1era$Metadata_Region)

nuc.tac1.sf1era <- read.csv("~/Dropbox/JEV/CellProfiler/tac1sf1esr1Nuclei.csv")
nuc.tac1.sf1era$Metadata_Mouse <- as.factor(nuc.tac1.sf1era$Metadata_Mouse)
nuc.tac1.sf1era$Metadata_Region <- as.factor(nuc.tac1.sf1era$Metadata_Region)


nuccyto.tac1.sf1era <- bind_cols(cytoplasm.tac1.sf1era, nuc.tac1.sf1era)
nuccyto.tac1.sf1era$KO <- ifelse(nuccyto.tac1.sf1era$Metadata_Mouse == 1879 | nuccyto.tac1.sf1era$Metadata_Mouse == 1880 | nuccyto.tac1.sf1era$Metadata_Mouse == 1886, "KO", "WT")
nuccyto.tac1.sf1era$KO <- as.factor(nuccyto.tac1.sf1era$KO)
nuccyto.tac1.sf1era$KO <- ordered(nuccyto.tac1.sf1era$KO, c("WT", "KO"))

nuccyto.tac1.sf1era.vmh <- filter(nuccyto.tac1.sf1era, Metadata_Region == "vmh")
nuccyto.tac1.sf1era.arc <- filter(nuccyto.tac1.sf1era, Metadata_Region == "arc")

#add WT animals from sf1 experiments to nkx experiments

nuccyto.tac1.sf1era.vmh.wt <- filter(nuccyto.tac1.sf1era.vmh, KO == "WT")
nuccyto.rprm.sf1era.vmh.wt <- filter(nuccyto.rprm.sf1era.vmh, KO == "WT")

nuccyto.tac1.nkxera.vmh.combined <- bind_rows(nuccyto.tac1.nkxera.vmh, nuccyto.tac1.sf1era.vmh.wt)
nuccyto.rprm.nkxera.vmh.combined <- bind_rows(nuccyto.rprm.nkxera.vmh, nuccyto.rprm.sf1era.vmh.wt)

ggplot(data = nuccyto.rprm.nkxera.vmh.combined, aes(x = Metadata_Mouse, y = Intensity_MedianIntensity_era)) +
  geom_point(position = "jitter")



z <- nuccyto.rprm
z <- nuccyto.rprm.sf1era
z <- nuccyto.rprm.sf1era.vmh
z <- nuccyto.rprm.sf1era.arc
z <- nuccyto.rprm.nkxera
z <- nuccyto.rprm.nkxera.vmh
z <- nuccyto.rprm.nkxera.arc
z <- nuccyto.tac1.sf1era
z <- nuccyto.tac1.sf1era.vmh
z <- nuccyto.tac1.sf1era.arc
z <- nuccyto.tac1.nkxera.vmh
z <- nuccyto.tac1.nkxera.arc
z <- nuccyto.rprm.sirna
z <- nuccyto.rprm.sirna.jev
z <- nuccyto.rprm.sirna.scv
z <- nuccyto.rprm.sirna.scv.erapos
z <- nuccyto.rprm.era
z <- nuccyto.tac1.era
z <- nuccyto.sst.era
z <- nuccyto.pdyn.era
z <- cells.tacrprm
z <- nuccyto.tac1.nkxera.males
z <- nuccyto.rprm.nkxera.males
z <- nuccyto.tac1.nkxera.vmh.combined
z <- nuccyto.rprm.nkxera.vmh.combined


### here begin stats ####

#1n per mouse - use

agg.nuc.era <- z %>% group_by(Metadata_Mouse, KO) %>% summarise(mean_intensity = mean(Intensity_MeanIntensity_era), median_intensity = mean(Intensity_MedianIntensity_era))
agg.cyto.proteinx <- z %>% group_by(Metadata_Mouse, KO) %>% summarise(mean_intensity = mean(Intensity_MeanIntensity_proteinx), median_intensity = mean(Intensity_MedianIntensity_proteinx))


agg.nuc.era <- z %>% group_by(Metadata_Mouse, MF) %>% summarise(mean_intensity = mean(Intensity_MeanIntensity_era), median_intensity = mean(Intensity_MedianIntensity_era))
agg.cyto.proteinx <- z %>% group_by(Metadata_Mouse, MF) %>% summarise(mean_intensity = mean(Intensity_MeanIntensity_proteinx), median_intensity = mean(Intensity_MedianIntensity_proteinx))



#test for normality

shapiro.test(agg.cyto.proteinx$median_intensity)
shapiro.test(agg.nuc.era$median_intensity)


# KO stats

t.test(median_intensity ~ KO, data = agg.cyto.proteinx)
t.test(median_intensity ~ KO, data = agg.nuc.era)

wilcox.test(median_intensity ~ KO, data = agg.nuc.era)
wilcox.test(median_intensity ~ KO, data = agg.cyto.proteinx)


#mvf stats


t.test(median_intensity ~ MF, data = agg.cyto.proteinx)
t.test(median_intensity ~ MF, data = agg.nuc.era)

wilcox.test(median_intensity ~ MF, data = agg.nuc.era)
wilcox.test(median_intensity ~ MF, data = agg.cyto.proteinx)


## Graphs here ####
#first set only work for tac/rprm coloc



ggplot(z, aes(x = Intensity_MedianIntensity_rprm, y = Intensity_MedianIntensity_tac1, color = Metadata_Mouse)) +
  geom_point() +
  geom_hline(yintercept = .075, linetype="dashed") +
  geom_vline(xintercept = .025, linetype="dashed") +
  ylim(0,.5) +
  xlim(0,.15)+
  labs(title = "tac1")


ggplot(z, aes(x = Intensity_MedianIntensity_rprm, y = Intensity_MedianIntensity_tac1)) +
  stat_binhex(bins = 40) +
  geom_hline(yintercept = .075, linetype="dashed") +
  geom_vline(xintercept = .025, linetype="dashed") +
  ylim(0,.5) +
  xlim(0,.15)+
  labs(title = "tac1")


#general graphs


p1 <- ggplot(z, aes(x = Metadata_Mouse, y = Intensity_MedianIntensity_proteinx)) + 
  geom_point(position = "jitter") +
  stat_summary(fun.data=mean_sdl,
               geom="pointrange", color="red") +
  coord_cartesian(ylim = c(0,.5))

p2 <- ggplot(z, aes(x = Metadata_Mouse, y = Intensity_MedianIntensity_era)) + 
  geom_point(position = "jitter") +
  stat_summary(fun.data=mean_sdl,
               geom="pointrange", color="red") +
  coord_cartesian(ylim = c(0,.5))

plot_grid(p1, p2)


p1 <- ggplot(z, aes(x = MF, y = Intensity_MedianIntensity_proteinx, color = Metadata_Mouse)) + 
  geom_point(position = "jitter") +
  stat_summary(fun.data=mean_sdl,
               geom="pointrange", color="black", size = 1) +
  coord_cartesian(ylim = c(0,.5)) + 
  scale_color_discrete() + 
  labs(color = "Mouse") + 
  labs(title = "Rprm in ARC, nkx:erKO, p = .9143")



ggplot(z, aes(x = KO, y = Intensity_MedianIntensity_era, color = Metadata_Mouse)) + 
  geom_point(position = "jitter") +
  stat_summary(fun.data=mean_sdl,
               geom="pointrange", color="black", size = 1) +
  coord_cartesian(ylim = c(0,.4)) + 
  scale_color_viridis_d(alpha = .5) + 
  labs(color = "Mouse", title = "Era in VMH nkx:eraKO, p = .0095")

plot_grid(p1,p2)

ggsave("~/Desktop/NatureFigs/nkxKOrprmarc.tiff", plot = p1, width = 5, height = 8, dpi = 320, units = "in")
ggsave("~/Desktop/NatureFigs/nkxKOeravmh.tiff", plot = p2, width = 5, height = 8, dpi = 320, units = "in")



p1 <- ggplot(z, aes(x = KO, y = Intensity_MedianIntensity_era)) +
  geom_point(position = "jitter", color = "darkgray", alpha = .7, size = 1) +
  geom_point(data = agg.nuc.era, aes(x = KO, y = median_intensity, fill = Metadata_Mouse), size = 3, alpha = .7, shape = 23, color = "black") +
  theme_classic() +
  labs(fill = "Mouse", title = "t test p = .0004") + 
  ylim(0,.2)


z <- nuccyto.tac1.nkxera.males
z <- nuccyto.rprm.nkxera.males
z <- nuccyto.tac1.nkxera.vmh.combined
z <- nuccyto.rprm.nkxera.vmh.combined
agg.cyto.proteinx <- z %>% group_by(Metadata_Mouse, KO) %>% summarise(mean_intensity = mean(Intensity_MeanIntensity_proteinx), median_intensity = mean(Intensity_MedianIntensity_proteinx))


#this is the plot seen in figures 5 and 7: all cells in gray, per animal means as larger diamonds. 

p1 <- ggplot(z, aes(x = KO, y = Intensity_MedianIntensity_proteinx)) +
  geom_point(position = "jitter", color = "darkgray", alpha = .7, size = 1) +
  geom_point(position = position_jitter(width = 0.1, height = 0), data = agg.cyto.proteinx, aes(x = KO, y = median_intensity, fill = Metadata_Mouse), size = 2.5, alpha = .7, shape = 23, color = "black") +
  theme_classic() + 
  labs(fill = "Mouse", title = "t-test pval = .01548") +
  theme(legend.position = "none") 




plot_grid(p1, p2)
save_plot(filename = "~/box/scHypo paper/nkxkofemalesrprm.tiff", plot = p1, base_asp = .55, base_height = 3)
save_plot(filename = "~/box/scHypo paper/nkxkomalesrprm-rprm.tiff", plot = p2, base_asp = .8)

p2 <- ggplot(z, aes(x = Intensity_MedianIntensity_era, y = Intensity_MedianIntensity_proteinx, color = KO)) + 
  geom_point(alpha = .4) + 
  theme_classic() + 
  scale_color_manual(values = c("red", "blue"))
plot_grid(p1, p2)

z <- nuccyto.tac1.nkxera.males
z <- nuccyto.rprm.nkxera.males
z <- nuccyto.tac1.nkxera.vmh.combined
z <- nuccyto.rprm.nkxera.vmh.combined
z <- filter(z, Intensity_MedianIntensity_era > .01 | Intensity_MedianIntensity_proteinx > .05)
agg.nuc.cyto <- z %>% group_by(Metadata_Mouse, KO) %>% summarise(era_intensity = mean(Intensity_MedianIntensity_era), proteinx_intensity = mean(Intensity_MedianIntensity_proteinx))


## Start two axis graphs ####

z <- nuccyto.rprm.era
z <- nuccyto.tac1.era
z <- nuccyto.sst.era
z <- nuccyto.pdyn.era

# line based on normal distribution and observed % pos in images for rough cutoff
# era in females roughly 50% positive = 0sd
# rprm in females roughly 32% positive = .5sd
# tac1 in females roughly 50% positive = 0sd
# sst in females 10% = 1.5sd
# pdyn in males 
# now that I am combining male and female in same analyses, will have to set lines based on either male or female staining separately

vmhvline = (mean(z$Intensity_MedianIntensity_era) + 0*sd(z$Intensity_MedianIntensity_era))
vmhhline = (mean(z$Intensity_MedianIntensity_proteinx) + 1.5*sd(z$Intensity_MedianIntensity_proteinx))
vmhvline = (mean(z$Intensity_MedianIntensity_era) + 0*sd(z$Intensity_MedianIntensity_era))


vmhvline = (mean(z$Intensity_MedianIntensity_era) + 0*sd(z$Intensity_MedianIntensity_era))
vmhvline = (mean(z$Intensity_MedianIntensity_era) + .5*sd(z$Intensity_MedianIntensity_era))
vmhvline = (mean(z$Intensity_MedianIntensity_era) + 1*sd(z$Intensity_MedianIntensity_era))

vmhhline = (mean(filter(z, MF == "Female")$Intensity_MedianIntensity_proteinx) + 0*sd(filter(z, MF == "Female")$Intensity_MedianIntensity_proteinx))
vmhvline = (mean(filter(z, MF == "Female")$Intensity_MedianIntensity_era) + 0*sd(filter(z, MF == "Female")$Intensity_MedianIntensity_era))

vmhhline = (mean(filter(z, MF == "Male")$Intensity_MedianIntensity_proteinx) + 0*sd(filter(z, MF == "Male")$Intensity_MedianIntensity_proteinx))
vmhvline = (mean(filter(z, MF == "Female")$Intensity_MedianIntensity_era) + 0*sd(filter(z, MF == "Female")$Intensity_MedianIntensity_era))

#rprm tac1 doubles
vmhvline = (mean(z$Intensity_MedianIntensity_rprm) + 0*sd(z$Intensity_MedianIntensity_rprm))
vmhhline = (mean(z$Intensity_MedianIntensity_tac1) + 0*sd(z$Intensity_MedianIntensity_tac1))



#counts for rprm/tac1 
q1 <- count(z, z$Intensity_MeanIntensity_era < vmhvline & 
              z$Intensity_MeanIntensity_proteinx >= vmhhline)
q2 <- count(z, z$Intensity_MeanIntensity_era >= vmhvline & 
              z$Intensity_MeanIntensity_proteinx >= vmhhline)
q3 <- count(z, z$Intensity_MeanIntensity_era >= vmhvline & 
              z$Intensity_MeanIntensity_proteinx < vmhhline)
q4 <- count(z, z$Intensity_MeanIntensity_era < vmhvline & 
              z$Intensity_MeanIntensity_proteinx < vmhhline)


#separate counts for m and f

q1 <- count(filter(z, MF == "Female"), filter(z, MF == "Female")$Intensity_MedianIntensity_era < vmhvline & 
              filter(z, MF == "Female")$Intensity_MedianIntensity_proteinx >= vmhhline)
q2 <- count(filter(z, MF == "Female"), filter(z, MF == "Female")$Intensity_MedianIntensity_era >= vmhvline & 
              filter(z, MF == "Female")$Intensity_MedianIntensity_proteinx >= vmhhline)
q3 <- count(filter(z, MF == "Female"), filter(z, MF == "Female")$Intensity_MedianIntensity_era >= vmhvline & 
              filter(z, MF == "Female")$Intensity_MedianIntensity_proteinx < vmhhline)
q4 <- count(filter(z, MF == "Female"), filter(z, MF == "Female")$Intensity_MedianIntensity_era < vmhvline & 
              filter(z, MF == "Female")$Intensity_MedianIntensity_proteinx < vmhhline)
q5 <- count(filter(z, MF == "Male"), filter(z, MF == "Male")$Intensity_MedianIntensity_era < vmhvline & 
              filter(z, MF == "Male")$Intensity_MedianIntensity_proteinx >= vmhhline)
q6 <- count(filter(z, MF == "Male"), filter(z, MF == "Male")$Intensity_MedianIntensity_era >= vmhvline & 
              filter(z, MF == "Male")$Intensity_MedianIntensity_proteinx >= vmhhline)
q7 <- count(filter(z, MF == "Male"), filter(z, MF == "Male")$Intensity_MedianIntensity_era >= vmhvline & 
              filter(z, MF == "Male")$Intensity_MedianIntensity_proteinx < vmhhline)
q8 <- count(filter(z, MF == "Male"), filter(z, MF == "Male")$Intensity_MedianIntensity_era < vmhvline & 
              filter(z, MF == "Male")$Intensity_MedianIntensity_proteinx < vmhhline)



# rprm tac1 counts. note different q1, etc. definitions than previous. Because I switched the x and y of rprm, this is out of order. redo if you like, but its right in prism
q1 <- count(cells.tacrprm, z$Intensity_MedianIntensity_rprm <= .03 & 
              z$Intensity_MedianIntensity_tac1 <= .075)
q2 <- count(cells.tacrprm, z$Intensity_MedianIntensity_rprm <= .03 & 
              z$Intensity_MedianIntensity_tac1 > .075)
q3 <- count(cells.tacrprm, z$Intensity_MedianIntensity_rprm > .03 & 
              z$Intensity_MedianIntensity_tac1 > .075)
q4 <- count(cells.tacrprm, z$Intensity_MedianIntensity_rprm > .03 & 
              z$Intensity_MedianIntensity_tac1 <= .075)


#to make two axis graphs in figure 4
p1 <- ggplot(data = z, aes(x = Intensity_MedianIntensity_era, y = Intensity_MedianIntensity_proteinx, color = MF)) + 
  geom_point(size = .4, alpha = .5, show.legend = F) +
  geom_hline(yintercept = vmhhline, linetype="dashed", color = "black") +
  geom_vline(xintercept = vmhvline, linetype="dashed", color = "black") +
  scale_color_manual(values = c("hotpink", "royalblue")) +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) 

#sst only: figure supp. 4
p1 <- ggplot(data = z, aes(x = Intensity_MedianIntensity_era, y = Intensity_MedianIntensity_proteinx)) + 
  geom_point(size = .4, alpha = .5, show.legend = F, color = "hotpink") +
  geom_hline(yintercept = vmhhline, linetype="dashed", color = "black") +
  geom_vline(xintercept = vmhvline, linetype="dashed", color = "black") +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  ylim(0,.5)

save_plot(filename = "~/box/scHypo paper/sstera.tiff", plot = p1, base_asp = 1, base_height = 2)



p1 <- ggplot(data = z, aes(x = Intensity_MedianIntensity_tac1, y = Intensity_MedianIntensity_rprm)) + 
  geom_point(size = 1, alpha = .5, color = "hotpink", show.legend = F) +
  geom_vline(xintercept = .075, linetype="dashed", color = "black") +
  geom_hline(yintercept = .03, linetype="dashed", color = "black") +
  xlim(0,.5) +
  ylim(0,.15) + 
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) 

save_plot(filename = "~/box/scHypo paper/rprmtac1.tiff", plot = p1, base_asp = 1, base_height = 2)


