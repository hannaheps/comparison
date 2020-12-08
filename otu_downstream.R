##97% clustered OTU Statistical Analyses for comparison MS##

library(qiime2R)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyverse)
library(gridExtra)
library(vegan)
library(decontam)
library(RColorBrewer)

#Upload all files from qiime2 into phyloseq
phy.97 <- qza_to_phyloseq("~/Desktop/comparison/97_OTUS/table_97_filtered.qza", 
                          "~/Desktop/comparison/97_OTUS/rooted-tree_all_97.qza", 
                          "~/Desktop/comparison/97_OTUS/tax_97.qza", 
                          "~/Desktop/comparison/metadata_comp.txt")

#Check the taxonomic classification is correct
rank_names(phy.97)
tax_table(phy.97)


phy.97 <- subset_taxa(phy.97, Kingdom != "D_0__Eukaryota")
phy.97 <- subset_taxa(phy.97, Kingdom != "Eukaryota")
phy.97<- subset_taxa(phy.97, Kingdom != "Unassigned")
#phy.97 <- subset_taxa(phy.97, Family != "Mitochondria")
phy.97 <- subset_samples(phy.97, shared == "yes")

##Remove Mitochondria taxa:
#first check which taxa are not bacteria by subsetting Family:Mitochondria and blasting resulting seqs
mito.97 <- subset_taxa(phy.97, Family == "Mitochondria")
mito_taxa.97 <- rownames(otu_table(mito.97))
write.csv(mito_taxa.97, "~/Desktop/comparison/mito_97.csv")
#Here check ASVs on BLAST for mitochondrial vs. bacterial reads
#Keep mitochondrial reads in CSV and re-upload as "bad taxa" to remove
bad.taxa <- read.csv("~/Desktop/comparison/97_OTUS/bad_mito_97.csv", header = FALSE)
bad.taxa <- levels(bad.taxa$V1)
#remove them from the dataset
all.taxa <- taxa_names(phy.97)
all.taxa <- all.taxa[!(all.taxa %in% bad.taxa)]
physeq.nm <- prune_taxa(all.taxa, phy.97)

sample.data <- as(sample_data(physeq.nm), "data.frame")

#Remove contaminants
#Inspect library size
sample.data$LibrarySize <- sample_sums(physeq.nm)
sample.data <- sample.data[order(sample.data$LibrarySize),]
sample.data$Index <- seq(nrow(sample.data))
ggplot(data = sample.data, aes(x=Index, y=LibrarySize, color = sample_type)) +
  geom_point()

#Next check for contaminants using prevalence and threshold 0.5 (more conservative)
sample_data(physeq.nm)$is.neg <- sample_data(physeq.nm)$sample_type == "control"
contamdf.prev <- isContaminant(physeq.nm, method = "prevalence", neg = "is.neg", threshold = 0.5)
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

#Remove contaminants
physeq.noncont <- prune_taxa(!contamdf.prev$contaminant, physeq.nm)

#Check final numbesr
tax_table(physeq.noncont)  #Total ASVs:3900
sample_data(physeq.noncont) #Total samples:80
#Count Sequences
totalreads <- sample_sums(physeq.noncont)
sum(totalreads) #Total reads:1062639

#Remove blanks because we already dealt with potential contaminants
physeq.noncont <- subset_samples(physeq.noncont, is.neg != "TRUE")

#Prune singletons (these are reads that are only found once)
physeq.prune <- prune_taxa(taxa_sums(physeq.noncont) > 1, physeq.noncont)

saveRDS(physeq.prune, "~/Desktop/comparison/97_OTUS/phy97_prune.RDS")
physeq.prune <- readRDS("~/Dropbox/MacBook_Transfer/comparison/97_OTUS/phy97_prune.RDS")
data.prune <- as(sample_data(physeq.prune), "data.frame")

#Unrarefied - cut out everything with reads below 998
physeq.nr <- prune_samples(sample_sums(physeq.prune)>=998, physeq.prune)
bad.samples <- c("KI15BFMD086", "KI15BFMD107", "KI15BFMD108", "KI15BFMD119", "KI15BFMD123", "KI15BFMD135", "KI15BFMD145", "KI15BFMD195", "KI15BFMD196", "KI15BFMD208", "KI15BFMD208", "KI15BFMD230", "KI15BFMD232_EMP")
all.samples <- sample_names(physeq.nr)
all.samples <- all.samples[!(all.samples %in% bad.samples)]
physeq.nr <- prune_samples(all.samples, physeq.nr)
data.nr <- as(sample_data(physeq.nr), "data.frame")

totalreads <- sample_sums(physeq.nr)
sum(totalreads) #Total reads:953396
tax_table(physeq.nr) #2174 OTUs
View(data.nr)

#Rarefied to 1000 reads
physeq.r <- rarefy_even_depth(physeq.prune, sample.size = 998, rngseed = 711) 
bad.samples.r <- c("KI15BFMD086", "KI15BFMD107", "KI15BFMD108", "KI15BFMD119", "KI15BFMD123", "KI15BFMD135", "KI15BFMD145", "KI15BFMD195", "KI15BFMD196", "KI15BFMD208", "KI15BFMD230", "KI15BFMD232_EMP")
all.samples.r <- sample_names(physeq.r)
all.samples.r <- all.samples.r[!(all.samples.r %in% bad.samples.r)]
physeq.r <- prune_samples(all.samples.r, physeq.r)
data.r <- as(sample_data(physeq.r), "data.frame")

View(data.nr)

#Relative abundance
physeq.ra <- transform_sample_counts(physeq.nr, function(x) x/ sum(x)) 
data.ra <- as(sample_data(physeq.ra), "data.frame")

#Relative abundance with rarefied data
physeq.r.ra <- transform_sample_counts(physeq.r, function(x) x/ sum(x)) 
data.r.ra <- as(sample_data(physeq.r), "data.frame")

#First check barplots by Phylum

p <- plot_bar(physeq.ra, x = "sample_label", fill = "Phylum")
#p <- plot_bar(physeq.r.ra, x = "sample_label", fill = "Phylum")

p$data$sample_label <- factor(x = p$data$sample_label, levels = c("KI15BFMD082", "KI15BFMD089", "KI15BFMD124", "KI15BFMD142",
                                                                  "KI15BFMD144", "KI15BFMD146", "KI15BFMD190", "KI15BFMD193",
                                                                  "KI15BFMD221", "KI15BFMD228", "KI15BFMD240", "KI15BFMD079",
                                                                  "KI15BFMD091", "KI15BFMD105", "KI15BFMD112", 
                                                                  "KI15BFMD115", "KI15BFMD122", "KI15BFMD141", "KI15BFMD189",
                                                                  "KI15BFMD201", "KI15BFMD229", "KI15BFMD235"))

nb.cols <- 33
mycolors <- colorRampPalette(brewer.pal(11, "RdYlBu"))(nb.cols)
p1 <- p + geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = 'stack', width = 0.7) +
  #guides(fill = FALSE, color = FALSE) +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) +
  ylab("Relative Abundance") +
  xlab("Sample Label") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


p1 + facet_grid(rows = vars(seq_platform)) 

ggsave("~/Desktop/barplot_phy_97.pdf")
#ggsave("~/Desktop/barplot_phy_97_rare.pdf")

#Look at top 10 phyla relative abundances by platform
phabund <- tax_glom(physeq.ra, "Phylum")
miseq_phabund <- subset_samples(phabund, seq_platform=="MiSeq")
hiseq_phabund <- subset_samples(phabund, seq_platform == "HiSeq")

TopNOTUs = names(sort(taxa_sums(miseq_phabund), TRUE)[1:10])
newdat = prune_taxa(TopNOTUs, miseq_phabund)
newdatmelt1 <- psmelt(newdat)
above_zero <- newdatmelt1[which(newdatmelt1$Abundance > 0),]
print(levels(above_zero$Genus))

sum.miseq <- ddply(above_zero, c("Phylum"), summarise,
             N = length(Abundance), 
             mean = mean(Abundance),
             sd = sd(Abundance), 
             se = sd/sqrt(N)
)

TopNOTUs = names(sort(taxa_sums(hiseq_phabund), TRUE)[1:10])
newdat = prune_taxa(TopNOTUs, hiseq_phabund)
newdatmelt1 <- psmelt(newdat)
above_zero <- newdatmelt1[which(newdatmelt1$Abundance > 0),]
print(levels(above_zero$Genus))

sum.hiseq <- ddply(above_zero, c("Phylum"), summarise,
                   N = length(Abundance), 
                   mean = mean(Abundance),
                   sd = sd(Abundance), 
                   se = sd/sqrt(N)
)

##Alpha diversity
#Observed Sp. Richness, Shannon, and Faith's PD using rarefied data (Chao1 uses unrarefied, see below)
#First need to input function estimate_richness_wPD (uses the estimate_richness command in phyloseq
#but modified to add faith's phylogenetic distance metric to it as "FaithPD")
#see estimate_richness_wFPD.R file
erich <- estimate_richness_wPD(physeq.r, measures = c("Observed", "Shannon", "FaithPD"))
erich$platform <- data.r$seq_platform
erich$sample <- data.r$sample_label
erich$species <- data.r$field_host_name

#Separate by host species and re-order data for use in paired t-tests
erich.mfol <- erich %>% filter(species == "Montipora foliosa")
erich.mfol <-  erich.mfol[order(erich.mfol$platform, erich.mfol$sample),]
erich.plob <- erich %>% filter(species == "Porites lobata")
erich.plob <-  erich.plob[order(erich.plob$platform, erich.plob$sample),]

#Paired t-tests for all Alpha diversity metrics
#need the following package for qq-plots to check normality
library(ggpubr)

##M. aequituberculata: Observed Sp. Richness
#Check normality
hist(erich.mfol$Observed)
ggqqplot(erich.mfol$Observed)
hist(sqrt(erich.mfol$Observed))
ggqqplot(sqrt(erich.mfol$Observed))
shapiro.test(sqrt(erich.mfol$Observed))
#results: W = 0.98501, p-value = 0.9752

#t-test
t.mfol.sqOb <- t.test(sqrt(Observed) ~ platform, data = erich.mfol, paired = TRUE)
#results: t = -0.13484, df = 10, p-value = 0.8954

erich.mfol$sqOb <- sqrt(erich.mfol$Observed)
#use the function summarySE
t.mfol.ob.sum <- summarySE(erich.mfol, measurevar = "sqOb", groupvars = "platform")
mfol.ob <- ggplot(t.mfol.ob.sum, aes(y = sqOb, x = platform)) + 
  geom_errorbar(aes(ymin = sqOb -se, ymax = sqOb +se), width = 0.1) +
  geom_point() +
  geom_line() +
  scale_y_continuous("Sqrt(Observed Richness)") + 
  theme_classic()

##P.lobata: Observed Species Richness
hist(erich.plob$Observed) 
ggqqplot(erich.plob$Observed) 
hist(log(erich.plob$Observed))
ggqqplot(log(erich.plob$Observed))
shapiro.test(log(erich.plob$Observed))
#results: W = 0.94146, p-value = 0.2122

t.plob.ob <- t.test(log(Observed) ~ platform, data = erich.plob, paired = TRUE)
#results: t = 2.0437, df = 10, p-value = 0.06822

erich.plob$logOb <- log(erich.plob$Observed)
#use the function summarySE
t.plob.ob.sum <- summarySE(erich.plob, measurevar = "logOb", groupvars = "platform")
plob.ob <- ggplot(t.plob.ob.sum, aes(y = logOb, x = platform)) + 
  geom_errorbar(aes(ymin = logOb -se, ymax = logOb +se), width = 0.1) +
  geom_point() +
  geom_line() +
  scale_y_continuous("log(Observed Richness)") + 
  theme_classic()

#M. aquituberculata: Shannon
hist(erich.mfol$Shannon)
ggqqplot(erich.mfol$Shannon)
shapiro.test(erich.mfol$Shannon)
#results:W = 0.98834, p-value = 0.9933

#t-test
t.mfol.sh <- t.test(Shannon ~ platform, data = erich.mfol, paired = TRUE)
#results: t = 0.057506, df = 10, p-value = 0.9553

#use the function summarySE
t.mfol.sh.sum <- summarySE(erich.mfol, measurevar = "Shannon", groupvars = "platform")
mfol.sh <- ggplot(t.mfol.sh.sum, aes(y = Shannon, x = platform)) + 
  geom_errorbar(aes(ymin = Shannon -se, ymax = Shannon +se), width = 0.1) +
  geom_point() +
  geom_line() +
  scale_y_continuous("Shannon Diversity Index") + 
  theme_classic()


#P. lobata: Shannon Diversity
hist(erich.plob$Shannon)
ggqqplot(erich.plob$Shannon)
shapiro.test(erich.plob$Shannon)
#results: W = 0.95476, p-value = 0.3911

#t-test
t.plob.sh <- t.test(Shannon ~ platform, data = erich.plob, paired = TRUE)
#results:t = 2.2674, df = 10, p-value = 0.04678

#use the function summarySE
t.plob.sh.sum <- summarySE(erich.plob, measurevar = "Shannon", groupvars = "platform")
plob.sh <- ggplot(t.plob.sh.sum, aes(y = Shannon, x = platform)) + 
  geom_errorbar(aes(ymin = Shannon -se, ymax = Shannon +se), width = 0.1) +
  geom_point() +
  geom_line() +
  scale_y_continuous("Shannon Diversity Index") + 
  theme_classic()

##Faith's PD

#M. aquituberculata FPD
hist(erich.mfol$FaithPD)
ggqqplot(erich.mfol$FaithPD)
hist(log(erich.mfol$FaithPD))
ggqqplot(log(erich.mfol$FaithPD))
shapiro.test(log(erich.mfol$FaithPD))
#results: W = 0.92606, p-value = 0.1016

#t-test
t.mfol.fpd <- t.test(log(FaithPD) ~ platform, data = erich.mfol, paired = TRUE)
#results: t = -0.69056, df = 10, p-value = 0.5056

erich.mfol$logFPD <- log(erich.mfol$FaithPD)
#use the function summarySE
t.mfol.fpd.sum <- summarySE(erich.mfol, measurevar = "logFPD", groupvars = "platform")
mfol.fpd <- ggplot(t.mfol.fpd.sum, aes(y = logFPD, x = platform)) + 
  geom_errorbar(aes(ymin = logFPD -se, ymax = logFPD +se), width = 0.1) +
  geom_point() +
  geom_line() +
  scale_y_continuous("Faith's Phylogenetic Diversity") + 
  theme_classic()

#P. lobata FPD
hist(erich.plob$FaithPD)
ggqqplot(erich.plob$FaithPD)
hist(log(erich.plob$FaithPD))
ggqqplot(log(erich.plob$FaithPD))
shapiro.test(log(erich.plob$FaithPD))
#results: W = 0.92675, p-value = 0.105

#t-test
t.plob.fpd <- t.test(log(FaithPD) ~ platform, data = erich.plob, paired = TRUE)
#results: t = 2.4852, df = 10, p-value = 0.03225

erich.plob$logFPD <- log(erich.plob$FaithPD)
#use the function summarySE
t.plob.fpd.sum <- summarySE(erich.plob, measurevar = "logFPD", groupvars = "platform")
plob.fpd <- ggplot(t.plob.fpd.sum, aes(y = logFPD, x = platform)) + 
  geom_errorbar(aes(ymin = logFPD -se, ymax = logFPD +se), width = 0.1) +
  geom_point() +
  geom_line() +
  scale_y_continuous("log(Faith's Phylogenetic Diversity)") + 
  theme_classic()

#Chao1

erich.nr <- estimate_richness(physeq.nr, measures = "Chao1")
erich.nr$platform <- data.nr$seq_platform
erich.nr$sample <- data.nr$sample_label
erich.nr$species <- data.nr$field_host_name

erich.nr.mfol <- erich.nr %>% filter(species == "Montipora foliosa")
erich.nr.mfol <-  erich.nr.mfol[order(erich.nr.mfol$platform, erich.nr.mfol$sample),]
erich.nr.plob <- erich.nr %>% filter(species == "Porites lobata")
erich.nr.plob <-  erich.nr.plob[order(erich.nr.plob$platform, erich.nr.plob$sample),]

#M. aquituberculata Chao1
hist(erich.nr.mfol$Chao1)
ggqqplot(erich.nr.mfol$Chao1)
hist(log(erich.nr.mfol$Chao1))
ggqqplot(log(erich.nr.mfol$Chao1))
shapiro.test(log(erich.nr.mfol$Chao1))
#results: W = 0.95407, p-value = 0.3793

#t-test
t.mfol.ch <- t.test(log(Chao1) ~ platform, data = erich.nr.mfol, paired = TRUE)
#results:t = 1.4041, df = 10, p-value = 0.1906

erich.nr.mfol$logCh <- log(erich.nr.mfol$Chao1)
#use the function summarySE
t.mfol.ch.sum <- summarySE(erich.nr.mfol, measurevar = "logCh", groupvars = "platform")
mfol.ch <- ggplot(t.mfol.ch.sum, aes(y = logCh, x = platform)) + 
  geom_errorbar(aes(ymin = logCh -se, ymax = logCh +se), width = 0.1) +
  geom_point() +
  geom_line() +
  scale_y_continuous("log(Chao1)") + 
  theme_classic()

#P. lobata Chao1
hist(erich.nr.plob$Chao1)
ggqqplot(erich.nr.plob$Chao1)
hist(log(erich.nr.plob$Chao1))
ggqqplot(log(erich.nr.plob$Chao1))
shapiro.test(log(erich.nr.plob$Chao1))
#results: W = 0.95679, p-value = 0.4272

#t-test
t.plob.ch <- t.test(log(Chao1) ~ platform, data = erich.nr.plob, paired = TRUE)
#results:t = 1.7554, df = 10, p-value = 0.1097 

erich.nr.plob$logCh <- log(erich.nr.plob$Chao1)
#use the function summarySE
t.plob.ch.sum <- summarySE(erich.nr.plob, measurevar = "logCh", groupvars = "platform")
plob.ch <- ggplot(t.plob.ch.sum, aes(y = logCh, x = platform)) + 
  geom_errorbar(aes(ymin = logCh -se, ymax = logCh +se), width = 0.1) +
  geom_point() +
  geom_line() +
  scale_y_continuous("log(Chao1)") + 
  theme_classic()


##Betadiversity
#Split data by host species

mfol.nr <- subset_samples(physeq.nr, field_host_name == "Montipora foliosa")
mfol.data.nr <- as(sample_data(mfol.nr), "data.frame")
mfol.ra <- transform_sample_counts(mfol.nr, function(x) x/ sum(x)) 
mfol.data.ra <- as(sample_data(mfol.ra), "data.frame")

plob.nr <- subset_samples(physeq.nr, field_host_name == "Porites lobata")
plob.data.nr <- as(sample_data(plob.nr), "data.frame")
plob.ra <- transform_sample_counts(plob.nr, function(x) x/ sum(x)) 
plob.data.ra <- as(sample_data(plob.ra), "data.frame")


## Binary Jaccard
#Montipora foliosa
jc.mfol <- phyloseq::distance(mfol.nr, method = "jaccard", binary = TRUE)
adonis(jc.mfol ~ seq_platform, strata = mfol.data.nr$sample_label, data = mfol.data.nr)
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    0.5371 0.53714  1.2356 0.05818  0.001 ***
#  Residuals    20    8.6945 0.43473         0.94182           
#Total        21    9.2317                 1.00000 
permutest(betadisper(jc.mfol, mfol.data.nr$seq_platform, type = "centroid"))
#          Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.0007436 0.00074358 1.5201    999  0.238
#Residuals 20 0.0097833 0.00048917  


#Porites lobata
jc.plob <- phyloseq::distance(plob.nr, method = "jaccard", binary = TRUE)
adonis(jc.plob ~ seq_platform, strata = plob.data.nr$sample_label, data = plob.data.nr)
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    0.6268 0.62677   1.487 0.0692  0.001 ***
#Residuals    20    8.4300 0.42150         0.9308           
#Total        21    9.0568                 1.0000  
permutest(betadisper(jc.plob, plob.data.nr$seq_platform, type = "centroid"))
#Df    Sum Sq    Mean Sq     F N.Perm Pr(>F)
#Groups     1 0.000141 0.00014095 0.0801    999  0.778
#Residuals 20 0.035203 0.00176016   

##Bray curtis - uses relative abundance to account for differences in sequencing depth
#Montipora foliosa
bc.mfol <- phyloseq::distance(mfol.ra, method = "bray")
adonis(bc.mfol ~ seq_platform, strata = mfol.data.ra$sample_label, data = mfol.data.ra)
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.4860 0.48596  1.1781 0.05563  0.002 **
#  Residuals    20    8.2497 0.41248         0.94437          
#Total        21    8.7356                 1.00000 
permutest(betadisper(bc.mfol, mfol.data.ra$seq_platform, type = "centroid"))
#          Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.0000306 0.00003057 0.0205    999  0.884
#Residuals 20 0.0298830 0.00149415

#Porites lobata
bc.plob <- phyloseq::distance(plob.ra, method = "bray")
adonis(bc.plob ~ seq_platform, strata = plob.data.ra$sample_label, data = plob.data.ra)
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    0.6720 0.67201  2.0676 0.09369  0.001 ***
#Residuals    20    6.5005 0.32502         0.90631           
#Total        21    7.1725                 1.00000 
permutest(betadisper(bc.plob, plob.data.ra$seq_platform, type = "centroid"))
#Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
#Groups     1 0.16747 0.167471 6.9726    999  0.008 **
#Residuals 20 0.48037 0.024018 


##Weighted unifrac
#M.fol
w.mfol <- phyloseq::distance(mfol.ra, method = "wunifrac")
adonis(w.mfol ~ seq_platform, strata = mfol.data.ra$sample_label, data = mfol.data.ra)
#         Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#seq_platform  1   0.03426 0.03426  1.4281 0.06664  0.098 .
#Residuals    20   0.47981 0.02399         0.93336         
#Total        21   0.51407                 1.00000  
permutest(betadisper(w.mfol, mfol.data.ra$seq_platform, type = "centroid"))
#        Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)   
#Groups     1 0.023109 0.0231085 8.4243    999  0.002 **
#  Residuals 20 0.054862 0.0027431   

#Porites lobata
w.plob <- phyloseq::distance(plob.ra, method = "wunifrac")
adonis(w.plob ~ seq_platform, strata = plob.data.ra$sample_label, data = plob.data.ra)
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1   0.01941 0.019413  1.2793 0.06012  0.003 **
#Residuals    20   0.30351 0.015175         0.93988          
#Total        21   0.32292                  1.00000 
permutest(betadisper(w.plob, plob.data.ra$seq_platform, type = "centroid"))
#    Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.002145 0.0021450 0.3985    999  0.546
#Residuals 20 0.107655 0.0053827

##Unweighted unifrac

un.mfol <- phyloseq::distance(mfol.nr, method = "uunifrac")
adonis(un.mfol ~ seq_platform, strata = mfol.data.nr$sample_label, data = mfol.data.nr)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#seq_platform  1    0.3819  0.3819  1.0946 0.05189  0.019 *
#  Residuals    20    6.9780  0.3489         0.94811         
#Total        21    7.3599                 1.00000  
permutest(betadisper(un.mfol, mfol.data.nr$seq_platform, type = "centroid"))
#        Df    Sum Sq    Mean Sq     F N.Perm Pr(>F)
#Groups     1 0.0005619 0.00056187 0.363    999  0.547
#Residuals 20 0.0309537 0.00154768  

#Porites lobata
un.plob <- phyloseq::distance(plob.nr, method = "uunifrac")
adonis(un.plob ~ seq_platform, strata = plob.data.nr$sample_label, data = plob.data.nr)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.5342 0.53423  1.3923 0.06509  0.005 **
#Residuals    20    7.6738 0.38369         0.93491          
#Total        21    8.2081                 1.00000    
permutest(betadisper(un.plob, plob.data.nr$seq_platform, type = "centroid"))
#Df   Sum Sq   Mean Sq     F N.Perm Pr(>F)   
#Groups     1 0.015205 0.0152050 17.77    999  0.001 ***
#Residuals 20 0.017114 0.0008557         

##Ordinations 
#Jaccard
#Montipora foliosa
ord.mfol.jc <- ordinate(mfol.nr, "NMDS", "jaccard", binary = TRUE, trymax = 100) 
stressplot(ord.mfol.jc)
scores.mfol.jc <- as.data.frame(scores(ord.mfol.jc))
scores.mfol.jc$platform <- mfol.data.nr$seq_platform
scores.mfol.jc$sample <- mfol.data.nr$sample_label
ggplot(scores.mfol.jc, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = platform, shape = platform), size = 4) +
  #geom_text(label = rownames(scores), nudge_x = 0.05, nudge_y = 0.05,check_overlap = T) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = platform), linetype = 2) +
  ggtitle("M. aequituberculata Jaccard (OTUs)") +
  theme_classic()
ggsave("~/Desktop/jc_mfol_nmds_97.pdf")

#Porites lobata
ord.plob.jc <- ordinate(plob.nr, "NMDS", "jaccard", binary = TRUE, trymax = 100) 
stressplot(ord.plob.jc)
scores.plob.jc <- as.data.frame(scores(ord.plob.jc))
scores.plob.jc$platform <- plob.data.nr$seq_platform
scores.plob.jc$sample <- plob.data.nr$sample_label
ggplot(scores.plob.jc, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = platform, shape = platform), size = 4) +
  #geom_text(label = rownames(scores), nudge_x = 0.05, nudge_y = 0.05,check_overlap = T) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = platform), linetype = 2) +
  ggtitle("P. lobata Jaccard (OTUs)") +
  theme_classic()
ggsave("~/Desktop/jc_plob_nmds_97.pdf")


#Bray Curtis
#Montipora foliosa
ord.mfol <- ordinate(mfol.ra, "NMDS", "bray", trymax = 100) 
stressplot(ord.mfol)
scores.mfol <- as.data.frame(scores(ord.mfol))
scores.mfol$platform <- mfol.data.ra$seq_platform
scores.mfol$sample <- mfol.data.ra$sample_label
ggplot(scores.mfol, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = platform, shape = platform), size = 4) +
  #geom_text(label = rownames(scores), nudge_x = 0.05, nudge_y = 0.05,check_overlap = T) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = platform), linetype = 2) +
  ggtitle("M. aequituberculata Bray Curtis (OTUs)") +
  theme_classic()
ggsave("~/Desktop/br_mfol_nmds_97.pdf")

#Porites lobata
ord.plob <- ordinate(plob.ra, "NMDS", "bray", trymax = 500) 
stressplot(ord.plob)
scores.plob <- as.data.frame(scores(ord.plob))
scores.plob$platform <- plob.data.ra$seq_platform
scores.plob$sample <- plob.data.ra$sample_label
ggplot(scores.plob, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = platform, shape = platform), size = 4) +
  #geom_text(label = rownames(scores), nudge_x = 0.05, nudge_y = 0.05,check_overlap = T) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = platform), linetype = 2) +
  ggtitle("P. lobata Bray Curtis (OTUs)") +
  theme_classic()
ggsave("~/Desktop/br_plob_nmds_97.pdf")


#Unweighted Unifrac

ord.un.mfol <- ordinate(mfol.ra, "NMDS", "uunifrac") #trymax is not accepted here
stressplot(ord.un.mfol)
scores.un.mfol <- as.data.frame(scores(ord.un.mfol))
scores.un.mfol$platform <- mfol.data.ra$seq_platform
scores.un.mfol$sample <- mfol.data.ra$sample_label
ggplot(scores.un.mfol, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = platform, shape = platform), size = 4) +
  #geom_text(label = rownames(scores), nudge_x = 0.05, nudge_y = 0.05,check_overlap = T) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = platform), linetype = 2) +
  ggtitle("M. aequituberculata Unweighted UniFrac (OTUs)") +
  theme_classic()
ggsave("~/Desktop/un_mfol_nmds_97.pdf")

#Porites lobata
ord.un.plob <- ordinate(plob.ra, "NMDS", "uunifrac") 
stressplot(ord.un.plob)
scores.un.plob <- as.data.frame(scores(ord.un.plob))
scores.un.plob$platform <- plob.data.ra$seq_platform
scores.un.plob$sample <- plob.data.ra$sample_label
ggplot(scores.un.plob, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = platform, shape = platform), size = 4) +
  #geom_text(label = rownames(scores), nudge_x = 0.05, nudge_y = 0.05,check_overlap = T) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = platform), linetype = 2) +
  ggtitle("P. lobata Unweighted UniFrac (OTUs)") +
  theme_classic()
ggsave("~/Desktop/un_plob_nmds_97.pdf")


##Weighted Unifrac
#Montipora foliosa
ord.wun.mfol <- ordinate(mfol.ra, "NMDS", "wunifrac") #trymax is not accepted here
stressplot(ord.wun.mfol)
scores.wun.mfol <- as.data.frame(scores(ord.wun.mfol))
scores.wun.mfol$platform <- mfol.data.ra$seq_platform
scores.wun.mfol$sample <- mfol.data.ra$sample_label
ggplot(scores.wun.mfol, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = platform, shape = platform), size = 4) +
  #geom_text(label = rownames(scores), nudge_x = 0.05, nudge_y = 0.05,check_overlap = T) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = platform), linetype = 2) +
  ggtitle("M. aequituberculata Weighted UniFrac (OTUs)") +
  theme_classic()
ggsave("~/Desktop/wun_mfol_nmds_97.pdf")

#Porites lobata
ord.wun.plob <- ordinate(plob.ra, "NMDS", "wunifrac") 
stressplot(ord.wun.plob)
scores.wun.plob <- as.data.frame(scores(ord.wun.plob))
scores.wun.plob$platform <- plob.data.ra$seq_platform
scores.wun.plob$sample <- plob.data.ra$sample_label
ggplot(scores.wun.plob, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = platform, shape = platform), size = 4) +
  #geom_text(label = rownames(scores.wun.plob)) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = platform), linetype = 2) +
  ggtitle("P. lobata Weighted UniFrac (OTUs)") +
  theme_classic()
ggsave("~/Desktop/wun_plob_nmds_97.pdf")



##Betadiversity by tax_glom - Family & Phylum

mfol.fam.ra <- mfol.ra %>% 
  tax_glom(taxrank = "Family") 
mfol.fam.ra.data <- as(sample_data(mfol.fam.ra), "data.frame")

bc.mfol.fam <- phyloseq::distance(mfol.fam.ra, method = "bray")
adonis(bc.mfol.fam ~ seq_platform, strata = mfol.fam.ra.data$sample_label, data = mfol.fam.ra.data)
##   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#seq_platform  1    0.3317 0.33172  1.1824 0.05582  0.063 .
#Residuals    20    5.6109 0.28054         0.94418         
#Total        21    5.9426                 1.00000  
permutest(betadisper(bc.mfol.fam, mfol.fam.ra.data$seq_platform, type = "centroid"))
#   Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.000197 0.0001972 0.0203    999  0.892
#Residuals 20 0.194234 0.0097117 

w.mfol.fam <- phyloseq::distance(mfol.fam.ra, method = "wunifrac")
adonis(w.mfol.fam ~ seq_platform, strata = mfol.fam.ra.data$sample_label, data = mfol.fam.ra.data)
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#seq_platform  1   0.11114 0.11114   1.067 0.05065  0.125
#Residuals    20   2.08318 0.10416         0.94935       
#Total        21   2.19432                 1.00000    

permutest(betadisper(w.mfol.fam, mfol.fam.ra.data$seq_platform, type = "centroid"))
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.000024 0.0000238 0.0037    999  0.954
#Residuals 20 0.127063 0.0063531   

mfol.fam.nr <- mfol.nr %>% 
  tax_glom(taxrank = "Family") 
mfol.fam.nr.data <- as(sample_data(mfol.fam.nr), "data.frame")

jc.mfol.fam <- phyloseq::distance(mfol.fam.nr, method = "jaccard", binary = TRUE)
adonis(jc.mfol.fam ~ seq_platform, strata = mfol.fam.nr.data$sample_label, data = mfol.fam.nr.data)
#        Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.4502 0.45021  1.5854 0.07345  0.004 **
#  Residuals    20    5.6796 0.28398         0.92655          
#Total        21    6.1298                 1.00000    
permutest(betadisper(jc.mfol.fam, mfol.fam.nr.data$seq_platform, type = "centroid"))
#    Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.007014 0.0070145 1.8274    999  0.206
#Residuals 20 0.076771 0.0038386 
un.mfol.fam <- phyloseq::distance(mfol.fam.nr, method = "uunifrac")
adonis(un.mfol.fam ~ seq_platform, strata = mfol.fam.nr.data$sample_label, data = mfol.fam.nr.data)
#       Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.3752 0.37517  1.9126 0.08728  0.002 **
#  Residuals    20    3.9232 0.19616         0.91272          
#Total        21    4.2983                 1.00000       
permutest(betadisper(un.mfol.fam, mfol.fam.nr.data$seq_platform, type = "centroid"))
#    Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.004146 0.0041462 0.5435    999  0.525
#Residuals 20 0.152570 0.0076285  

mfol.phy.ra <- mfol.ra %>%
  tax_glom(taxrank = "Phylum")
mfol.phy.ra.data <- as(sample_data(mfol.phy.ra), "data.frame")
bc.mfol.phyl <- phyloseq::distance(mfol.phy.ra, method = "bray")
adonis(bc.mfol.phyl ~ seq_platform, strata =mfol.phy.ra.data$sample_label, data = mfol.phy.ra.data)
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
#seq_platform  1   0.07661 0.076614 0.78946 0.03797  0.384
#Residuals    20   1.94093 0.097047         0.96203       
#Total        21   2.01755                  1.00000   
permutest(betadisper(bc.mfol.phyl, mfol.phy.ra.data$seq_platform, type = "centroid"))
#      Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.035867 0.035867 2.8032    999  0.116
#Residuals 20 0.255895 0.012795 

w.mfol.phy <- phyloseq::distance(mfol.phy.ra, method = "wunifrac")
adonis(w.mfol.phy ~ seq_platform, strata = mfol.phy.ra.data$sample_label, data = mfol.phy.ra.data)
#    Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
#seq_platform  1   0.07632 0.076316  1.0421 0.04953  0.313
#Residuals    20   1.46459 0.073230         0.95047       
#Total        21   1.54091                  1.00000  
permutest(betadisper(w.mfol.phy, mfol.phy.ra.data$seq_platform, type = "centroid"))
#        Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.015543 0.015543 1.4631    999  0.265
#Residuals 20 0.212472 0.010624  

mfol.phy.nr <- mfol.nr %>% 
  tax_glom(taxrank = "Phylum") 
mfol.phy.nr.data <- as(sample_data(mfol.phy.nr), "data.frame")

jc.mfol.phy <- phyloseq::distance(mfol.phy.nr, method = "jaccard", binary = TRUE)
adonis(jc.mfol.phy ~ seq_platform, strata = mfol.phy.nr.data$sample_label, data = mfol.phy.nr.data)
#   Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.2979 0.297900  3.3332 0.14285  0.004 **
#  Residuals    20    1.7875 0.089373         0.85715          
#Total        21    2.0854                  1.00000        
permutest(betadisper(jc.mfol.phy, mfol.phy.nr.data$seq_platform, type = "centroid"))
#  Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.010758 0.010759 1.1138    999  0.332
#Residuals 20 0.193179 0.009659  

un.mfol.phy <- phyloseq::distance(mfol.phy.nr, method = "uunifrac")
adonis(un.mfol.phy ~ seq_platform, strata = mfol.phy.nr.data$sample_label, data = mfol.phy.nr.data)
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#seq_platform  1   0.33577 0.33577  3.8511 0.16146  0.012 *
#  Residuals    20   1.74375 0.08719         0.83854         
#Total        21   2.07951                 1.00000           
permutest(betadisper(un.mfol.phy, mfol.phy.nr.data$seq_platform, type = "centroid"))
#    Df   Sum Sq   Mean Sq     F N.Perm Pr(>F)
#Groups     1 0.013808 0.0138081 1.558    999  0.227
#Residuals 20 0.177250 0.0088625 

##Porites
plob.fam.ra <- plob.ra %>% 
  tax_glom(taxrank = "Family") 
plob.fam.ra.data <- as(sample_data(plob.fam.ra), "data.frame")
bc.plob.fam <- phyloseq::distance(plob.fam.ra, method = "bray")
adonis(bc.plob.fam ~ seq_platform, strata = plob.fam.ra.data$sample_label, data = plob.fam.ra.data)

#         Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.6822 0.68224  2.6902 0.11856  0.001 ***
#Residuals    20    5.0719 0.25360         0.88144           
#Total        21    5.7542                 1.00000        
permutest(betadisper(bc.plob.fam, plob.fam.ra.data$seq_platform, type = "centroid"))
#Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
#Groups     1 0.12565 0.125653 4.7188    999  0.041 *
# Residuals 20 0.53257 0.026628               

w.plob.fam <- phyloseq::distance(plob.fam.ra, method = "wunifrac")
adonis(w.plob.fam ~ seq_platform, strata = plob.fam.ra.data$sample_label, data = plob.fam.ra.data)
#          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1   0.1504 0.15040  2.0418 0.09263  0.006 **
#Residuals    20    1.4732 0.07366         0.90737          
#Total        21    1.6236                 1.00000 
permutest(betadisper(w.plob.fam, plob.fam.ra.data$seq_platform, type = "centroid"))
#        Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
#Groups     1 0.026471 0.0264708 2.7464    999  0.107
#Residuals 20 0.192767 0.0096384        

plob.fam.nr <- plob.nr %>% 
  tax_glom(taxrank = "Family") 
plob.fam.nr.data <- as(sample_data(plob.fam.nr), "data.frame")

jc.plob.fam <- phyloseq::distance(plob.fam.nr, method = "jaccard", binary = TRUE)
adonis(jc.plob.fam ~ seq_platform, strata = plob.fam.nr.data$sample_label, data = plob.fam.nr.data)
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.6597 0.65971  2.0233 0.09187  0.001 ***
#Residuals    20    6.5210 0.32605         0.90813           
#Total        21    7.1807                 1.00000       
permutest(betadisper(jc.plob.fam, plob.fam.nr.data$seq_platform, type = "centroid"))
#     Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)   
#Groups     1 0.031993 0.031993 20.976    999  0.001 ***
#Residuals 20 0.030505 0.001525    

un.plob.fam <- phyloseq::distance(plob.fam.nr, method = "uunifrac")
adonis(un.plob.fam ~ seq_platform, strata = plob.fam.nr.data$sample_label, data = plob.fam.nr.data)
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.5410 0.54098   2.097 0.0949  0.002 **
#Residuals    20    5.1596 0.25798         0.9051          
#Total        21    5.7006                 1.0000   
permutest(betadisper(un.plob.fam, plob.fam.nr.data$seq_platform, type = "centroid"))
#      Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)    
#Groups     1 0.054975 0.054975 21.658    999  0.001 ***
#Residuals 20 0.050767 0.002538  

plob.phy.ra <- plob.ra %>%
  tax_glom(taxrank = "Phylum")
plob.phy.ra.data <- as(sample_data(plob.phy.ra), "data.frame")
bc.plob.phyl <- phyloseq::distance(plob.phy.ra, method = "bray")
adonis(bc.plob.phyl ~ seq_platform, strata = plob.phy.ra.data$sample_label, data = plob.phy.ra.data)
#         Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
#seq_platform  1   0.07297 0.072966  1.2314 0.058   0.02 *
#Residuals    20   1.18509 0.059255         0.942         
#Total        21   1.25806                  1.000       
permutest(betadisper(bc.plob.phyl, plob.phy.ra.data$seq_platform, type = "centroid"))
#      Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.00349 0.0034865 0.1427    999  0.605
#Residuals 20 0.48870 0.0244348    

w.plob.phy <- phyloseq::distance(plob.phy.ra, method = "wunifrac")
adonis(w.plob.phy ~ seq_platform, strata = plob.phy.ra.data$sample_label, data = plob.phy.ra.data)
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1   0.09016 0.090165  1.5592 0.07232  0.006 **
#Residuals    20   1.15656 0.057828         0.92768          
#Total        21   1.24672                  1.00000  
permutest(betadisper(w.plob.phy, plob.phy.ra.data$seq_platform, type = "centroid"))
#      Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.00149 0.0014851 0.0777    999  0.758
#Residuals 20 0.38249 0.0191244 

plob.phy.nr <- plob.nr %>% 
  tax_glom(taxrank = "Phylum") 
plob.phy.nr.data <- as(sample_data(plob.phy.nr), "data.frame")

jc.plob.phy <- phyloseq::distance(plob.phy.nr, method = "jaccard", binary = TRUE)
adonis(jc.plob.phy ~ seq_platform, strata = plob.phy.nr.data$sample_label, data = plob.phy.nr.data)
#          Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)   
#seq_platform  1    0.5623 0.56229  3.3254 0.14257  0.002 **
#Residuals    20    3.3818 0.16909         0.85743          
#Total        21    3.9441                 1.00000  
permutest(betadisper(jc.plob.phy, plob.phy.nr.data$seq_platform, type = "centroid"))
#       Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
#Groups     1 0.20383 0.203826 34.188    999  0.001 ***
#Residuals 20 0.11924 0.005962    

un.plob.phy <- phyloseq::distance(plob.phy.nr, method = "uunifrac")
adonis(un.plob.phy ~ seq_platform, strata = plob.phy.nr.data$sample_label, data = plob.phy.nr.data)
#       Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.5253 0.52533  2.7908 0.12245  0.001 ***
#Residuals    20    3.7648 0.18824         0.87755           
#Total        21    4.2901                 1.00000  
permutest(betadisper(un.plob.phy, plob.phy.nr.data$seq_platform, type = "centroid"))
#        Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
#Groups     1 0.18881 0.188810 26.39    999  0.001 ***
#Residuals 20 0.14309 0.007155 


#Filtered for truncation
##Keep only taxa that have > 0.5% relative abundance in each sample
physeq.raF = filter_taxa(physeq.ra, function(x) mean(x) < .005,TRUE)
rmtaxa = taxa_names(physeq.raF)
alltaxa = taxa_names(physeq.ra)

myTaxa = alltaxa[!alltaxa %in% rmtaxa]

physeq.nrF <- prune_taxa(myTaxa,physeq.nr)
tax_table(physeq.raF) #2144 taxa

mfol.nrF <- subset_samples(physeq.nrF, field_host_name == "Montipora foliosa")
mfol.data.nrF <- as(sample_data(mfol.nrF), "data.frame")
plob.nrF <- subset_samples(physeq.nrF, field_host_name == "Porites lobata")
plob.data.nrF <- as(sample_data(plob.nrF), "data.frame")

mfol.raF <- transform_sample_counts(mfol.nrF, function(x) x/ sum(x))
mfol.data.raF <- as(sample_data(mfol.raF), "data.frame")
plob.raF <- transform_sample_counts(plob.nrF, function(x) x/ sum(x)) 
plob.data.raF <- as(sample_data(plob.raF), "data.frame")

jc.mfol.F <- phyloseq::distance(mfol.nrF, method = "jaccard", binary = TRUE)
adonis(jc.mfol.F ~ seq_platform, strata = mfol.data.nrF$sample_label, data = mfol.data.nrF)
#       Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.8752 0.87518  3.9731 0.16573  0.003 **
#  Residuals    20    4.4055 0.22028         0.83427          
#Total        21    5.2807                 1.00000  
permutest(betadisper(jc.mfol.F, mfol.data.nrF$seq_platform, type = "centroid"))
#Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.014855 0.014856 1.2217    999  0.268
#Residuals 20 0.243200 0.012160    

un.mfol.F <- phyloseq::distance(mfol.nrF, method = "uunifrac")
adonis(un.mfol.F ~ seq_platform, strata = mfol.data.nrF$sample_label, data = mfol.data.nrF)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.5063 0.50630  3.1342 0.13548  0.002 **
#  Residuals    20    3.2308 0.16154         0.86452          
#Total        21    3.7371                 1.00000   
permutest(betadisper(un.mfol.F, mfol.data.nrF$seq_platform, type = "centroid"))
#        Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.024012 0.024012 1.7255    999  0.205
#Residuals 20 0.278312 0.013916                     

w.mfol.F <- phyloseq::distance(mfol.raF, method = "wunifrac")
adonis(w.mfol.F ~ seq_platform, strata = mfol.data.raF$sample_label, data = mfol.data.raF)
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#seq_platform  1    0.2751 0.27512  1.1236 0.05319  0.044 *
#  Residuals    20    4.8972 0.24486         0.94681         
#Total        21    5.1724                 1.00000      
permutest(betadisper(w.mfol.F, mfol.data.raF$seq_platform, type = "centroid"))
#  Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.010513 0.010513 0.8625    999  0.375
#Residuals 20 0.243782 0.012189  

jc.plob.F <- phyloseq::distance(plob.nrF, method = "jaccard", binary = TRUE)
adonis(jc.plob.F ~ seq_platform, strata = plob.data.nrF$sample_label, data = plob.data.nrF)
#        Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    0.9465 0.94653  4.3903 0.17291  0.001 ***
#  Residuals    21    4.5276 0.21560         0.82709           
#Total        22    5.4741                 1.00000 
permutest(betadisper(jc.plob.F, plob.data.nrF$seq_platform, type = "centroid"))
#   Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.003206 0.0032056 0.4283    999  0.495
#Residuals 21 0.157184 0.0074849 

u.plob.F <- phyloseq::distance(plob.nrF, method = "uunifrac")
adonis(u.plob.F ~ seq_platform, strata = plob.data.nrF$sample_label, data = plob.data.nrF)
#   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    0.8320 0.83199  4.1793 0.16598  0.001 ***
#  Residuals    21    4.1805 0.19907         0.83402           
#Total        22    5.0125                 1.00000  
permutest(betadisper(u.plob.F, plob.data.nrF$seq_platform, type = "centroid"))
#Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
#Groups     1 0.13537 0.135365 26.591    999  0.001 ***
#  Residuals 21 0.10690 0.005091  

w.plob.F <- phyloseq::distance(plob.raF, method = "wunifrac")
adonis(w.plob.F ~ seq_platform, strata = plob.data.raF$sample_label, data = plob.data.raF)
#     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    0.4374 0.43738  2.3088 0.09905  0.001 ***
#  Residuals    21    3.9783 0.18944         0.90095           
#Total        22    4.4157                 1.00000    
permutest(betadisper(w.plob.F, plob.data.raF$seq_platform, type = "centroid"))
#     Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
#Groups     1 0.13329 0.133290 5.0506    999  0.023 *
#  Residuals 21 0.55421 0.026391  

bc.F <- phyloseq::distance(mfol.raF, method = "bray")
adonis(bc.F ~ seq_platform, strata = mfol.data.raF$sample_label, data = mfol.data.raF)
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.4510 0.45103  1.3886 0.06492   0.01 **
#  Residuals    20    6.4962 0.32481         0.93508          
#Total        21    6.9472                 1.00000  
permutest(betadisper(bc.F, mfol.data.raF$seq_platform, type = "centroid"))
#          Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.005343 0.0053432 0.8905    999  0.364
#Residuals 20 0.120012 0.0060006  

bc.plob.F <- phyloseq::distance(plob.raF, method = "bray")
adonis(bc.plob.F ~ seq_platform, strata = plob.data.raF$sample_label, data = plob.data.raF)
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.7513 0.75125  2.6761 0.11303  0.003 **
#  Residuals    21    5.8952 0.28072         0.88697          
#Total        22    6.6464                 1.00000           
permutest(betadisper(bc.plob.F, plob.data.raF$seq_platform, type = "centroid"))
#          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)   
#Groups     1 0.22984 0.229840 8.7323    999  0.008 **
#  Residuals 21 0.55273 0.026321 


###Filtered for truncation###
#Keep only taxa that have > 1% relative abundance in each sample
physeq.raF = filter_taxa(physeq.ra, function(x) mean(x) < .01,TRUE)
rmtaxa = taxa_names(physeq.raF)
alltaxa = taxa_names(physeq.ra)

myTaxa = alltaxa[!alltaxa %in% rmtaxa]

physeq.nrF <- prune_taxa(myTaxa,physeq.nr)
tax_table(physeq.nrF) #11 taxa

mfol.nrF <- subset_samples(physeq.nrF, field_host_name == "Montipora foliosa")
mfol.data.nrF <- as(sample_data(mfol.nrF), "data.frame")
plob.nrF <- subset_samples(physeq.nrF, field_host_name == "Porites lobata")
plob.data.nrF <- as(sample_data(plob.nrF), "data.frame")

mfol.raF <- transform_sample_counts(mfol.nrF, function(x) x/ sum(x))
mfol.data.raF <- as(sample_data(mfol.raF), "data.frame")
plob.raF <- transform_sample_counts(plob.nrF, function(x) x/ sum(x)) 
plob.data.raF <- as(sample_data(plob.raF), "data.frame")

jc.mfol.F <- phyloseq::distance(mfol.nrF, method = "jaccard", binary = TRUE)
adonis(jc.mfol.F ~ seq_platform, strata = mfol.data.nrF$sample_label, data = mfol.data.nrF)
#   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#seq_platform  1    0.4486 0.44862  2.6166 0.11569  0.007 **
#Residuals    20    3.4290 0.17145         0.88431          
#Total        21    3.8777                 1.00000   
permutest(betadisper(jc.mfol.F, mfol.data.nrF$seq_platform, type = "centroid"))
#    Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.05105 0.051047 2.7522    999   0.13
#Residuals 20 0.37096 0.018548    
un.mfol.F <- phyloseq::distance(mfol.nrF, method = "uunifrac")
adonis(un.mfol.F ~ seq_platform, strata = mfol.data.nrF$sample_label, data = mfol.data.nrF)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#Seq_platform  1    0.42987 0.42987  3.6812 0.15545  0.005 **
#Residuals    20   2.33550 0.11677         0.84455          
#Total        21   2.76537                 1.00000   
permutest(betadisper(un.mfol.F, mfol.data.nrF$seq_platform, type = "centroid"))
#   Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
#Groups     1 0.06957 0.069570 3.2618    999  0.078 .
#Residuals 20 0.42657 0.021329    
w.mfol.F <- phyloseq::distance(mfol.raF, method = "wunifrac")
adonis(w.mfol.F ~ seq_platform, strata = mfol.data.raF$sample_label, data = mfol.data.raF)
#        Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#seq_platform  1    0.3155 0.31551  1.5058 0.07002  0.034 *
#Residuals    20    4.1907 0.20954         0.92998         
#Total        21    4.5063                 1.00000  
permutest(betadisper(w.mfol.F, mfol.data.raF$seq_platform, type = "centroid"))
#    Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.00212 0.0021195 0.222    999  0.621
#Residuals 20 0.19095 0.0095473 
jc.plob.F <- phyloseq::distance(plob.nrF, method = "jaccard", binary = TRUE)
adonis(jc.plob.F ~ seq_platform, strata = plob.data.nrF$sample_label, data = plob.data.nrF)
#       Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    .5877 0.58769  3.7305 0.1572  0.002 **
#Residuals    20    3.1507 0.15754         0.8428          
#Total        21    3.7384                 1.0000   
permutest(betadisper(jc.plob.F, plob.data.nrF$seq_platform, type = "centroid"))
#       Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.013432 0.013432 1.2632    999  0.265
#Residuals 20 0.212675 0.010634     
un.plob.F <- phyloseq::distance(plob.nrF, method = "uunifrac")
adonis(un.plob.F ~ seq_platform, strata = plob.data.nrF$sample_label, data = plob.data.nrF)
#      Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    0.7771 0.77708  5.5537 0.21734  0.002 **
#Residuals    20    2.7984 0.13992         0.78266          
#Total        21    3.5755                 1.00000
permutest(betadisper(un.plob.F, plob.data.nrF$seq_platform, type = "centroid"))
#       Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
#Groups     1 .17530 0.175301 24.464    999  0.001 ***
#Residuals 20 0.14332 0.007166     
w.plob.F <- phyloseq::distance(plob.raF, method = "wunifrac")
adonis(w.plob.F ~ seq_platform, strata = plob.data.raF$sample_label, data = plob.data.raF)
#         Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    0.3424 0.34245   1.988 0.09041  0.001 ***
#Residuals    20    3.4451 0.17225         0.90959           
#Total        21    3.7875                 1.00000   
permutest(betadisper(w.plob.F, plob.data.raF$seq_platform, type = "centroid"))
#      Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
#Groups     1 00.07883 0.078834 3.7599    999  0.069 .
#Residuals 20 0.41934 0.020967   
bc.F <- phyloseq::distance(mfol.raF, method = "bray")
adonis(bc.F ~ seq_platform, strata = mfol.data.raF$sample_label, data = mfol.data.raF)
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#seq_platform  1    0.3993 0.39932  1.4346 0.06693  0.012 *
#Residuals    20    5.5672 0.27836         0.93307         
#Total        21    5.9665                 1.00000        
permutest(betadisper(bc.F, mfol.data.raF$seq_platform, type = "centroid"))
#          Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.00413 0.0041303 0.6163    999  0.468
#Residuals 20 0.13403 0.0067017  
bc.plob.F <- phyloseq::distance(plob.raF, method = "bray")
adonis(bc.plob.F ~ seq_platform, strata = plob.data.raF$sample_label, data = plob.data.raF)
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.5185 0.51846  2.1098 0.09542  0.008 **
#Residuals    20    4.9147 0.24574         0.90458          
#Total        21    5.4332                 1.00000  
permutest(betadisper(bc.plob.F, plob.data.raF$seq_platform, type = "centroid"))
#          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
#Groups     1 0.13847 0.138473 5.4174    999  0.035 *
#Residuals 20 0.51122 0.025561

#Keep only taxa that have > 0.5% relative abundance in each sample
physeq.raF = filter_taxa(physeq.ra, function(x) mean(x) < .005,TRUE)
rmtaxa = taxa_names(physeq.raF)
alltaxa = taxa_names(physeq.ra)

myTaxa = alltaxa[!alltaxa %in% rmtaxa]

physeq.nrF <- prune_taxa(myTaxa,physeq.nr)
tax_table(physeq.nrF) #30 taxa

mfol.nrF <- subset_samples(physeq.nrF, field_host_name == "Montipora foliosa")
mfol.data.nrF <- as(sample_data(mfol.nrF), "data.frame")
plob.nrF <- subset_samples(physeq.nrF, field_host_name == "Porites lobata")
plob.data.nrF <- as(sample_data(plob.nrF), "data.frame")

mfol.raF <- transform_sample_counts(mfol.nrF, function(x) x/ sum(x))
mfol.data.raF <- as(sample_data(mfol.raF), "data.frame")
plob.raF <- transform_sample_counts(plob.nrF, function(x) x/ sum(x)) 
plob.data.raF <- as(sample_data(plob.raF), "data.frame")

jc.mfol.F <- phyloseq::distance(mfol.nrF, method = "jaccard", binary = TRUE)
adonis(jc.mfol.F ~ seq_platform, strata = mfol.data.nrF$sample_label, data = mfol.data.nrF)
#   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#seq_platform  1    0.8752 0.87518  3.9731 0.16573  0.002 **
#Residuals    20    4.4055 0.22028         0.83427          
#Total        21    5.2807                 1.00000  
permutest(betadisper(jc.mfol.F, mfol.data.nrF$seq_platform, type = "centroid"))
#    Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.014855 0.014856 1.2217    999  0.289
#Residuals 20 0.243200 0.012160    
un.mfol.F <- phyloseq::distance(mfol.nrF, method = "uunifrac")
adonis(un.mfol.F ~ seq_platform, strata = mfol.data.nrF$sample_label, data = mfol.data.nrF)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#Seq_platform  1    0.5063 0.50630  3.1342 0.13548  0.005 **
#Residuals    20    3.2308 0.16154         0.86452          
#Total        21    3.7371                 1.00000   
permutest(betadisper(un.mfol.F, mfol.data.nrF$seq_platform, type = "centroid"))
#   Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
#Groups     1 0.024012 0.024012 1.7255    999  0.193
#Residuals 20 0.278312 0.013916      
w.mfol.F <- phyloseq::distance(mfol.raF, method = "wunifrac")
adonis(w.mfol.F ~ seq_platform, strata = mfol.data.raF$sample_label, data = mfol.data.raF)
#        Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#seq_platform  1    0.2751 0.27512  1.1236 0.05319  0.054 .
#Residuals    20    4.8972 0.24486         0.94681         
#Total        21    5.1724                 1.00000 
permutest(betadisper(w.mfol.F, mfol.data.raF$seq_platform, type = "centroid"))
#    Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.010513 0.010513 0.8625    999  0.366
#Residuals 20 0.243782 0.012189  
jc.plob.F <- phyloseq::distance(plob.nrF, method = "jaccard", binary = TRUE)
adonis(jc.plob.F ~ seq_platform, strata = plob.data.nrF$sample_label, data = plob.data.nrF)
#       Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    1.0333 1.03330  4.8608 0.19552  0.002 **
#Residuals    20    4.2515 0.21258         0.80448          
#Total        21    5.2848                 1.00000   
permutest(betadisper(jc.plob.F, plob.data.nrF$seq_platform, type = "centroid"))
#       Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.00147 0.0014702 0.1629    999  0.698
#Residuals 20 0.18055 0.0090275    
un.plob.F <- phyloseq::distance(plob.nrF, method = "uunifrac")
adonis(un.plob.F ~ seq_platform, strata = plob.data.nrF$sample_label, data = plob.data.nrF)
#      Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    0.8918 0.89185  4.5416 0.18506  0.001 ***
#Residuals    20    3.9274 0.19637         0.81494           
#Total        21    4.8193                 1.00000           

permutest(betadisper(un.plob.F, plob.data.nrF$seq_platform, type = "centroid"))
#       Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
#Groups     1 0.135951 0.135951 29.138    999  0.001 ***
#Residuals 20 0.093314 0.004666    
w.plob.F <- phyloseq::distance(plob.raF, method = "wunifrac")
adonis(w.plob.F ~ seq_platform, strata = plob.data.raF$sample_label, data = plob.data.raF)
#         Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    0.3846 0.38465  1.9364 0.08827  0.001 ***
#Residuals    20    3.9727 0.19864         0.91173           
#Total        21    4.3574                 1.00000   
permutest(betadisper(w.plob.F, plob.data.raF$seq_platform, type = "centroid"))
#      Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
#Groups     1 0.09654 0.096545 3.6505    999  0.076 .
#Residuals 20 0.52895 0.026447    
bc.F <- phyloseq::distance(mfol.raF, method = "bray")
adonis(bc.F ~ seq_platform, strata = mfol.data.raF$sample_label, data = mfol.data.raF)
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#seq_platform  1    0.4510 0.45103  1.3886 0.06492  0.006 **
#Residuals    20    6.4962 0.32481         0.93508          
#Total        21    6.9472                 1.00000        
permutest(betadisper(bc.F, mfol.data.raF$seq_platform, type = "centroid"))
#          Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.005343 0.0053432 0.8905    999  0.366
#Residuals 20 0.120012 0.0060006
bc.plob.F <- phyloseq::distance(plob.raF, method = "bray")
adonis(bc.plob.F ~ seq_platform, strata = plob.data.raF$sample_label, data = plob.data.raF)
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.6664 0.66641  2.2742 0.1021  0.004 **
#Residuals    20    5.8605 0.29303         0.8979          
#Total        21    6.5269                 1.0000  
permutest(betadisper(bc.plob.F, plob.data.raF$seq_platform, type = "centroid"))
#          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
#Groups     1 0.18544 0.18544 6.8682    999  0.014 *
#Residuals 20 0.54000 0.02700  


library(DESeq2)
#Perform the analysis (You need to use raw counts for DESeq, not Relative Abundances)

#Mfol
platform <-  phyloseq_to_deseq2(mfol.nr, ~ seq_platform)
gm_mean <-  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <-  apply(counts(platform), 1, gm_mean)
platform <-  estimateSizeFactors(platform, geoMeans = geoMeans)
platform <-  DESeq(platform, fitType="local")
test <-  DESeq2::DESeq(platform, test="Wald", fitType="parametric")
res <-  DESeq2::results(test, cooksCutoff = FALSE)
alpha <- 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq.nr)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)
View(sigtab)
res1 <- DESeq2::results(test, contrast = c("seq_platform", "HiSeq","MiSeq"))
res1 <- res[which(res$padj < alpha),]
res1 <-  cbind(as(res1, "data.frame"), as(tax_table(physeq.nr)[rownames(res1), ], "matrix"))

write.csv(res1, '~/Dropbox/MacBook_Transfer/comparison/mfol_sigtab_97.csv')

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
p.mfol.d2 <- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

#Plob
platform <-  phyloseq_to_deseq2(plob.nr, ~ seq_platform)
gm_mean <-  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <-  apply(counts(platform), 1, gm_mean)
platform <-  estimateSizeFactors(platform, geoMeans = geoMeans)
#platform <-  DESeq(platform, fitType="local")
test <-  DESeq2::DESeq(platform, test="Wald", fitType="parametric")
res <-  DESeq2::results(test, cooksCutoff = FALSE)
alpha <- 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq.nr)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)
View(sigtab)
res1 <- DESeq2::results(test, contrast = c("seq_platform", "HiSeq","MiSeq"))
res1 <- res[which(res$padj < alpha),]
res1 <-  cbind(as(res1, "data.frame"), as(tax_table(physeq.nr)[rownames(res1), ], "matrix"))

write.csv(res1, '~/Dropbox/MacBook Transfer/comparison/plob_sigtab_97.csv')

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
p.plob.d2 <- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggsave("~/Desktop/deseq_plob.pdf")
