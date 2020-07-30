##ASV Statistical Analyses for comparison MS###

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
physeq <- qza_to_phyloseq("~/Desktop/comparison/table_merged_filtered.qza", 
                          "~/Desktop/comparison/rooted-tree.qza", 
                          "~/Desktop/comparison/tax-without-spaces.qza", 
                          "~/Desktop/comparison/metadata_comp.txt")


#Check the taxonomic classification is correct
rank_names(physeq)
tax_table(physeq)

#Remove eukaryotes and Unassigned reads
physeq <- subset_taxa(physeq, Kingdom != "D_0__Eukaryota")
physeq <- subset_taxa(physeq, Kingdom != "Eukaryota")
physeq <- subset_taxa(physeq, Kingdom != "Unassigned")
physeq <- subset_samples(physeq, shared == "yes")

#Check for mitochondria & blast seqs to make sure you are not removing bacterial taxa
mito <- subset_taxa(physeq, Family == "Mitochondria")
mito_taxa <- rownames(otu_table(mito))
write.csv(mito_taxa, "~/Desktop/comparison/mito_asvs.csv")
#Here check ASVs on BLAST for mitochondrial vs. bacterial reads
#Keep mitochondrial reads in CSV and re-upload as "bad taxa" to remove
bad.taxa <- read.csv("~/Desktop/comparison/bad_mito_asvs.csv", header = FALSE)
bad.taxa <- levels(bad.taxa$V1)
#remove them from the dataset
all.taxa <- taxa_names(physeq)
all.taxa <- all.taxa[!(all.taxa %in% bad.taxa)]
physeq.nm <- prune_taxa(all.taxa, physeq)

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
tax_table(physeq.noncont)  #Total ASVs:11057
sample_data(physeq.noncont) #Total samples:82
#Count Sequences
totalreads <- sample_sums(physeq.noncont)
sum(totalreads) #Total reads:1566695

#Remove blanks because we already dealt with potential contaminants
physeq.noncont <- subset_samples(physeq.noncont, is.neg != "TRUE")

#Prune singletons (these are reads that are only found once)
physeq.prune <- prune_taxa(taxa_sums(physeq.noncont) > 1, physeq.noncont)

saveRDS(physeq.prune, "~/Dropbox/MacBook_Transfer/comparison/physeq_prune_asv.RDS")
physeq.prune <- readRDS("~/Dropbox/MacBook_Transfer/comparison/physeq_prune_asv.RDS")

#Unrarefied - cut out every sample with below 1000 reads
physeq.nr <- prune_samples(sample_sums(physeq.prune)>=1000, physeq.prune)
bad.samples <- c("KI15BFMD107", "KI15BFMD108", "KI15BFMD119", "KI15BFMD123", "KI15BFMD135", "KI15BFMD145", "KI15BFMD195", "KI15BFMD196", "KI15BFMD208", "KI15BFMD232_EMP")
all.samples <- sample_names(physeq.nr)
all.samples <- all.samples[!(all.samples %in% bad.samples)]
physeq.nr <- prune_samples(all.samples, physeq.nr)
data.nr <- as(sample_data(physeq.nr), "data.frame")

tax_table(physeq.nr)  #Total ASVs:5512
sample_data(physeq.nr) #Total samples:48
#Count Sequences
totalreads <- sample_sums(physeq.nr)
sum(totalreads) #1444493

#Rarefied to 1000 reads
physeq.r <- rarefy_even_depth(physeq.prune, sample.size = 1000, rngseed = 711) 
bad.samples.r <- c("KI15BFMD107", "KI15BFMD108", "KI15BFMD119", "KI15BFMD123", "KI15BFMD135", "KI15BFMD145", "KI15BFMD195", "KI15BFMD196", "KI15BFMD208", "KI15BFMD232_EMP")
all.samples.r <- sample_names(physeq.r)
all.samples.r <- all.samples.r[!(all.samples.r %in% bad.samples.r)]
physeq.r <- prune_samples(all.samples.r, physeq.r)
data.r <- as(sample_data(physeq.r), "data.frame")

#Relative abundance
physeq.ra <- transform_sample_counts(physeq.nr, function(x) x/ sum(x)) 
data.ra <- as(sample_data(physeq.ra), "data.frame")




#First check barplots by Phylum
#Re order factor levels 

p <- plot_bar(physeq.ra, x = "sample_label", fill = "Phylum")

p$data$sample_label <- factor(x = p$data$sample_label, levels = c("KI15BFMD082", "KI15BFMD089", "KI15BFMD124", "KI15BFMD142",
                                                                  "KI15BFMD144", "KI15BFMD146", "KI15BFMD190", "KI15BFMD193",
                                                                  "KI15BFMD221", "KI15BFMD228", "KI15BFMD240", "KI15BFMD079",
                                                                  "KI15BFMD086", "KI15BFMD091", "KI15BFMD105", "KI15BFMD112", 
                                                                  "KI15BFMD115", "KI15BFMD122", "KI15BFMD141", "KI15BFMD189",
                                                                  "KI15BFMD201", "KI15BFMD229", "KI15BFMD230", "KI15BFMD235"))

nb.cols <- 41
mycolors <- colorRampPalette(brewer.pal(11, "RdYlBu"))(nb.cols)
p1 <- p + geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = 'stack', width = 0.7) +
  #guides(fill = FALSE, color = FALSE) +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) +
  ylab("Relative Abundance") +
  xlab("Sample Label") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


p1 + facet_grid(rows = vars(seq_platform)) 


ggsave("~/Desktop/barplot_phy_asvs.pdf")



#Alpha diversity
#Observed Sp. Richness and Shannon
erich <- estimate_richness_wPD(physeq.r, measures = c("Observed", "Shannon", "FaithPD"))
erich$platform <- data.r$seq_platform
erich$sample <- data.r$sample_label
erich$species <- data.r$field_host_name

erich.mfol <- erich %>% filter(species == "Montipora foliosa")

erich.plob <- erich %>% filter(species == "Porites lobata")


#Model fit for Alpha diversity 

library(nlme)
library(sjPlot)
library(effects)
qq.line = function(x) {
  # following four lines from base R's qqline()
  y <- quantile(x[!is.na(x)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  return(c(int = int, slope = slope))
}

#mfol - Observed
sqOb.model.a <- lme(log(Observed) ~ platform, random = ~1|sample, data = erich.mfol, method = "REML")

plot(sqOb.model.a)
qq.line(resid(sqOb.model.a))
plot_grid(plot_model(sqOb.model.a, type = "diag"))
plot(allEffects(sqOb.model.a))
plot_model(sqOb.model.a, type = "eff", terms = "platform")
anova(sqOb.model.a)
summary(sqOb.model.a)

newdata = data.frame(platform = levels(erich.mfol$platform))
Xmat = model.matrix(~platform, data = newdata)
coefs = fixef(sqOb.model.a)
fit = as.vector(coefs %*% t(Xmat))
se = sqrt(diag(Xmat %*% vcov(sqOb.model.a) %*% t(Xmat)))
q = qt(0.975, df = nrow(sqOb.model.a$data) - length(coefs) -
         2)
newdata = cbind(newdata, fit = fit, lower = fit - q * se, upper = fit +
                  q * se)
mfol.ob <- ggplot(newdata, aes(y = fit, x = platform)) + geom_pointrange(aes(ymin = lower,
                                                                             ymax = upper)) + scale_y_continuous("log(Observed Richness)") + theme_classic()


#plob - Observed
sqOb.model.b <- lme(log(Observed) ~ platform, random = ~1|sample, data = erich.plob, method = "REML")

plot(sqOb.model.b)
qq.line(resid(sqOb.model.b))
plot_grid(plot_model(sqOb.model.b, type = "diag"))
plot(allEffects(sqOb.model.b))
plot_model(sqOb.model.b, type = "eff", terms = "platform")
anova(sqOb.model.b)
summary(sqOb.model.b)

newdata = data.frame(platform = levels(erich.plob$platform))
Xmat = model.matrix(~platform, data = newdata)
coefs = fixef(sqOb.model.b)
fit = as.vector(coefs %*% t(Xmat))
se = sqrt(diag(Xmat %*% vcov(sqOb.model.b) %*% t(Xmat)))
q = qt(0.975, df = nrow(sqOb.model.b$data) - length(coefs) -
         2)
newdata = cbind(newdata, fit = fit, lower = fit - q * se, upper = fit +
                  q * se)
plob.ob <- ggplot(newdata, aes(y = fit, x = platform)) + geom_pointrange(aes(ymin = lower,
                                                                             ymax = upper)) + scale_y_continuous("log(Observed Richness)") + theme_classic()


#Mfol - Shannon
model.c <- lme(Shannon ~ platform, random = ~1|sample, data = erich.mfol, method = "REML")

plot(model.c)
qq.line(resid(model.c))
plot_grid(plot_model(model.c, type = "diag"))
plot(allEffects(model.c))
plot_model(model.c, type = "eff", terms = "platform")
anova(model.c)
summary(model.c)

newdata = data.frame(platform = levels(erich.mfol$platform))
Xmat = model.matrix(~platform, data = newdata)
coefs = fixef(model.c)
fit = as.vector(coefs %*% t(Xmat))
se = sqrt(diag(Xmat %*% vcov(model.c) %*% t(Xmat)))
q = qt(0.975, df = nrow(model.c$data) - length(coefs) -
         2)
newdata = cbind(newdata, fit = fit, lower = fit - q * se, upper = fit +
                  q * se)
mfol.sh <- ggplot(newdata, aes(y = fit, x = platform)) + geom_pointrange(aes(ymin = lower,
                                                                             ymax = upper)) + scale_y_continuous("Shannon Diversity Index") + theme_classic()



#Plob - Shannon

model.d <- lme(Shannon ~ platform, random = ~1|sample, data = erich.plob, method = "REML")

plot(model.d)
qq.line(resid(model.d))
plot_grid(plot_model(model.d, type = "diag"))
plot(allEffects(model.d))
plot_model(model.d, type = "eff", terms = "platform")
anova(model.d)
summary(model.d)

newdata = data.frame(platform = levels(erich.plob$platform))
Xmat = model.matrix(~platform, data = newdata)
coefs = fixef(model.d)
fit = as.vector(coefs %*% t(Xmat))
se = sqrt(diag(Xmat %*% vcov(model.d) %*% t(Xmat)))
q = qt(0.975, df = nrow(model.d$data) - length(coefs) -
         2)
newdata = cbind(newdata, fit = fit, lower = fit - q * se, upper = fit +
                  q * se)
plob.sh <- ggplot(newdata, aes(y = fit, x = platform)) + geom_pointrange(aes(ymin = lower,
                                                                             ymax = upper)) + scale_y_continuous("Shannon Diversity Index") + theme_classic()


#Chao1

erich.nr <- estimate_richness(physeq.nr, measures = "Chao1")
erich.nr$platform <- data.nr$seq_platform
erich.nr$sample <- data.nr$sample_label
erich.nr$species <- data.nr$field_host_name

erich.nr.mfol <- erich.nr %>% filter(species == "Montipora foliosa")
erich.nr.plob <- erich.nr %>% filter(species == "Porites lobata")

hist(erich.nr$Chao1)



#Mfol- Chao1
model.e <- lme(log(Chao1) ~ platform, random = ~1|sample, data = erich.nr.mfol, method = "REML")

plot(model.e)
qq.line(resid(model.e))
plot_grid(plot_model(model.e, type = "diag"))
plot(allEffects(model.e))
plot_model(model.e, type = "eff", terms = "platform")
anova(model.e)
summary(model.e)

newdata = data.frame(platform = levels(erich.nr.mfol$platform))
Xmat = model.matrix(~platform, data = newdata)
coefs = fixef(model.e)
fit = as.vector(coefs %*% t(Xmat))
se = sqrt(diag(Xmat %*% vcov(model.e) %*% t(Xmat)))
q = qt(0.975, df = nrow(model.e$data) - length(coefs) -
         2)
newdata = cbind(newdata, fit = fit, lower = fit - q * se, upper = fit +
                  q * se)
mfol.ch <- ggplot(newdata, aes(y = fit, x = platform)) + geom_pointrange(aes(ymin = lower,
                                                                             ymax = upper)) + scale_y_continuous("log(Chao1 Diversity Index)") + theme_classic()

#P.lob - Chao 1
model.f <- lme(log(Chao1) ~ platform, random = ~1|sample, data = erich.nr.plob, method = "REML")

plot(model.f)
qq.line(resid(model.f))
plot_grid(plot_model(model.f, type = "diag"))
plot(allEffects(model.f))
plot_model(model.f, type = "eff", terms = "platform")
anova(model.f)
summary(model.f)

newdata = data.frame(platform = levels(erich.nr.plob$platform))
Xmat = model.matrix(~platform, data = newdata)
coefs = fixef(model.f)
fit = as.vector(coefs %*% t(Xmat))
se = sqrt(diag(Xmat %*% vcov(model.f) %*% t(Xmat)))
q = qt(0.975, df = nrow(model.f$data) - length(coefs) -
         2)
newdata = cbind(newdata, fit = fit, lower = fit - q * se, upper = fit +
                  q * se)
plob.ch <- ggplot(newdata, aes(y = fit, x = platform)) + geom_pointrange(aes(ymin = lower,
                                                                             ymax = upper)) + scale_y_continuous("log(Chao1 Diversity Index)") + theme_classic()


##Faith's PD

#MFol FaithsD
model.7 <- lme(FaithPD ~ platform, random = ~1|sample, data = erich.mfol, method = "REML")
model.7a <- lme(log(FaithPD) ~ platform, random = ~1|sample, data = erich.mfol, method = "REML")
model.7b <- lme(sqrt(FaithPD) ~ platform, random = ~1|sample, data = erich.mfol, method = "REML")
AIC(model.7, model.7a, model.7b) #use model.7a 

plot(model.7a)
qq.line(resid(model.7a))
plot_grid(plot_model(model.7a, type = "diag"))
plot(allEffects(model.7a))
plot_model(model.7a, type = "eff", terms = "platform")
anova(model.7a)
summary(model.7a)

#numDF denDF   F-value p-value
#(Intercept)     1    10 261.93513  <.0001
#platform        1    10   0.71028  0.4191

newdata = data.frame(platform = levels(erich.mfol$platform))
Xmat = model.matrix(~platform, data = newdata)
coefs = fixef(model.7a)
fit = as.vector(coefs %*% t(Xmat))
se = sqrt(diag(Xmat %*% vcov(model.7a) %*% t(Xmat)))
q = qt(0.975, df = nrow(model.7a$data) - length(coefs) -
         2)
newdata = cbind(newdata, fit = fit, lower = fit - q * se, upper = fit +
                  q * se)
mfol.fpd <- ggplot(newdata, aes(y = fit, x = platform)) + geom_pointrange(aes(ymin = lower,
                                                                              ymax = upper)) + scale_y_continuous("log(Faith's PD)") + theme_classic()


#Plob FaithsD
model.8 <- lme(FaithPD ~ platform, random = ~1|sample, data = erich.plob, method = "REML")
model.8a <- lme(log(FaithPD) ~ platform, random = ~1|sample, data = erich.plob, method = "REML")
model.8b <- lme(sqrt(FaithPD) ~ platform, random = ~1|sample, data = erich.plob, method = "REML")
AIC(model.8, model.8a, model.8b) #use model.8a 

plot(model.8a)
qq.line(resid(model.8a))
plot_grid(plot_model(model.8a, type = "diag"))
plot(allEffects(model.8a))
plot_model(model.8a, type = "eff", terms = "platform")
anova(model.8a)
summary(model.8a)

#numDF denDF   F-value p-value
#(Intercept)     1    12 250.54081  <.0001
#platform        1    12   2.50525  0.1395

newdata = data.frame(platform = levels(erich.plob$platform))
Xmat = model.matrix(~platform, data = newdata)
coefs = fixef(model.8a)
fit = as.vector(coefs %*% t(Xmat))
se = sqrt(diag(Xmat %*% vcov(model.8a) %*% t(Xmat)))
q = qt(0.975, df = nrow(model.8a$data) - length(coefs) -
         2)
newdata = cbind(newdata, fit = fit, lower = fit - q * se, upper = fit +
                  q * se)
plob.fpd <- ggplot(newdata, aes(y = fit, x = platform)) + geom_pointrange(aes(ymin = lower,
                                                                              ymax = upper)) + scale_y_continuous("log(Faith's PD)") + theme_classic()

###Betadiversity###
##All results from PERMANOVA are copied below tests after hashtags##

#Split data by species

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
#adonis(formula = jc.mfol ~ seq_platform, data = mfol.data.nr,      strata = mfol.data.nr$sample_label) 

#Blocks:  strata 
#Permutation: free
#Number of permutations: 2047

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.5378  0.5378  1.1709 0.05531  0.002 **
#Residuals    20    9.1861  0.4593         0.94469          
#Total        21    9.7239                 1.00000
permutest(betadisper(jc.mfol, mfol.data.nr$seq_platform, type = "centroid"))
#Response: Distances
#Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.0000249 2.493e-05 0.1294    999  0.731
#Residuals 20 0.0038519 1.926e-04  

#weighted Unifrac (presence absence + abundance + phylogeny)
wu.mfol <- phyloseq::distance(mfol.ra, method = "wunifrac")
adonis(wu.mfol ~ seq_platform, strata = mfol.data.ra$sample_label, data = mfol.data.ra)
#Df SumsOfSqs   MeanSqs F.Model     R2 Pr(>F)
#seq_platform  1  0.003830 0.0038305 0.85294 0.0409  0.293
#Residuals    20  0.089819 0.0044910         0.9591       
#Total        21  0.093650                   1.0000    
permutest(betadisper(wu.mfol, mfol.data.ra$seq_platform, type = "centroid"))
#Df   Sum Sq    Mean Sq     F N.Perm Pr(>F)
#Groups     1 0.000192 0.00019204 0.192    999  0.685
#Residuals 20 0.020002 0.00100010     


#Porites lobata
jc.plob <- phyloseq::distance(plob.nr, method = "jaccard", binary = TRUE)
adonis(jc.plob ~ seq_platform, strata = plob.data.nr$sample_label, data = plob.data.nr)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    0.6409 0.64094  1.4444 0.05677  0.001 ***
#  Residuals    24   10.6502 0.44376         0.94323           
#Total        25   11.2911                 1.00000           
permutest(betadisper(jc.plob, plob.data.nr$seq_platform, type = "centroid"))
#Response: Distances
#Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.0008659 0.00086589 0.8327    999  0.371
#Residuals 24 0.0249557 0.00103982  

##Weighted Unifrac
wu.plob <- phyloseq::distance(plob.ra, method = "wunifrac")
adonis(wu.plob ~ seq_platform, strata = plob.data.ra$sample_label, data = plob.data.ra)
# Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1  0.008345 0.0083453 0.71903 0.02909  0.002 **
#  Residuals    24  0.278552 0.0116063         0.97091          
#Total        25  0.286897                   1.00000 
permutest(betadisper(wu.plob, plob.data.ra$seq_platform, type = "centroid"))
#Response: Distances
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.001526 0.0015259 0.2821    999  0.604
#Residuals 24 0.129820 0.0054092  


##Bray curtis - uses relative abundance to account for differences in sequencing depth
#Montipora foliosa
bc.mfol <- phyloseq::distance(mfol.ra, method = "bray")
adonis(bc.mfol ~ seq_platform, strata = mfol.data.ra$sample_label, data = mfol.data.ra)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.4907 0.49074  1.1402 0.05394  0.007 **
#  Residuals    20    8.6078 0.43039         0.94606          
#Total        21    9.0986                 1.00000  
permutest(betadisper(bc.mfol, mfol.data.ra$seq_platform, type = "centroid"))
#Response: Distances
#Df    Sum Sq    Mean Sq     F N.Perm Pr(>F)
#Groups     1 0.0000003 0.00000031 4e-04    999  0.986
#Residuals 20 0.0175232 0.00087616 

#Unweighted unifrac
un.mfol <- phyloseq::distance(mfol.nr, method = "uunifrac")
adonis(un.mfol ~ seq_platform, strata = mfol.data.nr$sample_label, data = mfol.data.nr)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.7417 0.74171  2.1399 0.09666  0.006 **
#  Residuals    20    6.9321 0.34660         0.90334          
#Total        21    7.6738                 1.00000          

permutest(betadisper(un.mfol, mfol.data.nr$seq_platform, type = "centroid"))
#Response: Distances
#Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)   
#Groups     1 0.058683 0.058683 10.386    999  0.003 **
#  Residuals 20 0.113011 0.005651  

#Porites lobata
bc.plob <- phyloseq::distance(plob.ra, method = "bray")
adonis(bc.plob ~ seq_platform, strata = plob.data.ra$sample_label, data = plob.data.ra)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    1.1826 1.18256   3.408 0.12434  0.004 **
#  Residuals    24    8.3278 0.34699         0.87566          
#Total        25    9.5104                 1.00000  
permutest(betadisper(bc.plob, plob.data.ra$seq_platform, type = "centroid"))
#Response: Distances
#Df  Sum Sq  Mean Sq     F N.Perm Pr(>F)  
#Groups     1 0.18479 0.184793 7.585    999  0.012 *
#  Residuals 24 0.58471 0.024363 

un.plob <- phyloseq::distance(plob.nr, method = "uunifrac")
adonis(un.plob ~ seq_platform, strata = plob.data.nr$sample_label, data = plob.data.nr)
#       Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    0.5186 0.51862  1.4754 0.05791  0.001 ***
#  Residuals    24    8.4363 0.35151         0.94209           
#Total        25    8.9549                 1.00000             
permutest(betadisper(un.plob, plob.data.nr$seq_platform, type = "centroid"))
#Response: Distances
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.001417 0.0014169 0.2161    999  0.645
#Residuals 24 0.157353 0.0065564  


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
  ggtitle("M. aequituberculata Jaccard") +
  theme_classic()
ggsave("~/Desktop/jc_mfol_asv_nmds.pdf")

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
  ggtitle("P. lobata Jaccard") +
  theme_classic()
ggsave("~/Desktop/jc_plob_asv_nmds.pdf")

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
  ggtitle("M. aequituberculata Bray Curtis") +
  theme_classic()
ggsave("~/Desktop/br_mfol_asv_nmds.pdf")

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
  ggtitle("P. lobata Bray Curtis") +
  theme_classic()
ggsave("~/Desktop/br_plob_asv_nmds.pdf")


##Unweighted Unifrac
#Montipora foliosa
ord.un.mfol <- ordinate(mfol.ra, "NMDS", "uunifrac") #trymax is not accepted here
stressplot(ord.un.mfol)
scores.un.mfol <- as.data.frame(scores(ord.un.mfol))
scores.un.mfol$platform <- mfol.data.ra$seq_platform
scores.un.mfol$sample <- mfol.data.ra$sample_label
ggplot(scores.un.mfol, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = platform, shape = platform), size = 4) +
  #geom_text(label = rownames(scores), nudge_x = 0.05, nudge_y = 0.05,check_overlap = T) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = platform), linetype = 2) +
  ggtitle("M. aequituberculata Unweighted UniFrac") +
  theme_classic()
ggsave("~/Desktop/un_mfol_asv_nmds.pdf")

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
  ggtitle("P. lobata Unweighted UniFrac") +
  theme_classic()
ggsave("~/Desktop/un_plob_asv_nmds.pdf")


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
  ggtitle("M. aequituberculata Weighted UniFrac") +
  theme_classic()
ggsave("~/Desktop/wun_mfol_asv_nmds.pdf")

#Porites lobata
ord.wun.plob <- ordinate(plob.ra, "NMDS", "wunifrac") 
stressplot(ord.wun.plob)
scores.wun.plob <- as.data.frame(scores(ord.wun.plob))
scores.wun.plob$platform <- plob.data.ra$seq_platform
scores.wun.plob$sample <- plob.data.ra$sample_label
ggplot(scores.wun.plob, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = platform, shape = platform), size = 4) +
  #geom_text(label = rownames(scores), nudge_x = 0.05, nudge_y = 0.05,check_overlap = T) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = platform), linetype = 2) +
  ggtitle("P. lobata Weighted UniFrac") +
  theme_classic()
ggsave("~/Desktop/wun_plob_asv_nmds.pdf")


##Check betadiversity on Family level

#M fol

mfol.fam.ra <- mfol.ra %>% 
  tax_glom(taxrank = "Family") 
mfol.fam.ra.data <- as(sample_data(mfol.fam.ra), "data.frame")

bc.mfol.fam <- phyloseq::distance(mfol.fam.ra, method = "bray")
adonis(bc.mfol.fam ~ seq_platform, strata = mfol.fam.ra.data$sample_label, data = mfol.fam.ra.data)
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#seq_platform  1    0.3531 0.35309  1.2421 0.05847  0.037 *
#  Residuals    20    5.6854 0.28427         0.94153         
#Total        21    6.0385                 1.00000 
permutest(betadisper(bc.mfol.fam, mfol.fam.ra.data$seq_platform, type = "centroid"))
#        Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.000547 0.0005470 0.0655    999  0.781
#Residuals 20 0.167041 0.0083521 

w.mfol.fam <- phyloseq::distance(mfol.fam.ra, method = "wunifrac")
adonis(w.mfol.fam ~ seq_platform, strata = mfol.fam.ra.data$sample_label, data = mfol.fam.ra.data)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#seq_platform  1    0.1785  0.1785  1.0785 0.05117  0.162
#Residuals    20    3.3101  0.1655         0.94883       
#Total        21    3.4886                 1.00000  
permutest(betadisper(w.mfol.fam, mfol.fam.ra.data$seq_platform, type = "centroid"))
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.006115 0.0061148 0.7127    999   0.38
#Residuals 20 0.171586 0.0085793   

mfol.fam.nr <- mfol.nr %>% 
  tax_glom(taxrank = "Family") 
mfol.fam.nr.data <- as(sample_data(mfol.fam.nr), "data.frame")

jc.mfol.fam <- phyloseq::distance(mfol.fam.nr, method = "jaccard", binary = TRUE)
adonis(jc.mfol.fam ~ seq_platform, strata = mfol.fam.nr.data$sample_label, data = mfol.fam.nr.data)
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.4548 0.45484  1.5133 0.07034  0.002 **
#  Residuals    20    6.0115 0.30057         0.92966          
#Total        21    6.4663                 1.00000  
permutest(betadisper(jc.mfol.fam, mfol.fam.nr.data$seq_platform, type = "centroid"))
#        Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.001669 0.0016685 0.6264    999  0.456
#Residuals 20 0.053272 0.0026636 

un.mfol.fam <- phyloseq::distance(mfol.fam.nr, method = "unifrac")
adonis(un.mfol.fam ~ seq_platform, strata = mfol.fam.nr.data$sample_label, data = mfol.fam.nr.data)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.3378 0.33784  1.5157 0.07045  0.003 **
#  Residuals    20    4.4578 0.22289         0.92955          
#Total        21    4.7956                 1.00000  
permutest(betadisper(un.mfol.fam, mfol.fam.nr.data$seq_platform, type = "centroid"))
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.000850 0.0008499 0.1875    999  0.751
#Residuals 20 0.090665 0.0045332 



plob.fam.ra <- plob.ra %>% 
  tax_glom(taxrank = "Family") 
plob.fam.ra.data <- as(sample_data(plob.fam.ra), "data.frame")

bc.plob.fam <- phyloseq::distance(plob.fam.ra, method = "bray")
adonis(bc.plob.fam ~ seq_platform, strata = plob.fam.ra.data$sample_label, data = plob.fam.ra.data)
#   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    0.9876 0.98755  3.2616 0.11964  0.001 ***
#  Residuals    24    7.2667 0.30278         0.88036           
#Total        25    8.2543                 1.00000  
permutest(betadisper(bc.plob.fam, plob.fam.ra.data$seq_platform, type = "centroid"))
#         Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)   
#Groups     1 0.23896 0.238961 9.4064    999  0.002 **
#  Residuals 24 0.60970 0.025404 

w.plob.fam <- phyloseq::distance(plob.fam.ra, method = "wunifrac")
adonis(w.plob.fam ~ seq_platform, strata = plob.fam.ra.data$sample_label, data = plob.fam.ra.data)
#Df SumsOfSqs   MeanSqs F.Model     R2 Pr(>F)   
#seq_platform  1  0.007394 0.0073943  2.3533 0.0893  0.004 **
#  Residuals    24  0.075409 0.0031421         0.9107          
#Total        25  0.082803                   1.0000  
permutest(betadisper(w.plob.fam, plob.fam.ra.data$seq_platform, type = "centroid"))
#      Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)  
#Groups     1 0.0019175 0.00191746 4.7047    999  0.045 *
#  Residuals 24 0.0097814 0.00040756        

plob.fam.nr <- plob.nr %>% 
  tax_glom(taxrank = "Family") 
plob.fam.nr.data <- as(sample_data(plob.fam.nr), "data.frame")

jc.plob.fam <- phyloseq::distance(plob.fam.nr, method = "jaccard", binary = TRUE)
adonis(jc.plob.fam ~ seq_platform, strata = plob.fam.nr.data$sample_label, data = plob.fam.nr.data)
#     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.6039 0.60393   1.775 0.06887  0.003 **
#  Residuals    24    8.1658 0.34024         0.93113          
#Total        25    8.7697                 1.00000   
permutest(betadisper(jc.plob.fam, plob.fam.nr.data$seq_platform, type = "centroid"))
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)    
#Groups     1 0.021761 0.0217605 15.932    999  0.001 ***
#  Residuals 24 0.032780 0.0013658      


un.plob.fam <- phyloseq::distance(plob.fam.nr, method = "unifrac")
adonis(un.plob.fam ~ seq_platform, strata = plob.fam.nr.data$sample_label, data = plob.fam.nr.data)
#     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.4218 0.42182  2.2349 0.08519  0.004 **
#  Residuals    24    4.5298 0.18874         0.91481          
#Total        25    4.9516                 1.00000   
permutest(betadisper(un.plob.fam, plob.fam.nr.data$seq_platform, type = "centroid"))
#      Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
#Groups     1 0.019617 0.0196172 5.4089    999  0.035 *
#  Residuals 24 0.087043 0.0036268        


##By PHylum

mfol.phy.ra <- mfol.ra %>% 
  tax_glom(taxrank = "Phylum") 
mfol.phy.ra.data <- as(sample_data(mfol.phy.ra), "data.frame")

bc.mfol.phy <- phyloseq::distance(mfol.phy.ra, method = "bray")
adonis(bc.mfol.phy ~ seq_platform, strata = mfol.phy.ra.data$sample_label, data = mfol.phy.ra.data)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#seq_platform  1   0.09778 0.09778 0.92318 0.04412  0.295
#Residuals    20   2.11833 0.10592         0.95588       
#Total        21   2.21611                 1.00000   
permutest(betadisper(bc.mfol.phy, mfol.phy.ra.data$seq_platform, type = "centroid"))
#    Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.016646 0.016646 1.6765    999  0.226
#Residuals 20 0.198580 0.009929  

w.mfol.phy <- phyloseq::distance(mfol.phy.ra, method = "wunifrac")
adonis(w.mfol.phy ~ seq_platform, strata = mfol.phy.ra.data$sample_label, data = mfol.phy.ra.data)
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
#seq_platform  1   0.06191 0.061912 0.84573 0.04057  0.422
#Residuals    20   1.46412 0.073206         0.95943       
#Total        21   1.52603                  1.00000    
permutest(betadisper(w.mfol.phy, mfol.phy.ra.data$seq_platform, type = "centroid"))
#      Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.027319 0.0273188 2.8582    999  0.121
#Residuals 20 0.191161 0.0095581   

mfol.phy.nr <- mfol.nr %>% 
  tax_glom(taxrank = "Phylum") 
mfol.phy.nr.data <- as(sample_data(mfol.phy.nr), "data.frame")

jc.mfol.phy <- phyloseq::distance(mfol.phy.nr, method = "jaccard", binary = TRUE)
adonis(jc.mfol.phy ~ seq_platform, strata = mfol.phy.nr.data$sample_label, data = mfol.phy.nr.data)
#         Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#seq_platform  1    0.2565 0.25650  2.3673 0.10584  0.021 *
#  Residuals    20    2.1670 0.10835         0.89416         
#Total        21    2.4235                 1.00000    
permutest(betadisper(jc.mfol.phy, mfol.phy.nr.data$seq_platform, type = "centroid"))
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.000193 0.0001935 0.0392    999  0.859
#Residuals 20 0.098820 0.0049410   

un.mfol.phy <- phyloseq::distance(mfol.phy.nr, method = "unifrac")
adonis(un.mfol.phy ~ seq_platform, strata = mfol.phy.nr.data$sample_label, data = mfol.phy.nr.data)
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
#seq_platform  1   0.23653 0.236532  2.8626 0.12521  0.024 *
#  Residuals    20   1.65255 0.082627         0.87479         
#Total        21   1.88908                  1.00000         
permutest(betadisper(un.mfol.phy, mfol.phy.nr.data$seq_platform, type = "centroid"))
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.000269 0.0002685 0.0499    999  0.848
#Residuals 20 0.107726 0.0053863 

plob.phy.ra <- plob.ra %>% 
  tax_glom(taxrank = "Phylum") 
plob.phy.ra.data <- as(sample_data(plob.fam.ra), "data.frame")

bc.plob.phy <- phyloseq::distance(plob.phy.ra, method = "bray")
adonis(bc.plob.phy ~ seq_platform, strata = plob.phy.ra.data$sample_label, data = plob.phy.ra.data)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    0.6428 0.64283  4.8871 0.16918  0.001 ***
#  Residuals    24    3.1569 0.13154         0.83082           
#Total        25    3.7997                 1.00000   
permutest(betadisper(bc.plob.phy, plob.phy.ra.data$seq_platform, type = "centroid"))
#Df  Sum Sq Mean Sq      F N.Perm Pr(>F)    
#Groups     1 0.45774 0.45774 19.573    999  0.001 ***
#  Residuals 24 0.56129 0.02339  


w.plob.phy <- phyloseq::distance(plob.phy.ra, method = "wunifrac")
adonis(w.plob.phy ~ seq_platform, strata = plob.phy.ra.data$sample_label, data = plob.phy.ra.data)
#     Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
#seq_platform  1   0.16258 0.162584  2.5249 0.09519  0.012 *
#  Residuals    24   1.54543 0.064393         0.90481         
#Total        25   1.70801                  1.00000   
permutest(betadisper(w.plob.phy, plob.phy.ra.data$seq_platform, type = "centroid"))
#    Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.03950 0.039501 1.9525    999  0.192
#Residuals 24 0.48555 0.020231  

plob.phy.nr <- plob.nr %>% 
  tax_glom(taxrank = "Phylum") 
plob.phy.nr.data <- as(sample_data(plob.phy.nr), "data.frame")

jc.plob.phy <- phyloseq::distance(plob.phy.nr, method = "jaccard", binary = TRUE)
adonis(jc.plob.phy ~ seq_platform, strata = plob.phy.nr.data$sample_label, data = plob.phy.nr.data)
#Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)   
#seq_platform  1    0.4440 0.44405  2.3563 0.0894   0.01 **
#  Residuals    24    4.5227 0.18845         0.9106          
#Total        25    4.9668                 1.0000  
permutest(betadisper(jc.plob.phy, plob.phy.nr.data$seq_platform, type = "centroid"))
#Df  Sum Sq  Mean Sq     F N.Perm Pr(>F)   
#Groups     1 0.11703 0.117034 22.39    999  0.002 **
#  Residuals 24 0.12545 0.005227 

un.plob.phy <- phyloseq::distance(plob.phy.nr, method = "unifrac")
adonis(bc.plob.phy ~ seq_platform, strata = plob.phy.nr.data$sample_label, data = plob.phy.nr.data)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    0.6428 0.64283  4.8871 0.16918  0.001 ***
#  Residuals    24    3.1569 0.13154         0.83082           
#Total        25    3.7997                 1.00000   
permutest(betadisper(bc.plob.phy, plob.phy.nr.data$seq_platform, type = "centroid"))
#Df  Sum Sq Mean Sq      F N.Perm Pr(>F)    
#Groups     1 0.45774 0.45774 19.573    999  0.001 ***
#  Residuals 24 0.56129 0.02339     

#Differential Abundance Analysis

library(DESeq2)
#Perform the analysis (You need to use raw counts for DESeq, not Relative Abundances)
#Montipora foliosa
platform <-  phyloseq_to_deseq2(mfol.nr, ~ seq_platform)
gm_mean <-  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <-  apply(counts(platform), 1, gm_mean)
platform <-  estimateSizeFactors(platform, geoMeans = geoMeans)
platform <-  DESeq(platform, fitType="local")
test.mfol <-  DESeq2::DESeq(platform, test="Wald", fitType="parametric")
res.mfol <-  DESeq2::results(test.mfol, cooksCutoff = FALSE)
alpha <- 0.05
sigtab.mfol = res[which(res.mfol$padj < alpha), ]
sigtab.mfol = cbind(as(sigtab.mfol, "data.frame"), as(tax_table(mfol.nr)[rownames(sigtab.mfol), ], "matrix"))
head(sigtab.mfol)
dim(sigtab.mfol)
View(sigtab.mfol)
res1.mfol <- DESeq2::results(test.mfol, contrast = c("seq_platform", "HiSeq","MiSeq"))
res1.mfol <- res[which(res.mfol$padj < alpha),]
res1.mfol <-  cbind(as(res1.mfol, "data.frame"), as(tax_table(mfol.nr)[rownames(res1.mfol), ], "matrix"))

write.csv(sigtab.mfol, "~/Desktop/comparison/sigtab_mfol.csv")

#theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
x = tapply(sigtab.mfol$log2FoldChange, sigtab.mfol$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab.mfol$Genus = factor(as.character(sigtab.mfol$Genus), levels=names(x))
mfol.deseq <- ggplot(sigtab.mfol, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  ggtitle("Differential Abundance - M. foliosa") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) 


#Porites lobata
#Montipora foliosa
platform <-  phyloseq_to_deseq2(plob.nr, ~ seq_platform)
gm_mean <-  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <-  apply(counts(platform), 1, gm_mean)
platform <-  estimateSizeFactors(platform, geoMeans = geoMeans)
platform <-  DESeq(platform, fitType="local")
test.plob <-  DESeq2::DESeq(platform, test="Wald", fitType="parametric")
res.plob <-  DESeq2::results(test.plob, cooksCutoff = FALSE)
alpha <- 0.05
sigtab.plob = res[which(res.plob$padj < alpha), ]
sigtab.plob = cbind(as(sigtab.plob, "data.frame"), as(tax_table(plob.nr)[rownames(sigtab.plob), ], "matrix"))
head(sigtab.plob)
dim(sigtab.plob)
View(sigtab.plob)
res1.plob <- DESeq2::results(test.plob, contrast = c("seq_platform", "HiSeq","MiSeq"))
res1.plob <- res[which(res.plob$padj < alpha),]
res1.plob <-  cbind(as(res1.plob, "data.frame"), as(tax_table(plob.nr)[rownames(res1.plob), ], "matrix"))

write.csv(sigtab.plob, "~/Desktop/comparison/sigtab_plob.csv")

#theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
x = tapply(sigtab.plob$log2FoldChange, sigtab.plob$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab.plob$Genus = factor(as.character(sigtab.plob$Genus), levels=names(x))
plob.deseq <- ggplot(sigtab.plob, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  ggtitle("Differential Abundance - P. lobata") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) 




##ANCOM###

#Make sure to load the functions using the script ANCOM.script
#Split data by field host

physeq.nr.mfol <- subset_samples(physeq.nr, field_host_name == "Montipora foliosa")
physeq.nr.plob <- subset_samples(physeq.nr, field_host_name == "Porites lobata")
#First filter out any taxa that have less than 10 reads across all samples
physeq.nr.mfol <- prune_taxa(taxa_sums(physeq.nr.mfol) > 10, physeq.nr.mfol)
data.nr.mfol <- as(sample_data(physeq.nr.mfol), "data.frame")
physeq.nr.plob <- prune_taxa(taxa_sums(physeq.nr.plob) > 10, physeq.nr.plob)
data.nr.plob <- as(sample_data(physeq.nr.plob), "data.frame")

OTU.mfol <- otu_table(physeq.nr.mfol)
OTU.plob <- otu_table(physeq.nr.plob)
OTU <- otu_table(physeq.nr)

##For all 

feature_table <- OTU
meta_data <- data.nr
meta_data$Sample.ID <- as.character(meta_data$sample)
sample_var <- "Sample.ID"
group_var <- "seq_platform"
out_cut <- 0.05
zero_cut <- 0.90
lib_cut <- 1000
neg_lb <- TRUE

pre_pro <- feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                     out_cut, zero_cut, lib_cut, neg_lb)
feature_table <- pre_pro$feature_table
meta_data <- pre_pro$meta_data
struc_zero <- pre_pro$structure_zeros

#Step 2 run the ANCOM

main_var <- "seq_platform"
p_adj_method <- "BH"
alpha <- 0.05
adj_formula <- NULL
rand_formula <- NULL

res <- ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, alpha,
             adj_formula, rand_formula)

head(res$out)
View(res)
write_csv(res$out, "~/Desktop/comparison/res_ancom.csv")

#Step 3 volcano plot
#number of taxa except structural zeros
n_taxa <- ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
#Cutoff values for declaring differentially abundant taxa
cut_off <-  c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) <-  c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

#Annotation & Figure
dat_ann <- data.frame(x = min(res$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")
fig <- res$fig +
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
fig  

#find taxonomy of significantly different features(ASVs)
ancom_results_mfol <- res$out
ancom_mfol_true_0.7 <- ancom_results_mfol %>% filter(detected_0.7 == "TRUE")
ancom_results_true_0.6 <- ancom_results %>% filter(detected_0.6 == "TRUE")

tax_table <- data.frame(tax_table(physeq.nr))


tax_match = match(ancom_mfol_true_0.7$taxa_id, rownames(tax_table))
tax_match_res = tax_table[tax_match,]
head(tax_match_res)
View(tax_match_res)
write.csv(tax_match_res, "~/Desktop/comparison/tax_match.csv")
#40cca0895b09695fccc2a9b98eed3b8d Rickettsiales 
#775a42986d3850462cc4afe534671365 Vibrio


##How about RA of endozoicomonas and Vibrio?

endo.mfol <- subset_taxa(mfol.ra, Genus == "Endozoicomonas")
endo.mfol.data <- as(sample_data(endo.mfol), "data.frame")
endo.plob <- subset_taxa(plob.ra, Genus == "Endozoicomonas")
endo.plob.data <- as(sample_data(endo.plob), "data.frame")

endo.mfol.melt <- psmelt(endo.mfol)

sum.mfol.endo <- ddply(endo.mfol.melt, c("seq_platform", "Species"), summarise,
                       N = length(Abundance), 
                       mean = mean(Abundance),
                       sd = sd(Abundance), 
                       se = sd/sqrt(N)
)

ggplot(sum.mfol.endo, aes(x = seq_platform, y = mean)) +
  geom_bar(stat = "identity", aes(fill = Species, colour = Species)) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) +
  xlab("\nSequencing Platform") +
  ylab("Relative Abundance\n") +
  ggtitle("Endozoicomonas - M. foliosa") 



endo.plob.melt <- psmelt(endo.plob)

sum.plob.endo <- ddply(endo.plob.melt, c("seq_platform", "Species"), summarise,
                       N = length(Abundance), 
                       mean = mean(Abundance),
                       sd = sd(Abundance), 
                       se = sd/sqrt(N)
)

ggplot(sum.plob.endo, aes(x = seq_platform, y = mean)) +
  geom_bar(stat = "identity", aes(fill = Species, colour = Species)) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) +
  xlab("\nSequencing Platform") +
  ylab("Relative Abundance\n") +
  ggtitle("Endozoicomonas - P. lobata") 



vib.mfol <- subset_taxa(mfol.ra, Genus == "Vibrio")
vib.plob <- subset_taxa(plob.ra, Genus == "Vibrio")
vib.mfol.melt <- psmelt(vib.mfol)
vib.plob.melt <- psmelt(vib.plob)

sum.mfol.vib <- ddply(vib.mfol.melt, c("seq_platform", "Species"), summarise,
                      N = length(Abundance), 
                      mean = mean(Abundance),
                      sd = sd(Abundance), 
                      se = sd/sqrt(N)
)

ggplot(sum.mfol.vib, aes(x = seq_platform, y = mean)) +
  geom_bar(stat = "identity", aes(fill = Species, colour = Species)) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) +
  xlab("\nSequencing Platform") +
  ylab("Relative Abundance\n") +
  ggtitle("Vibrio - M.fol") 


sum.plob.vib <- ddply(vib.plob.melt, c("seq_platform", "Species"), summarise,
                      N = length(Abundance), 
                      mean = mean(Abundance),
                      sd = sd(Abundance), 
                      se = sd/sqrt(N)
)

ggplot(sum.plob.vib, aes(x = seq_platform, y = mean)) +
  geom_bar(stat = "identity", aes(fill = Species, colour = Species)) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) +
  xlab("\nSequencing Platform") +
  ylab("Relative Abundance\n") +
  ggtitle("Vibrio - P. lobata") 




##What happens if you remove endozoicomonas from the entire dataset? 
mfol.nr.ne <- subset_taxa(mfol.nr, Family != "Endozoicomonadaceae")
data.mfol.nr.ne <- as(sample_data(mfol.nr.ne), "data.frame")
mfol.ra.ne <- subset_taxa(mfol.ra, Family != "Endozoicomonadaceae")
data.mfol.ra.ne <- as(sample_data(mfol.ra.ne), "data.frame")
plob.nr.ne <- subset_taxa(plob.nr, Family != "Endozoicomonadaceae")
data.plob.nr.ne <- as(sample_data(plob.nr.ne), "data.frame")
plob.ra.ne <- subset_taxa(plob.ra, Family != "Endozoicomonadaceae")
data.plob.ra.ne <- as(sample_data(plob.ra.ne), "data.frame")

#Jaccard 
jc.mfol.ne <- phyloseq::distance(mfol.nr.ne, method = "jaccard", binary = TRUE)
adonis(jc.mfol.ne ~ seq_platform, strata = data.mfol.nr.ne$sample_label, data = data.mfol.nr.ne)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.5213 0.52132  1.1603 0.05483  0.002 **
#  Residuals    20    8.9864 0.44932         0.94517          
#Total        21    9.5077                 1.00000  
permutest(betadisper(jc.mfol.ne, data.mfol.nr.ne$seq_platform, type = "centroid"))
#Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.0003230 0.00032304 1.3946    999  0.276
#Residuals 20 0.0046328 0.00023164 

bc.mfol.ne <- phyloseq::distance(mfol.ra.ne, method = "bray")
adonis(bc.mfol.ne ~ seq_platform, strata = data.mfol.ra.ne$sample_label, data = data.mfol.ra.ne)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.4807 0.48067  1.1361 0.05375   0.01 **
#  Residuals    20    8.4613 0.42307         0.94625          
#Total        21    8.9420                 1.00000 
permutest(betadisper(bc.mfol.ne, data.mfol.ra.ne$seq_platform, type = "centroid"))
#Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.0000017 0.00000167 0.0017    999  0.966
#Residuals 20 0.0196978 0.00098489   


jc.plob.ne <- phyloseq::distance(plob.nr.ne, method = "jaccard", binary = TRUE)
adonis(jc.plob.ne ~ seq_platform, strata = data.plob.nr.ne$sample_label, data = data.plob.nr.ne)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.6331 0.63313  1.4504 0.05699  0.002 **
#  Residuals    24   10.4761 0.43651         0.94301          
#Total        25   11.1093                 1.00000          
#---
permutest(betadisper(jc.plob.ne, data.plob.nr.ne$seq_platform, type = "centroid"))
#        Df   Sum Sq    Mean Sq     F N.Perm Pr(>F)
#Groups     1 0.000493 0.00049305 0.401    999   0.51
#Residuals 24 0.029512 0.00122966   

bc.plob.ne <- phyloseq::distance(plob.ra.ne, method = "bray")
adonis(bc.plob.ne ~ seq_platform, strata = data.plob.ra.ne$sample_label, data = data.plob.ra.ne)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    0.9863 0.98628  2.7985 0.10443  0.001 ***
#  Residuals    24    8.4582 0.35243         0.89557           
#Total        25    9.4445                 1.00000           
#---
permutest(betadisper(bc.plob.ne, data.plob.ra.ne$seq_platform, type = "centroid"))
#Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)   
#Groups     1 0.22919 0.229192 9.6371    999  0.004 **
#  Residuals 24 0.57077 0.023782  
