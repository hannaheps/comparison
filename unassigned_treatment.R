###Removing extra few Unassigned Eukaryotic Reads###
#This script is to check if removing the small subset of unassigned bacterial reads that blasted to eukaryotes
#changes the outcome of the paper

##Lots of unassigned reads in the ASV dataset - check by blast (same as mitos above) & remove any that have eukaryotic hits
#This is not the same in OTU dataset.
physeq.unassigned <- physeq.nr
tax.clean <- data.frame(tax_table(physeq.unassigned))
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- "Unassigned"
tax_table(physeq.unassigned) <- as.matrix(tax.clean)
tax_table(physeq.unassigned)
physeq.unassigned <- subset_taxa(physeq.unassigned, Phylum == "Unassigned")
tax_table(physeq.unassigned) #1331 taxa are unassigned
unassigned <- rownames(otu_table(physeq.unassigned))
write.csv(unassigned, "unassigned_ids.csv")
#With this csv, run blast
#Remove all eukaryotic taxa
euk.taxa <- read.csv("euk_asvs.csv", header = TRUE)
euk.taxa <- euk.taxa$feature.id
#remove them from the dataset
all.taxa <- taxa_names(physeq.nr)
all.taxa <- all.taxa[!(all.taxa %in% euk.taxa)]
physeq.nr <- prune_taxa(all.taxa, physeq.nr)

physeq.r <- rarefy_even_depth(physeq.prune, sample.size = 1000, rngseed = 711) 
#Remove bad samples
bad.samples.r <- c("KI15BFMD107", "KI15BFMD108", "KI15BFMD119", "KI15BFMD123", "KI15BFMD135", "KI15BFMD145", "KI15BFMD195", "KI15BFMD196", "KI15BFMD208", "KI15BFMD232_EMP")
all.samples.r <- sample_names(physeq.r)
all.samples.r <- all.samples.r[!(all.samples.r %in% bad.samples.r)]
physeq.r <- prune_samples(all.samples.r, physeq.r)
#remove eukaryotic unassigned reads
all.taxa.r <- taxa_names(physeq.r)
all.taxa.r <- all.taxa.r[!(all.taxa.r %in% euk.taxa)]
physeq.r <- prune_taxa(all.taxa.r, physeq.r)
saveRDS(physeq.r, "physeq_r.RDS")
physeq.r <- readRDS("physeq_r.RDS")
data.r <- as(sample_data(physeq.r), "data.frame")
#Relative abundance
physeq.ra <- transform_sample_counts(physeq.nr, function(x) x/ sum(x)) 
data.ra <- as(sample_data(physeq.ra), "data.frame")

#Relative abundance with rarefied data
physeq.r.ra <- transform_sample_counts(physeq.r, function(x) x/ sum(x)) 
data.r.ra <- as(sample_data(physeq.r), "data.frame")


###barplot###
percent.trial <- physeq.nr %>% 
  tax_glom(taxrank = "Phylum", NArm=FALSE) %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) 
perc.melt <- psmelt(percent.trial)

perc.melt$sample_label <- factor(perc.melt$sample_label, levels = c("KI15BFMD082", "KI15BFMD089", "KI15BFMD124", "KI15BFMD142",
                                                                    "KI15BFMD144", "KI15BFMD146", "KI15BFMD190", "KI15BFMD193",
                                                                    "KI15BFMD221", "KI15BFMD228", "KI15BFMD240", "KI15BFMD079",
                                                                    "KI15BFMD086", "KI15BFMD091", "KI15BFMD105", "KI15BFMD112", 
                                                                    "KI15BFMD115", "KI15BFMD122", "KI15BFMD141", "KI15BFMD189",
                                                                    "KI15BFMD201", "KI15BFMD229", "KI15BFMD230", "KI15BFMD235"))

sum.ra <- ddply(perc.melt, c("Phylum", "seq_platform", "sample_label"), summarise,
                N = length(Abundance),
                mean = mean(Abundance),
                sd = sd(Abundance), 
                se = sd/sqrt(N)
)

sum.ra$Phylum = factor(as.character(sum.ra$Phylum), levels=names(x))
sum.ra$col_phylum <- sum.ra$Phylum
#View(sum.ra$col_phylum)
#Set everything that is super low to NA so that we can call them "other"
sum.ra$col_phylum <- as.character(sum.ra$col_phylum)
sum.ra$col_phylum <- ifelse(is.na(sum.ra$col_phylum), 
                            'Unassigned', sum.ra$col_phylum)
sum.ra$col_phylum <- as.factor(sum.ra$col_phylum)

sum.ra$col_phylum[sum.ra$col_phylum != "Proteobacteria" &
                    sum.ra$col_phylum != "Firmicutes" & 
                    sum.ra$col_phylum != "Bacteroidetes" &
                    sum.ra$col_phylum != "Cyanobacteria" &
                    sum.ra$col_phylum != "Verrucomicrobia" &
                    sum.ra$col_phylum != "Euryarchaeota" &
                    sum.ra$col_phylum != "Actinobacteria" &
                    sum.ra$col_phylum != "Epsilonbacteraeota" &
                    sum.ra$col_phylum != "Planctomycetes" &
                    sum.ra$col_phylum != "Deinococcus-Thermus" &
                    sum.ra$col_phylum != "Unassigned"] <- NA


levels(sum.ra$col_phylum)
# add new factor
sum.ra$col_phylum <- factor(sum.ra$col_phylum, levels = c(levels(sum.ra$col_phylum), "Other"))
# convert NAs to other
sum.ra$col_phylum [is.na(sum.ra$col_phylum)] = "Other"
sum.ra$col_phylum <- factor(x = sum.ra$col_phylum, levels = c("Proteobacteria", "Firmicutes", "Bacteroidetes",
                                                              "Cyanobacteria", "Verrucomicrobia", "Euryarchaeota",
                                                              "Actinobacteria", "Epsilonbacteraeota", "Planctomycetes",
                                                              "Deinococcus-Thermus", "Other", "Unassigned"))
#Make a colour scheme
ale.colors <-  c( "#80B1D3","#FFFFB3","#ABA3E6","#FB8072","#8DD3C7", "#FDB462", "#B3DE69",
                  "#FCCDE5", "#BC80BD", "#FFED6F", "#CCEBC5", "#A3A5A8")

#plot!
p <- ggplot(sum.ra, aes(x = sample_label, y = mean, fill = col_phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=ale.colors) +
  ylab("Relative Abundance") +
  xlab("Sample Label") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#To re-order the sample names so it goes M. aquituberculata first, then P. lobata
p$data$sample_label <- factor(x = p$data$sample_label, levels = c("KI15BFMD082", "KI15BFMD089", "KI15BFMD124", "KI15BFMD142",
                                                                  "KI15BFMD144", "KI15BFMD146", "KI15BFMD190", "KI15BFMD193",
                                                                  "KI15BFMD221", "KI15BFMD228", "KI15BFMD240", "KI15BFMD079",
                                                                  "KI15BFMD086", "KI15BFMD091", "KI15BFMD105", "KI15BFMD112", 
                                                                  "KI15BFMD115", "KI15BFMD122", "KI15BFMD141", "KI15BFMD189",
                                                                  "KI15BFMD201", "KI15BFMD229", "KI15BFMD230", "KI15BFMD235"))
#facet grid by platform
p + facet_grid(rows = vars(seq_platform)) 

ggsave("barplot_asv_noeuk.pdf", plot = last_plot())

###Betadiversity###

##This is run twice - once for unrarefied and once for rarefied data
#Split non unrarefied data by species
mfol.nr <- subset_samples(physeq.nr, field_host_name == "Montipora foliosa")
mfol.data.nr <- as(sample_data(mfol.nr), "data.frame")
mfol.ra <- transform_sample_counts(mfol.nr, function(x) x/ sum(x)) 
mfol.data.ra <- as(sample_data(mfol.ra), "data.frame")

plob.nr <- subset_samples(physeq.nr, field_host_name == "Porites lobata")
plob.data.nr <- as(sample_data(plob.nr), "data.frame")
plob.ra <- transform_sample_counts(plob.nr, function(x) x/ sum(x)) 
plob.data.ra <- as(sample_data(plob.ra), "data.frame")

#Split rarefied data by species
mfol.r <- subset_samples(physeq.r, field_host_name == "Montipora foliosa")
mfol.data.r <- as(sample_data(mfol.r), "data.frame")
mfol.r.ra <- transform_sample_counts(mfol.r, function(x) x/ sum(x)) 
mfol.data.r.ra <- as(sample_data(mfol.r.ra), "data.frame")

plob.r <- subset_samples(physeq.r, field_host_name == "Porites lobata")
plob.data.r <- as(sample_data(plob.r), "data.frame")
plob.r.ra <- transform_sample_counts(plob.r, function(x) x/ sum(x)) 
plob.data.r.ra <- as(sample_data(plob.r.ra), "data.frame")

#Unrarefied
## Binary Jaccard
#Montipora foliosa
jc.mfol <- phyloseq::distance(mfol.nr, method = "jaccard", binary = TRUE)
adonis(jc.mfol ~ seq_platform, strata = mfol.data.nr$sample_label, data = mfol.data.nr)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.5365 0.53655  1.1688 0.05521  0.002 **
#  Residuals    20    9.1813 0.45906         0.94479          
#Total        21    9.7178                 1.00000 
permutest(betadisper(jc.mfol, mfol.data.nr$seq_platform, type = "centroid"))
#Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.0000214 2.1446e-05 0.1118    999  0.749
#Residuals 20 0.0038373 1.9187e-04     

#weighted Unifrac (presence absence + abundance + phylogeny)
wu.mfol <- phyloseq::distance(mfol.ra, method = "wunifrac")
adonis(wu.mfol ~ seq_platform, strata = mfol.data.ra$sample_label, data = mfol.data.ra)
#Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)
#seq_platform  1  0.003914 0.0039141 0.83574 0.04011  0.291
#Residuals    20  0.093669 0.0046835         0.95989       
#Total        21  0.097583                   1.00000  
permutest(betadisper(wu.mfol, mfol.data.ra$seq_platform, type = "centroid"))
# Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.0002975 0.00029752 0.2882    999  0.605
#Residuals 20 0.0206442 0.00103221        


#Porites lobata
jc.plob <- phyloseq::distance(plob.nr, method = "jaccard", binary = TRUE)
adonis(jc.plob ~ seq_platform, strata = plob.data.nr$sample_label, data = plob.data.nr)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1     0.623 0.62302  1.4056 0.05533  0.001 ***
#  Residuals    24    10.637 0.44323         0.94467           
#Total        25    11.261                 1.00000            
permutest(betadisper(jc.plob, plob.data.nr$seq_platform, type = "centroid"))
#Response: Distances
#Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.0009018 0.00090183 0.8685    999  0.349
#Residuals 24 0.0249216 0.00103840  

##Weighted Unifrac
wu.plob <- phyloseq::distance(plob.ra, method = "wunifrac")
adonis(wu.plob ~ seq_platform, strata = plob.data.ra$sample_label, data = plob.data.ra)
#  Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
#seq_platform  1  0.010051 0.010051 0.80771 0.03256  0.012 *
#  Residuals    24  0.298662 0.012444         0.96744         
#Total        25  0.308713                  1.00000         
#---
permutest(betadisper(wu.plob, plob.data.ra$seq_platform, type = "centroid"))
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.005762 0.0057619 1.1725    999  0.339
#Residuals 24 0.117941 0.0049142    


##Bray curtis - uses relative abundance to account for differences in sequencing depth
#Montipora foliosa
bc.mfol <- phyloseq::distance(mfol.ra, method = "bray")
adonis(bc.mfol ~ seq_platform, strata = mfol.data.ra$sample_label, data = mfol.data.ra)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.4903 0.49031  1.1409 0.05397  0.009 **
#  Residuals    20    8.5953 0.42976         0.94603          
#Total        21    9.0856                 1.00000
permutest(betadisper(bc.mfol, mfol.data.ra$seq_platform, type = "centroid"))
#Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.0000013 0.00000135 0.0015    999  0.968
#Residuals 20 0.0183962 0.00091981   

#Unweighted unifrac
un.mfol <- phyloseq::distance(mfol.nr, method = "uunifrac")
adonis(un.mfol ~ seq_platform, strata = mfol.data.nr$sample_label, data = mfol.data.nr)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#seq_platform  1    0.7425 0.74250  2.1445 0.09684  0.012 *
##  Residuals    20    6.9248 0.34624         0.90316         
#Total        21    7.6673                 1.00000         
        

permutest(betadisper(un.mfol, mfol.data.nr$seq_platform, type = "centroid"))
#Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)   
#Groups     1 0.061131 0.061131 10.741    999  0.007 **
#  Residuals 20 0.113826 0.005691     

#Porites lobata
bc.plob <- phyloseq::distance(plob.ra, method = "bray")
adonis(bc.plob ~ seq_platform, strata = plob.data.ra$sample_label, data = plob.data.ra)
#f SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    0.8554 0.85535  2.4607 0.09299  0.001 ***
#  Residuals    24    8.3426 0.34761         0.90701           
#Total        25    9.1979                 1.00000           
#---
permutest(betadisper(bc.plob, plob.data.ra$seq_platform, type = "centroid"))
#Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)   
#Groups     1 0.18983 0.189835 8.0732    999   0.01 **
#  Residuals 24 0.56434 0.023514                        


un.plob <- phyloseq::distance(plob.nr, method = "uunifrac")
adonis(un.plob ~ seq_platform, strata = plob.data.nr$sample_label, data = plob.data.nr)
#        Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.5398 0.53982   1.552 0.06074  0.002 **
#  Residuals    24    8.3478 0.34782         0.93926          
#Total        25    8.8876                 1.00000           
permutest(betadisper(un.plob, plob.data.nr$seq_platform, type = "centroid"))
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.004713 0.0047130 0.6482    999  0.421
#Residuals 24 0.174510 0.0072713  



#for ease of re-running
mfol.nr <- subset_samples(physeq.r, field_host_name == "Montipora foliosa")
mfol.data.nr <- as(sample_data(mfol.r), "data.frame")
mfol.ra <- transform_sample_counts(mfol.r, function(x) x/ sum(x)) 
mfol.data.ra <- as(sample_data(mfol.r.ra), "data.frame")

plob.nr <- subset_samples(physeq.r, field_host_name == "Porites lobata")
plob.data.nr <- as(sample_data(plob.r), "data.frame")
plob.ra <- transform_sample_counts(plob.r, function(x) x/ sum(x)) 
plob.data.ra <- as(sample_data(plob.r.ra), "data.frame")


#Rarefied
## Binary Jaccard
#Montipora foliosa
jc.mfol <- phyloseq::distance(mfol.nr, method = "jaccard", binary = TRUE)
adonis(jc.mfol ~ seq_platform, strata = mfol.data.nr$sample_label, data = mfol.data.nr)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.5260 0.52603  1.1554 0.05461  0.002 **
#  Residuals    20    9.1058 0.45529         0.94539          
#Total        21    9.6318                 1.00000   
permutest(betadisper(jc.mfol, mfol.data.nr$seq_platform, type = "centroid"))
# Df     Sum Sq    Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.00010184 0.00010184 0.6812    999  0.468
#Residuals 20 0.00298979 0.00014949         

#weighted Unifrac (presence absence + abundance + phylogeny)
wu.mfol <- phyloseq::distance(mfol.ra, method = "wunifrac")
adonis(wu.mfol ~ seq_platform, strata = mfol.data.ra$sample_label, data = mfol.data.ra)
#Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)
#seq_platform  1  0.003266 0.0032660 0.51374 0.02504  0.503
#Residuals    20  0.127146 0.0063573         0.97496       
#Total        21  0.130412                   1.00000   
permutest(betadisper(wu.mfol, mfol.data.ra$seq_platform, type = "centroid"))
# Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.0004054 0.00040535 0.3613    999  0.567
#Residuals 20 0.0224357 0.00112178       


#Porites lobata
jc.plob <- phyloseq::distance(plob.nr, method = "jaccard", binary = TRUE)
adonis(jc.plob ~ seq_platform, strata = plob.data.nr$sample_label, data = plob.data.nr)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    0.6167 0.61669  1.4168 0.05574  0.001 ***
#  Residuals    24   10.4467 0.43528         0.94426           
#Total        25   11.0634                 1.00000               
permutest(betadisper(jc.plob, plob.data.nr$seq_platform, type = "centroid"))
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.003007 0.0030073 1.8427    999  0.174
#Residuals 24 0.039168 0.0016320   

##Weighted Unifrac
wu.plob <- phyloseq::distance(plob.ra, method = "wunifrac")
adonis(wu.plob ~ seq_platform, strata = plob.data.ra$sample_label, data = plob.data.ra)
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
#seq_platform  1   0.01300 0.013003  0.9624 0.03855  0.011 *
#  Residuals    24   0.32425 0.013510         0.96145         
#Total        25   0.33725                  1.00000          
#---
permutest(betadisper(wu.plob, plob.data.ra$seq_platform, type = "centroid"))
#Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.011024 0.011024 2.6148    999   0.14
#Residuals 24 0.101183 0.004216  

##Bray curtis - uses relative abundance to account for differences in sequencing depth
#Montipora foliosa
bc.mfol <- phyloseq::distance(mfol.ra, method = "bray")
adonis(bc.mfol ~ seq_platform, strata = mfol.data.ra$sample_label, data = mfol.data.ra)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#seq_platform  1    0.4954 0.49541  1.1502 0.05438  0.007 **
##  Residuals    20    8.6142 0.43071         0.94562          
#Total        21    9.1096                 1.00000       
permutest(betadisper(bc.mfol, mfol.data.ra$seq_platform, type = "centroid"))
#Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.0000155 0.00001545 0.0158    999   0.89
#Residuals 20 0.0196026 0.00098013    

#Unweighted unifrac
un.mfol <- phyloseq::distance(mfol.nr, method = "uunifrac")
adonis(un.mfol ~ seq_platform, strata = mfol.data.nr$sample_label, data = mfol.data.nr)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
##seq_platform  1    0.4139 0.41386  1.1737 0.05543  0.082 .
#Residuals    20    7.0519 0.35260         0.94457         
#Total        21    7.4658                 1.00000    

permutest(betadisper(un.mfol, mfol.data.nr$seq_platform, type = "centroid"))
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.006501 0.0065015 1.5862    999  0.205
#Residuals 20 0.081975 0.0040987       

#Porites lobata
bc.plob <- phyloseq::distance(plob.ra, method = "bray")
adonis(bc.plob ~ seq_platform, strata = plob.data.ra$sample_label, data = plob.data.ra)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#seq_platform  1    0.8465 0.84646    2.43 0.09194  0.001 ***
#  Residuals    24    8.3601 0.34834         0.90806           
#Total        25    9.2066                 1.00000  
permutest(betadisper(bc.plob, plob.data.ra$seq_platform, type = "centroid"))
#Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)   
#Groups     1 0.19614 0.196137 8.2881    999  0.008 **
#  Residuals 24 0.56796 0.023665                        


un.plob <- phyloseq::distance(plob.nr, method = "uunifrac")
adonis(un.plob ~ seq_platform, strata = plob.data.nr$sample_label, data = plob.data.nr)
#        Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#seq_platform  1    0.4583 0.45833  1.3582 0.05356  0.092 .
#Residuals    24    8.0990 0.33746         0.94644         
#Total        25    8.5574                 1.00000             
permutest(betadisper(un.plob, plob.data.nr$seq_platform, type = "centroid"))
#Df   Sum Sq   Mean Sq     F N.Perm Pr(>F)
#Groups     1 0.000001 0.0000008 1e-04    999  0.994
#Residuals 24 0.214988 0.0089578    




