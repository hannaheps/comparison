###Plotting functions for both ASVs and OTUs for relative abundance bar plots
##To change colour scheme according to paper revisions

##1. ASVs
#Use the "sum.ra" object from the asv_downstream.R file
#First, order the phyla and name only significant phyla (all else as OTHER)
x = tapply(sum.ra$mean, sum.ra$Phylum, function(x) max(x))
x = sort(x, TRUE) #sort phyla 
View(x) #add top 10 to code below to set NAs
sum.ra$Phylum = factor(as.character(sum.ra$Phylum), levels=names(x))
sum.ra$col_phylum <- sum.ra$Phylum
View(sum.ra$col_phylum)
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

ggsave("barplot_asvs_revised.pdf", plot = last_plot())



##2. Plotting the OTUs
##Same deal - take sum.ra object from otu_downstream.R file
#Basically same code as above, but you have already set your colours so no need to do this again:

#First, order the phyla and name only significant phyla (all else as OTHER)
x = tapply(sum.ra$mean, sum.ra$Phylum, function(x) max(x))
x = sort(x, TRUE) #sort phyla 
View(x) #add top 10 to code below to set NAs
sum.ra$Phylum = factor(as.character(sum.ra$Phylum), levels=names(x))
sum.ra$col_phylum <- sum.ra$Phylum
View(sum.ra$col_phylum)
#Set everything that is super low to NA so that we can call them "other"
sum.ra$col_phylum <- as.character(sum.ra$col_phylum)
sum.ra$col_phylum <- ifelse(is.na(sum.ra$col_phylum), 
                            'Unassigned', sum.ra$col_phylum) #first change NA to "Unassigned"
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
                    sum.ra$col_phylum != "Unassigned" ] <- NA


levels(sum.ra$col_phylum)
# add new factor
sum.ra$col_phylum <- factor(sum.ra$col_phylum, levels = c(levels(sum.ra$col_phylum), "Other"))
# convert NAs to other
sum.ra$col_phylum [is.na(sum.ra$col_phylum)] = "Other"
sum.ra$col_phylum <- factor(x = sum.ra$col_phylum, levels = c("Proteobacteria", "Firmicutes", "Bacteroidetes",
                                                              "Cyanobacteria", "Verrucomicrobia", "Euryarchaeota",
                                                              "Actinobacteria", "Epsilonbacteraeota", "Planctomycetes",
                                                              "Deinococcus-Thermus", "Other", "Unassigned"))

levels(sum.ra$col_phylum)
# add new factor
sum.ra$col_phylum <- factor(sum.ra$col_phylum, levels = c(levels(sum.ra$col_phylum), "Other"))
# convert NAs to other
sum.ra$col_phylum [is.na(sum.ra$col_phylum)] = "Other"

p1 <- ggplot(sum.ra, aes(x = sample_label, y = mean, fill = col_phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=ale.colors) +
  ylab("Relative Abundance") +
  xlab("Sample Label") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#To re-order the sample names so it goes M. aquituberculata first, then P. lobata
p1$data$sample_label <- factor(x = p1$data$sample_label, levels = c("KI15BFMD082", "KI15BFMD089", "KI15BFMD124", "KI15BFMD142",
                                                                  "KI15BFMD144", "KI15BFMD146", "KI15BFMD190", "KI15BFMD193",
                                                                  "KI15BFMD221", "KI15BFMD228", "KI15BFMD240", "KI15BFMD079",
                                                                  "KI15BFMD086", "KI15BFMD091", "KI15BFMD105", "KI15BFMD112", 
                                                                  "KI15BFMD115", "KI15BFMD122", "KI15BFMD141", "KI15BFMD189",
                                                                  "KI15BFMD201", "KI15BFMD229", "KI15BFMD230", "KI15BFMD235"))
#facet grid by platform
p1 + facet_grid(rows = vars(seq_platform)) 

ggsave("barplot_otus_revised.pdf", plot = last_plot())
