###NMDS by sample for revised Supplementary Figures 2 & 3##

##Re-run the following code with sample as colour using the physeq.nr object
#from both ASV and OTU downstream R files. 
sample_data(physeq.nr) #Total samples:48 #44 for OTUs
#Count Sequences
totalreads <- sample_sums(physeq.nr)
totalreads$host <- data.nr$field_host_name
totalreads$seq_platform <- data.nr$seq_platform
View(totalreads)
library(microbiome)
reads <- readcount(physeq.nr)
sample_data(physeq.nr)$reads <- reads
data.nr.reads <- as(sample_data(physeq.nr), "data.frame")
sum.reads <- ddply(data.nr.reads, c("field_host_name", "seq_platform"), summarise,
                N = length(reads),
                mean = mean(reads),
                sd = sd(reads), 
                se = sd/sqrt(N),
                median = median(reads),
                min = min(reads),
                max = max(reads)
)


ord.mfol.jc <- ordinate(mfol.nr, "NMDS", "jaccard", binary = TRUE, trymax = 100) 
stressplot(ord.mfol.jc)
scores.mfol.jc <- as.data.frame(scores(ord.mfol.jc))
scores.mfol.jc$platform <- mfol.data.nr$seq_platform
scores.mfol.jc$sample <- mfol.data.nr$sample_label
ggplot(scores.mfol.jc, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = sample, shape = platform), size = 4) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, linetype = platform)) +
  ggtitle("M. aequituberculata Jaccard") +
  guides(fill = FALSE, color = FALSE, shape = FALSE, linetype = FALSE) +
  theme_classic()
ggsave("~/Desktop/new_figs/otu_jc_mfol_sample.pdf")

ord.plob.jc <- ordinate(plob.nr, "NMDS", "jaccard", binary = TRUE, trymax = 100) 
stressplot(ord.plob.jc)
scores.plob.jc <- as.data.frame(scores(ord.plob.jc))
scores.plob.jc$platform <- plob.data.nr$seq_platform
scores.plob.jc$sample <- plob.data.nr$sample_label
ggplot(scores.plob.jc, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = sample, shape = platform), size = 4) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, linetype = platform)) +
  ggtitle("P. lobata Jaccard") +
  #guides(fill = FALSE, color = FALSE, shape = FALSE, linetype = FALSE) +
  geom_text(label = rownames(scores.plob.jc), nudge_x = 0.05, nudge_y = 0.05,check_overlap = T) +
  theme_classic()
ggsave("~/Desktop/new_figs/legend.pdf")

ord.mfol <- ordinate(mfol.ra, "NMDS", "bray", trymax = 100) 
stressplot(ord.mfol)
scores.mfol <- as.data.frame(scores(ord.mfol))
scores.mfol$platform <- mfol.data.ra$seq_platform
scores.mfol$sample <- mfol.data.ra$sample_label
ggplot(scores.mfol, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = sample, shape = platform), size = 4) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, linetype = platform)) +
  ggtitle("M. aequituberculata Bray Curtis") +
  guides(fill = FALSE, color = FALSE, shape = FALSE, linetype = FALSE) +
  theme_classic()
ggsave("~/Desktop/new_figs/otu_br_mfol_sample.pdf")

ord.plob <- ordinate(plob.ra, "NMDS", "bray", trymax = 500) 
stressplot(ord.plob)
scores.plob <- as.data.frame(scores(ord.plob))
scores.plob$platform <- plob.data.ra$seq_platform
scores.plob$sample <- plob.data.ra$sample_label
ggplot(scores.plob, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = sample, shape = platform), size = 4) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, linetype = platform)) +
  ggtitle("P. lobata Bray Curtis") +
  guides(fill = FALSE, color = FALSE, shape = FALSE, linetype = FALSE) +
  theme_classic()
ggsave("~/Desktop/new_figs/otu_br_plob_sample.pdf")

ord.un.mfol <- ordinate(mfol.ra, "NMDS", "uunifrac") #trymax is not accepted here
stressplot(ord.un.mfol)
scores.un.mfol <- as.data.frame(scores(ord.un.mfol))
scores.un.mfol$platform <- mfol.data.ra$seq_platform
scores.un.mfol$sample <- mfol.data.ra$sample_label
ggplot(scores.un.mfol, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = sample, shape = platform), size = 4) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, linetype = platform)) +
  ggtitle("M. aequituberculata Unweighted UniFrac") +
  guides(fill = FALSE, color = FALSE, shape = FALSE, linetype = FALSE) +
  theme_classic()
ggsave("~/Desktop/new_figs/otu_un_mfol_sample.pdf")

ord.un.plob <- ordinate(plob.ra, "NMDS", "uunifrac") 
stressplot(ord.un.plob)
scores.un.plob <- as.data.frame(scores(ord.un.plob))
scores.un.plob$platform <- plob.data.ra$seq_platform
scores.un.plob$sample <- plob.data.ra$sample_label
ggplot(scores.un.plob, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = sample, shape = platform), size = 4) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, linetype = platform)) +
  ggtitle("P. lobata Unweighted UniFrac") +
  #guides(fill = FALSE, color = FALSE, shape = FALSE, linetype = FALSE) +
  theme_classic()
ggsave("~/Desktop/new_figs/otu_un_plob_sample.pdf")


ord.wun.mfol <- ordinate(mfol.ra, "NMDS", "wunifrac") #trymax is not accepted here
stressplot(ord.wun.mfol)
scores.wun.mfol <- as.data.frame(scores(ord.wun.mfol))
scores.wun.mfol$platform <- mfol.data.ra$seq_platform
scores.wun.mfol$sample <- mfol.data.ra$sample_label
ggplot(scores.wun.mfol, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = sample, shape = platform), size = 4) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, linetype = platform)) +
  ggtitle("M. aequituberculata Weighted UniFrac") +
  #guides(fill = FALSE, color = FALSE, shape = FALSE, linetype = FALSE) +
  theme_classic()
ggsave("~/Desktop/new_figs/legend_plob.pdf")

#Porites lobata
ord.wun.plob <- ordinate(plob.ra, "NMDS", "wunifrac") 
stressplot(ord.wun.plob)
scores.wun.plob <- as.data.frame(scores(ord.wun.plob))
scores.wun.plob$platform <- plob.data.ra$seq_platform
scores.wun.plob$sample <- plob.data.ra$sample_label
ggplot(scores.wun.plob, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = sample, shape = platform), size = 4) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, linetype = platform)) +
  ggtitle("P. lobata Weighted UniFrac") +
  #guides(fill = FALSE, color = FALSE, shape = FALSE, linetype = FALSE) +
  theme_classic()
ggsave("~/Desktop/new_figs/otu_wun_plob_sample.pdf")
