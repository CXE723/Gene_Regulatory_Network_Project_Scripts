### Gene-gene Co-correlation networks for Cytoscape ###

### Packages needed.
library(dplyr)
library(corrr)
library(psych)
library(pheatmap)
library(corrplot)
library(ggplot2)
library(ggcorrplot)
library(tidyr)

###------------------------------------------------------------------------------
### Load in matrix for gene-motif counts.
CD34 <- CD34_motifs_counts_3
CD34 <- read.csv("CD34_motifs_counts_3.csv", sep = ",")

# Formatting so suitable for correlation analysis. 
CD34 <- CD34 %>% mutate_if(is.character, as.factor)
CD34 <- CD34 %>% mutate_if(is.factor, as.numeric)

# Correlation analysis of Genes in CD34+.
cor_CD34 <- cor(CD34, y = NULL, use = "pairwise.complete.obs", method = "pearson")

# Convert to data frame and remove NA values. 
cor_CD34 <- as.data.frame(as.table(cor_CD34))
cor_CD34 <- subset(cor_CD34, Var1!='?..gene' & Var2!='?..gene')
cor_CD34 <- na.omit(cor_CD34)

# Subset for correlation frequency greater than 0.5. 
cor_CD34_2 <- subset(cor_CD34, abs(Freq) > 0.5)
cor_CD34_2 <- cor_CD34_2[order(-abs(cor_CD34_2$Freq)),]

#Save
write.csv(cor_CD34, file = "CD34_correlation_on_motifCont.csv", sep = ',')
write.csv(cor_CD34_2, file = "CD34_correlation_on_motifCont_05.csv", sep = ',')

#Plot heatmap
pheatmap(cor_CD34_2, border_color = NA, color = colorRampPalette(c('white', 'red', "blue"))(256), fontsize = 18, lwd = 2)

###------------------------------------------------------------------------------
### Repeated for t(3;21).

t321 <- read.csv("t321_motifs_cor.csv", sep = ",")

colnames(t321)

t321 <- t321 %>% mutate_if(is.character, as.factor)
t321 <- t321 %>% mutate_if(is.factor, as.numeric)

cor_t321 <- cor(t321, y = NULL, use = "pairwise.complete.obs", method = "pearson")

cor_t321 <- as.data.frame(as.table(cor_t321))
cor_t321 <- subset(cor_t321, Var1!='?..Gene.Name' & Var2!='?..Gene.Name')
cor_t321 <- na.omit(cor_t321)

cor_t321_2 <- subset(t321, abs(Freq) > 0.5)
cor_t321_2 <- cor_t321_2[order(-abs(cor_t321_2$Freq)),]

write.csv(cor_t321, file = "t321_correlation_on_motifCont.csv", sep = ',')
write.csv(cor_t321_2, file = "t321_correlation_on_motifCont_05.csv", sep = ',')

###------------------------------------------------------------------------------
### Testing alternate heatmaps.

CD34_correlation <- read.delim("CD34_correlation_on_motifCont.txt", sep = "\t")
t321_correlation <- read.delim("t321_correlation_on_motifCont.txt", sep = "\t")

# CD34+
CD34_cor_map <- pivot_longer(data = CD34_correlation,
                             cols = -c(1),
                             names_to = "Gene",
                             values_to = "Frequency")
head(CD34_cor_map)

CD34.heatmap <- ggplot(data = CD34_cor_map, mapping = aes(x = X,
                                                          y = Gene,
                                                          fill= Frequency)) +
  geom_tile() +
  xlab(label = "Gene") +
  theme(
    axis.text.x = element_text(angle = 90, size = 6),
    axis.text.y = element_text(size = 6))

     
CD34.heatmap

# t(3;21)
t321_cor_map <- pivot_longer(data = t321_correlation,
                             cols = -c(1),
                             names_to = "Gene",
                             values_to = "Frequency")
head(t321_cor_map)

t321.heatmap <- ggplot(data = t321_cor_map, mapping = aes(x = X,
                                                          y = Gene,
                                                          fill= Frequency)) +
  geom_tile() +
  xlab(label = "Gene") +
  theme(
    axis.text.x = element_text(angle = 90, size = 2),
    axis.text.y = element_text(size = 2))

t321.heatmap

###------------------------------------------------------------------------------