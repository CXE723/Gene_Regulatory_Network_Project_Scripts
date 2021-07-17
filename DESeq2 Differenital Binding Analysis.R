### DESeq2 Analysis for Differential Binding/ Extracting PeaK Lists ###

### Packages needed. 
library(ggplot2)
library(pheatmap)
library(lDESeq2)
library(limma)
###
###------------------------------------------------------------------------------
# DESeq Analysis 

# Load in annotated peak list and round values. 
rawCounts <- read.csv('Annotated_Peaks_Counts.csv', header = TRUE, comment.char = '#', row.names = 1)
rawCounts <- round(rawCounts, digits = 0)
head(rawCounts)

# Create a dataframe of the column data.
colnames(rawCounts)<- c('CD34_Rep1', 'CD34_Rep2', 't321_Rep1', 't321_Rep2')
head(rawCounts)

colData<- data.frame(sampleID = colnames(rawCounts), 
                     condition = factor(c('CD34', 'CD34', 't321', 't321')), 
                     replicate = factor(c('Rep1', 'Rep2', 'Rep1', 'Rep2')))
colData

# Format data so that is suitable for DESeq analysis and run DESeq. 
dds<- DESeqDataSetFromMatrix(countData = rawCounts, colData = colData, design = ~ condition)

dds<- DESeq(dds)

# Quality Control by clustering from DESeq2 Data and plot as a heatmap. 
logCounts<- assay(rlog(dds))
pCor<- cor(logCounts, method = 'pearson')
pheatmap(pCor, border_color = NA, color = colorRampPalette(c('white', 'red'))(256), fontsize = 18, lwd = 2)


# Principal Component analysis of the log counts.
pca<- prcomp(t(logCounts))

# Combining the results from the PCA and column annotation data.
pca_results<- cbind(as.data.frame(pca$x), colData)
pca_results

# Calculating the proportion of the variance contained in the PCs.
propVar<- pca$sdev ^ 2 / sum(pca$sdev ^ 2) * 100
propVar

# Plot a PCA with the proportion of variance added for PC1 and PC2.
g<- ggplot(pca_results, aes(x = PC1, y = PC2, shape = replicate, color = condition)) + geom_point(size = 6) +
  xlab(paste('PC1 (', round(propVar[1], 2), '%)', sep = '')) + # 
  ylab(paste('PC2 (', round(propVar[2], 2), '%)', sep = '')) + # 
  theme(axis.title = element_text(size = 24), axis.text = element_text(size = 24))
g

### -----------------------------------------------------------------------------
# Extracting Significant Peaks

# Results from DESEq2.
resTable<- results(dds, contrast = c('condition', 'CD34', 't321'))

# Convert to a data frame and order by log2fold-change.
resTable<- as.data.frame(resTable)
resTable<- resTable[order(resTable$log2FoldChange, decreasing = TRUE),]

head(resTable)

# Get CD34/t321 specific sites from the results and paste the number of specific peaks found for t(3;21) and CD34+.
# Significant peaks identified by log2foldchange and adjusted p-values.
resTable.CD34_specific<- resTable[resTable$log2FoldChange > 1 & resTable$padj < 0.1,]
head(resTable.CD34_specific)
paste('Found', nrow(resTable.CD34_specific), 'CD34 specific peaks', sep = ' ')

resTable.t321_specific<- resTable[resTable$log2FoldChange < -1 & resTable$padj < 0.1,]
head(resTable.t321_specific)
paste('Found', nrow(resTable.t321_specific), 't321 specific peaks', sep = ' ')

# Plot a Volcano Plot.

head(resTable)

# Volcano plot using log2fold-change versus the adjusted p-values.
plot(x = resTable$log2FoldChange, y = -log2(resTable$padj), pch = 16, col = 'grey60', xlab = 'log2 fold-change', ylab = '-log2(adj. p-value)', cex.axis = 1.5, cex.lab = 1.5)

# Add the points for the t(3;21) and CD34+ specific results.
points(x = resTable.CD34_specific$log2FoldChange, y = -log2(resTable.CD34_specific$padj), pch = 16, col = 'red')
points(x = resTable.t321_specific$log2FoldChange, y = -log2(resTable.t321_specific$padj), pch = 16, col = 'blue')

# Lines which show the fold-change and p-value cutoffs. 
abline(v = c(-1,1), lwd = 2)
abline(h = -log2(0.1), lwd = 2)


# Extract only the counts for the differential peaks from the logCounts.
diffPeak_counts<- logCounts[row.names(logCounts) %in% c(row.names(resTable.CD34_specific), row.names(resTable.t321_specific)),]

# Transform the counts to the Z-Scale
zCounts<- t(scale(t(diffPeak_counts)))

# Plot in heatmap to see clustering.
pheatmap(zCounts, color = colorRampPalette(c('blue', 'white', 'red'))(256), show_rownames = FALSE, lwd = 2, fontsize = 18)

#------------------------------------------------------------------------------------------
### Adding strand information back to specific peak lists for HOMER.

# Checking data match in full peaks vs specific lists. 
# List of t321 specific peaks
t321_peaks <- read.csv('t321_specific.csv')
colnames(t321_peaks)
t321_peaks <- t321_peaks[c(1)] 
t321_peaks <- t321_peaks %>% 
  rename(
    peak_id = X
  )
head(t321_peaks)

full_peaks <- full_peaks[c(1,5)]
full_peaks <- full_peaks %>%
  rename( 
    peak_id = ?..
  )
head(full_peaks)

# Check all t321 in list of full peaks list. 
all(t321_peaks$peak_id %in% full_peaks$peak_id)

# Subset data including the strand information. 
strand_data_t321 <- subset(full_peaks, peak_id %in% t321_peaks$peak_id)
head(strand_data_t321)

write.csv(as.data.frame(strand_data_t321),file="t321_strand_data.csv")

### Repeat for CD34+ list.
# List of CD34 specific peaks
CD34_peaks <- read.csv('CD34_specific.csv')
head(CD34_peaks)
colnames(CD34_peaks)

CD34_peaks <- CD34_peaks[c(1)] 
CD34_peaks <- CD34_peaks %>% 
  rename(
    peak_id = X
  )
head(CD34_peaks)

# Check all CD34+ in full peaks list.
all(CD34_peaks$peak_id %in% full_peaks$peak_id)

strand_data_CD34 <- subset(full_peaks, peak_id %in% CD34_peaks$peak_id)
head(strand_data_CD34)

write.csv(as.data.frame(strand_data_CD34),file="CD34_strand_data.csv")

# Save to BED files to be suitable format for HOMER.

t321_peaklist <- read.csv("t321_specific_peaklist.csv")
CD34_peaklist <- read.csv("CD34_specific_peaklist.csv")

head(t321_peaklist)
  
write.table(t321_peaklist, file = 't321_spec_peaklist.bed', sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(CD34_peaklist, file = 'CD34_spec_peaklist.bed', sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

#---------------------------------------------------------------------------------
### Extracting Unique Peak IDs for HOMER fro full list and specific to t(3;21) and CD34+. 

bed<- strsplit2(row.names(resTable), split = '-')
head(bed)

# Write results to file
write.table(bed, file = 'CD34_t321_results_ordered.bed', sep = '-', quote = FALSE, row.names = FALSE, col.names = FALSE)

bed.CD34_specific<- strsplit2(row.names(resTable.CD34_specific), split = '-')
write.table(bed.CD34_specific, file = 'CD34_specific.bed', sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

bed.t321_specific<- strsplit2(row.names(resTable.t321_specific), split = '-')
write.table(bed.t321_specific, file = 't321_specific.bed', sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

###------------------------------------------------------------------------------
# Formatting peak files and subset data required for Wellington Footprint Analysis. 

t321_fp_peakslist <- t321_peaklist[c(1,2,3)]
write.csv(as.data.frame(t321_fp_peakslist),file="t321_fp_peakslist.csv")
t321_fp_peaklist <- read.csv("t321_fp_peakslist.csv")
head(t321_fp_peaklist)
write.table(t321_fp_peaklist, file = 't321_fp_list.bed', sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

CD34_fp_peakslist <- CD34_peaklist[c(1,2,3)]
write.csv(as.data.frame(CD34_fp_peakslist),file="CD34_fp_peakslist.csv")
CD34_fp_peaklist <- read.csv("CD34_fp_peakslist.csv")
head(CD34_fp_peaklist)
write.table(CD34_fp_peaklist, file = 'CD34_fp_peakslist.bed', sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

###------------------------------------------------------------------------------