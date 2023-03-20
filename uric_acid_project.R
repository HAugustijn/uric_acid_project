# Install/load packages
# If 'metagenomeSeq' is not yet installed locally, use Biocmanager::install('metagenomeSeq')
packages = c("metagenomeSeq","biomformat","ComplexHeatmap","viridisLite", 
             "RColorBrewer","dplyr", "stringr", "ggpubr")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    print("not installed") 
    install.packages(x, dependencies = TRUE)
  library(x, character.only = TRUE)
  suppressWarnings(suppressMessages(library(x, character.only = TRUE, quietly = T)))
  }
})

# Load biom file
biom_file <- "BiG-MAP.map.meta.dec.biom"
MR <- loadBiom(biom_file)

# Normalize the data using cumulative sum scaling (CSS)
p <- cumNormStat(MR)
MR_sample <- cumNorm(MR, p=p)
countsdf <- data.frame(MRcounts(MR_sample, norm=T, log=T))
countsdf[countsdf == 'NA'] <- NA

# Load and order the metadata
metadata <- pData(MR_sample)
metadata[metadata == 'NA'] <- NA

# Order the data on the HostDiet, CurrentAntibiotics and StudyDay
countsdf_mt <- data.frame(t(rbind(countsdf, t(metadata))))
countsdf_mt$StudyDay <- as.numeric(countsdf_mt$StudyDay)
countsdf_mt <- countsdf_mt[order(countsdf_mt$CurrentAntibiotics, countsdf_mt$HostDiet, countsdf_mt$StudyDay),]
countsdf_mt <- na.omit(countsdf_mt)
countsdf_mt[ ,c('HostDiet', 'CurrentAntibiotics', 'HostAge', 'HostBMI', 'HostHeight', 
                'Weight', 'StudyDay', 'host.ID', 'SampleType', 'SubjectID')] <- list(NULL)
metadata <- na.omit(metadata)
metadata <- metadata[order(metadata$CurrentAntibiotics, metadata$HostDiet, metadata$StudyDay),]

# Rename the columns to ID and species name
colnames(countsdf_mt) <- gsub("\\..SMASH(.*)","", colnames(countsdf_mt))
colnames(countsdf_mt) <- gsub("(.*)\\.GC_DNA..Entryname.uric_acid..OS.","", colnames(countsdf_mt))
colnames(countsdf_mt) <- gsub("\\.","", colnames(countsdf_mt))
colnames(countsdf_mt) <- gsub("\\__"," ", colnames(countsdf_mt))
colnames(countsdf_mt) <- gsub("\\_"," ", colnames(countsdf_mt))
metadata$CurrentAntibiotics <- gsub("(^(?)A|B|C)","", metadata$CurrentAntibiotics)
metadata$CurrentAntibiotics <- gsub("\\_"," ", metadata$CurrentAntibiotics)
countsdf_mt <- countsdf_mt[,order(colnames(countsdf_mt))]

countsmatrix <- apply(as.matrix.noquote(countsdf_mt, labels=TRUE),2,as.numeric)
rownames(countsmatrix) <- rownames(countsdf_mt)
countsmatrix <- t(countsmatrix)


# Set annotation
ann <- data.frame(metadata$CurrentAntibiotics, metadata$HostDiet)
colnames(ann) <- c('CurrentAntibiotics', 'HostDiet')
colours <- list('HostDiet' = c('EEN' = '#3c5aa8', 
                               'Vegan' = '#7392c9', 
                               'Omnivore' = '#15386b'),
      'CurrentAntibiotics' = c('Pre Antibiotics' = '#3c5aa8', 
                               'Antibiotics Treatment' = '#7392c9', 
                               'Post Antibiotics' = '#15386b'))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))

# Create heatmap
hmap <- Heatmap(
  countsmatrix,
  name = "Abundance",
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  row_names_gp = grid::gpar(fontsize = 6),
  top_annotation=colAnn)

draw(hmap, heatmap_legend_side="left", annotation_legend_side="right")

#------------------------------------------------------------------------------#

# Sum the RPKM values for each sample
sum_data <- data.frame(colSums(countsmatrix))
sum_meta <- data.frame(t(rbind(t(sum_data), t(metadata))))
sum_meta$colSums.countsmatrix. <- as.numeric(sum_meta$colSums.countsmatrix.)
sum_meta$HostAge <- as.numeric(sum_meta$HostAge)
sum_meta$HostBMI <- as.numeric(sum_meta$HostBMI)
colnames(sum_meta)[1] <- "SumMeta"

# Correcting for confounding variables
m_full <- lm(SumMeta ~ HostAge + HostBMI, data = sum_meta)
sum_meta$residuals <- residuals(m_full)

# Comparing the RPKM values across the different diets using a pairwise comparisons p-value and anova
my_comparisons <- list( c("EEN", "Omnivore"), c("EEN", "Vegan"), c("Vegan", "Omnivore") )

ggboxplot(sum_meta, x = "HostDiet", y = "residuals",
          color = "HostDiet",
          fill = "HostDiet", alpha=0.2,
          palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          facet.by = "CurrentAntibiotics") +
  xlab("") + ylab("") +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     label.y = c(220, 270, 200)) +
  stat_compare_means(method = "anova", label.y = -300) +
  geom_jitter(alpha=0.08, width=0.05)




