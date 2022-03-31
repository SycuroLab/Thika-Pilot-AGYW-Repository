#pheatmap_1.0.12
#install.packages("pheatmap", dependencies=TRUE)
library(pheatmap)

#RColorBrewer_1.1-2
#install.packages("RColorBrewer", dependencies=TRUE)
library(RColorBrewer)

# vegan_2.5-7
#install.packages("vegan", dependencies=TRUE)
library(vegan)

#install.packages("ctc", dependencies=TRUE)
library(ctc)

output_dir <- "heatmap_output_files"
if (!dir.exists(output_dir)){
    dir.create(output_dir)
} else {
    print("Dir already exists!")
}

dataframe <- data.matrix(read.csv("heatmap_input_files/THIKA_metaphlan_merged_abundance_counts_bacteria.txt", header = TRUE, row.names = 1,sep = "\t",check.names = FALSE))
#dataframe

metadataframe <- data.matrix(read.csv("heatmap_input_files/ThikaGHS_pilot17_heatmap.csv", header = TRUE, row.names = 1,sep = ","))
#metadataframe
#q()

# define the annotation
annotation_col = data.frame(
total_reads = as.numeric(metadataframe[,"total_reads"]),
age = metadataframe[,"age"],
days_since_lmp = metadataframe[,"days_since_lmp"],
sti_yeast = metadataframe[,"sti_yeast"],
nugent_score = metadataframe[,"nugent_score"],
sexual_debut = as.character(metadataframe[,"sexual_debut"])
#pula_gene = as.character(metadataframe[,"pula_gene"])
)

as.character(metadataframe[,"sti_yeast"])

#annotation_col

# change the color of annotation to what you want: (eg: "navy", "darkgreen")
sexual_debut <- c("#FF0000", "#0d8548")
names(sexual_debut) <- c(1, 2)

total_reads <- colorRampPalette(c("#d0b4dc", "#945cb4", "#1b047c"))(100)
names(total_reads) <- annotation_col$total_reads

age <- colorRampPalette(c("#ddefdd", "#9ace9a", "#5a8b4f"))(100)
names(age) <- annotation_col$age

nugent_score <- colorRampPalette(c("#C1C1C1", "#616161", "#464646"))(100)
names(nugent_score) <- annotation_col$nugent_score

days_since_lmp <- colorRampPalette(c("#7bccc8", "#46b4af", "#4682b4"))(100)
names(days_since_lmp) <- annotation_col$days_since_lmp

#sti_yeast <- c(yeast = "white", trichomoniasis="blue")
#names(sti_yeast) <- annotation_col$sti_yeast

sti_yeast <- c("#FFFFFF","#4682b4", "#F6BE00")
names(sti_yeast) <- c(1, 2, 3)

#pula_gene <- c("#4682b4", "#F6BE00")
#names(pula_gene) <- c(1, 2)

annotation_col

#anno_colors <- list(total_reads = total_reads, age = age, days_since_lmp=days_since_lmp, sti_yeast=sti_yeast,nugent_score = nugent_score, sexual_debut = sexual_debut, pula_gene=pula_gene)
#anno_colors <- list(total_reads = total_reads, age = age, days_since_lmp=days_since_lmp, sti_yeast=sti_yeast, nugent_score = nugent_score, sexual_debut = sexual_debut, pula_gene=pula_gene)

anno_colors <- list(total_reads = total_reads, age = age, days_since_lmp=days_since_lmp, sti_yeast=sti_yeast, nugent_score = nugent_score, sexual_debut = sexual_debut)
#dist <- dist(t(dataframe), method = "euclidean")   # find distance matrix

# Bray-Curtis dissimilarity matrix

dist <- vegdist(t(dataframe), method = "bray")

# average = UPGMA
hclust <- hclust(dist, method="average")
tiff(paste(output_dir, "THIKA_AGYW_pilot_hclust_average_dendrogram.tiff", sep = "/"), width=6000, height=3600, unit="px", res=300)
plot(hclust)
dev.off()

# ward.D
hclust_ward_d <- hclust(dist, method="ward.D")
tiff(paste(output_dir, "THIKA_AGYW_pilot_hclust_ward_d_dendrogram.tiff", sep = "/"), width=6000, height=3600, unit="px", res=300)
plot(hclust_ward_d)
dev.off()

# ward.D2
hclust_ward_d2 <- hclust(dist, method="ward.D2")
tiff(paste(output_dir, "THIKA_AGYW_pilot_hclust_ward_d2_dendrogram.tiff", sep = "/"), width=6000, height=3600, unit="px", res=300)
plot(hclust_ward_d2)
dev.off()


newick <- hc2Newick(hclust, flat=TRUE)
write(hc2Newick(hclust, flat=TRUE),file=paste(output_dir,  "THIKA_AGYW_pilot_hclust_dendrogram.newick.txt", sep = "/"))

# Plot the heatmap using hclust with "average" method.
tiff(paste(output_dir, "THIKA_AGYW_pilot_hclust_average_pheatmap.tiff", sep = "/"), width=6000, height=3600, unit="px", res=300)
pheatmap_object=pheatmap(dataframe,cluster_rows = FALSE, annotation_col = annotation_col,annotation_colors=anno_colors, treeheight_column = 1000,
cluster_cols = hclust, cellwidth = 50, cellheight=10, angle_col = 0,fontsize_number = 50, border_color = "black", na_col = "white")
dev.off()

# Plot the heatmap using hclust with "ward.D" method.
tiff(paste(output_dir, "THIKA_AGYW_pilot_hclust_ward_d_pheatmap.tiff", sep = "/"), width=6000, height=3600, unit="px", res=300)
pheatmap_object=pheatmap(dataframe,cluster_rows = FALSE, annotation_col = annotation_col,annotation_colors=anno_colors, treeheight_column = 1000,
cluster_cols = hclust_ward_d, cellwidth = 50, cellheight=10, angle_col = 0,fontsize_number = 50, border_color = "black", na_col = "white")
dev.off()

# Plot the heatmap using hclust with "ward.D2" method.
tiff(paste(output_dir, "THIKA_AGYW_pilot_hclust_ward_d2_pheatmap.tiff", sep = "/"), width=6000, height=3600, unit="px", res=300)
pheatmap_object=pheatmap(dataframe,cluster_rows = FALSE, annotation_col = annotation_col,annotation_colors=anno_colors, treeheight_column = 1000,
cluster_cols = hclust_ward_d2, cellwidth = 50, cellheight=10, angle_col = 0,fontsize_number = 50, border_color = "black", na_col = "white")
dev.off()
