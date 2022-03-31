
#ggplot2_3.3.3
#install.packages("vegan", dependencies=TRUE)
library(ggplot2)

output_dir <- "stacked_barplot_output_files"
if (!dir.exists(output_dir)){
    dir.create(output_dir)
} else {
    print("Dir already exists!")
}

data <- read.csv("/Users/kevin.muirhead/OneDrive - University of Calgary/Kevin_Digital_Notebook/tkg_pilot_project/2021_pullulanase_paper/figures/thika_pilot_pheatmap_2021-11-23/files/THIKA_metaphlan_merged_bacteria_cohort_ra_thresh_0.01_for_R_stacked_barplot.tsv", header = TRUE,sep = "\t")
data <- read.csv("stacked_barplot_input_files/THIKA_metaphlan_merged_bacteria_cohort_ra_thresh_0.01_for_R_stacked_barplot.tsv", header = TRUE,sep = "\t")


data$taxonomy = gsub("_", " ", data$taxonomy)
data$sample_id = as.character(data$sample_id)

# Order of the bray-curtis distance and "average" hclust dendrogram hardcoded.
#order = c("16","3","17","5","7","1","12","13","8","10","15","14", "6", "11", "2", "9", "4")

# Order of the bray-curtis distance and "ward.D" hclust dendrogram hardcoded.
order = c("11","15","13","14","10","8","7","2","12","6","1","9", "4", "17", "5", "16", "3")

tiff(paste(output_dir, "THIKA_AGYW_pilot_hclust_ward_d_stacked_barplot.tiff", sep = "/"), width=6000, height=2400, unit="px", res=300)

data$taxonomy <- reorder(data$taxonomy, data$relative_abundance)
data$taxonomy <- factor(data$taxonomy, levels=levels(data$taxonomy))

# Stacked + percent
ggplot(data, aes(fill=taxonomy, y=relative_abundance, x=factor(sample_id,order), width = 1)) +
    geom_bar(colour="black",position="fill", stat="identity") + xlab("Sample Name") + ylab("Relative Abundance (%)") + theme(text = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_rect(colour = "black", fill=NA, size=0.25), legend.position="bottom") + labs(fill = "Taxonomy") + scale_fill_manual(values = c("#3a93a0", "#ff216b", "#27235a", "#d6aea1", "#006a46", "#786b2a", "#561a7c", "#7a1035", "#c5a2d4", "#08ddfc", "#9c7cb0", "#31aec0", "#7aaaec", "#ffa900"), guide = guide_legend(reverse = TRUE)) + scale_y_continuous(labels = scales::percent_format())

