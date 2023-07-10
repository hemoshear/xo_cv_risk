library(ggplot2)
library(reshape2)
library(ggbeeswarm)
library(stringr)
library(snakecase)

out_dir <- "results/"

dir.create(out_dir, showWarnings = F)

meta <- readRDS("data/fbx_oxp_count_metadata.RDS")
meta$treatment <- factor(as.character(meta$treatment),
                         levels=c("ICS VEH",
                                  "500uM OXP",
                                  "100uM FBX"))

pathway_counts_ec <- readRDS("data/fbx_oxp_reactome_pathway_counts.ec.RDS")
pathway_counts_smc <- readRDS("data/fbx_oxp_reactome_pathway_counts.smc.RDS")

pathway_dep_ec <- readRDS("data/fbx_oxp_reactome_pathway_dep.ec.RDS")
pathway_dep_smc <- readRDS("data/fbx_oxp_reactome_pathway_dep.smc.RDS")

pathway_dep_ec <- do.call("rbind", lapply(names(pathway_dep_ec), function(contrast) {
  data <- pathway_dep_ec[[contrast]]
  data$pathway <- rownames(data)
  data$contrast <- contrast
  
  return(data)
}))

pathway_dep_smc <- do.call("rbind", lapply(names(pathway_dep_smc), function(contrast) {
  data <- pathway_dep_smc[[contrast]]
  data$pathway <- rownames(data)
  data$contrast <- contrast
  
  return(data)
}))

# list of pathways
pathway_list_ec <- c("eNOS activation",
                     "ROS, RNS production in phagocytes",
                     "TNFs bind their physiological receptors")

pathway_list_smc <- c("RHO GTPases activate PKNs",
                      "Calmodulin induced events",
                      "Synthesis of Leukotrienes (LT) and Eoxins (EX)")


pathway_box_plot <- function(pathway,
                             dep,
                             counts,
                             metadata) {
  
  
  counts_long <- melt(counts)  
  colnames(counts_long) <- c("pathway", "sample_id", "value")
  
  counts_long <- counts_long[counts_long$pathway == pathway,]
    
  counts_long <- merge(counts_long, metadata, by="sample_id")
  
  dep <- dep[dep$pathway == pathway,]
  dep$treatment <- sapply(dep$contrast, function(x) {unlist(strsplit(x, "-"))[1]})
  dep$treatment <- factor(dep$treatment, levels=levels(metadata$treatment))

  dep$plot_value <- paste0("logFC: ", round(dep$logFC, 3), 
                           "\nFDR: ", formatC(dep$adj.P.Val, format="e", digits=2))  
  
  counts_long <- merge(counts_long,
                  dep[,c("treatment","plot_value")],
                  all.x=T)
  
  counts_long[is.na(counts_long$plot_value),]$plot_value <- ""
  
  counts_long <- counts_long[order(counts_long$treatment),]
  counts_long$plot_value <- factor(counts_long$plot_value,
                                     levels = unique(counts_long$plot_value))
  
  ggplot(counts_long, aes(x=plot_value,
                            y=value,
                            color=treatment)) +
    geom_boxplot(outlier.shape = NA) +
    geom_quasirandom(size=1, width=0.1) +
    theme_bw() +
    scale_color_manual(values=c("black","purple","green4")) +
    theme(legend.position = "bottom") +
    scale_x_discrete(position="top") +
    theme(text = element_text(size = 7),
          axis.text.x = element_text(angle=30, hjust = 0.1, vjust=0),
          plot.title=element_text(hjust=0.5, size=10)) +
    labs(x=NULL, color=NULL, title=str_wrap(pathway, width = 30), y="GSVA Score")
  
}

for (pathway in pathway_list_ec) {
  
  outfile <- paste0(out_dir, to_snake_case(pathway), ".ec_box_plot.png")
  pathway_box_plot(pathway,
                   pathway_dep_ec,
                   pathway_counts_ec,
                   meta[meta$cell_type == "EC",])
  ggsave(outfile, width=3, height=3)
  
  
}

for (pathway in pathway_list_smc) {
  
  outfile <- paste0(out_dir, to_snake_case(pathway), ".smc_box_plot.png")
  pathway_box_plot(pathway,
                   pathway_dep_smc,
                   pathway_counts_smc,
                   meta[meta$cell_type == "SMC",])
  ggsave(outfile, width=3, height=3)
  
  
}

