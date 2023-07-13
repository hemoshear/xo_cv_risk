library(ggplot2)
library(plyr)
library(stringr)
library(reshape2)
library(ggsci)

# create results directory
out_dir <- "results/"

# Read in and reformat data -----------------------------------------------

# NOTE: This requires that both `pathway_correlation_analysis.R` and
# `pathway_correlation_analysis.scrambled.R` have been run

# get "normal" data
test_correlation_results_ec <- readRDS(paste0(out_dir, "test_gst_correlations_ec.RDS"))
test_correlation_results_smc <- readRDS(paste0(out_dir, "test_gst_correlations_smc.RDS"))

# reformat into dataframes
ec_pathways <- do.call("rbind", lapply(as.list(names(test_correlation_results_ec)), function(treatment){
  
  data <- test_correlation_results_ec[[treatment]][,1:11]
  
  data$treatment <- treatment
  
  data$Wilcox.p_value <- as.numeric(data$Wilcox.p_value)
  data$fisher.p_value <- as.numeric(data$fisher.p_value)
  
  data <- data[!is.na(data$Wilcox.p_value),]
  
  data$Wilcox.FDR <- p.adjust(data$Wilcox.p_value, method="fdr")
  data$fisher.FDR <- p.adjust(data$fisher.p_value, method="fdr")
  
  return(data)
  
}))

ec_pathways$cell <- "EC"


smc_pathways <- do.call("rbind", lapply(as.list(names(test_correlation_results_smc)), function(treatment){
  
  data <- test_correlation_results_smc[[treatment]][,1:11]
  
  data$treatment <- treatment
  
  data$Wilcox.p_value <- as.numeric(data$Wilcox.p_value)
  data$fisher.p_value <- as.numeric(data$fisher.p_value)
  
  data <- data[!is.na(data$Wilcox.p_value),]
  
  data$Wilcox.FDR <- p.adjust(data$Wilcox.p_value, method="fdr")
  data$fisher.FDR <- p.adjust(data$fisher.p_value, method="fdr")
  
  
  return(data)
  
}))

smc_pathways$cell <- "SMC"


total_pathways <- rbind(ec_pathways,
                        smc_pathways)

# get scrambled data
test_correlation_results_ec_s <- readRDS(paste0(out_dir, "test_gst_correlations_ec.scrambled.RDS"))
test_correlation_results_smc_s <- readRDS(paste0(out_dir, "test_gst_correlations_smc.scrambled.RDS"))

# reformat into dataframes
ec_pathways_s <- do.call("rbind", lapply(as.list(names(test_correlation_results_ec_s)), function(treatment){
  
  data <- test_correlation_results_ec_s[[treatment]][,1:11]
  
  data$treatment <- treatment
  
  data$Wilcox.p_value <- as.numeric(data$Wilcox.p_value)
  data$fisher.p_value <- as.numeric(data$fisher.p_value)
  
  data <- data[!is.na(data$Wilcox.p_value),]
  
  data$Wilcox.FDR <- p.adjust(data$Wilcox.p_value, method="fdr")
  data$fisher.FDR <- p.adjust(data$fisher.p_value, method="fdr")
  
  return(data)
  
}))

ec_pathways_s$cell <- "EC"


smc_pathways_s <- do.call("rbind", lapply(as.list(names(test_correlation_results_smc_s)), function(treatment){
  
  data <- test_correlation_results_smc_s[[treatment]][,1:11]
  
  data$treatment <- treatment
  
  data$Wilcox.p_value <- as.numeric(data$Wilcox.p_value)
  data$fisher.p_value <- as.numeric(data$fisher.p_value)
  
  data <- data[!is.na(data$Wilcox.p_value),]
  
  data$Wilcox.FDR <- p.adjust(data$Wilcox.p_value, method="fdr")
  data$fisher.FDR <- p.adjust(data$fisher.p_value, method="fdr")
  
  
  return(data)
  
}))

smc_pathways_s$cell <- "SMC"


total_pathways_s <- rbind(ec_pathways_s,
                        smc_pathways_s)

total_pathways_s$treatment <- paste0(total_pathways_s$treatment, "\nScrambled")


# determine significant pathways -------------------------------------------

total_sig_pathways <- total_pathways[total_pathways$Wilcox.p_value < 0.05 &
                                       total_pathways$fisher.p_value < 0.05,]

# significant for FBX
fbx_sig_pathways <- total_sig_pathways[total_sig_pathways$treatment == "FBX_100" &
                                         total_sig_pathways$cell == "SMC",]$pathway

# significant for FBX-1
fbx_1_sig_pathways <- total_sig_pathways[total_sig_pathways$treatment == "FBX-1_100" &
                                         total_sig_pathways$cell == "SMC",]$pathway


fbx_fbx_1_sig_pathways <- unique(c(fbx_sig_pathways,
                                 fbx_1_sig_pathways))

# Make one big log p plot ----------------------------------------------------------

treatments <- c("FBX_100","FBX-1_100","TPS_100",
                "OXP_500", "ATORV_ELD0105_D7", "FBX_100\nScrambled")

# merge all treatments
total_pathways$treatment <- as.character(total_pathways$treatment)
total_pathways_s$treatment <- as.character(total_pathways_s$treatment)

overlap_columns <- colnames(total_pathways)[colnames(total_pathways) %in% colnames(total_pathways_s)]

total_pathways_all <- rbind(total_pathways[,overlap_columns],
                            total_pathways_s[,overlap_columns])

# subset data for plotting
subset_p <- total_pathways_all[total_pathways_all$pathway %in% fbx_fbx_1_sig_pathways &
                                 total_pathways_all$treatment %in% treatments &
                                 total_pathways_all$cell == "SMC",]
subset_p$treatment <- factor(subset_p$treatment, levels=treatments)

subset_p$pathway_plot <- str_wrap(subset_p$pathway, width=50)

subset_p <- subset_p[order(subset_p$fisher.p_value),]
subset_p$pathway_plot <- factor(subset_p$pathway_plot, levels=unique(subset_p$pathway_plot))

# colors for each treatment
subset_colors <- c("FBX_100"="green4",
                   "FBX-1_100"="darkgoldenrod2",
                   "TPS_100"="dodgerblue2",
                   "OXP_500"="purple",
                   "ATORV_ELD0105_D7"="grey",
                   "FBX_100\nScrambled"="green")

# the plot
ggplot(subset_p, aes(x=pathway_plot,
                     y=-log(fisher.p_value),
                     group=treatment,
                     fill=treatment)) +
  geom_bar(stat="identity",
           position="dodge") +
  theme_bw() +
  scale_fill_manual(values=subset_colors) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  coord_flip() + scale_x_discrete(limits=rev(levels(subset_p$pathway_plot))) +
  geom_hline(yintercept=-log(0.05), color="black",linetype=2) +
  theme(legend.position = "bottom") +
  labs(x=NULL, y="-log(p) of Fisher's Exact Test\nfor Correlation with\nCV Risk Compounds", fill="Compound", shape="Count Type")
ggsave(paste0(out_dir, "neg_log_fisher.bar_plot.png"), width=11, height=17)









