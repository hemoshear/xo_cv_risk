library(parallel)
library(limma)
library(reshape2)
library(magrittr)

out_dir <- "results/"

dir.create(out_dir, showWarnings = F)

# Read in data ------------------------------------------------

# HRDB
hrdb_deg_ec1 <- readRDS("data/hrdb_deg_ec_1.RDS")
hrdb_deg_ec2 <- readRDS("data/hrdb_deg_ec_2.RDS")
hrdb_deg_ec <- c(hrdb_deg_ec1,
                 hrdb_deg_ec2)

hrdb_deg_smc1 <- readRDS("data/hrdb_deg_smc_1.RDS")
hrdb_deg_smc2 <- readRDS("data/hrdb_deg_smc_2.RDS")
hrdb_deg_smc <- c(hrdb_deg_smc1,
                  hrdb_deg_smc2)

# The test compounds
test_compounds_deg_ec <- readRDS("data/test_compounds_deg_ec.RDS")
test_compounds_deg_smc <- readRDS("data/test_compounds_deg_smc.RDS")

# scramble the data!
set.seed(42)

test_compounds_deg_ec <- lapply(test_compounds_deg_ec, function(x) {
  x$entrez_id <- sample(x$entrez_id)
  return(x)
})

test_compounds_deg_smc <- lapply(test_compounds_deg_smc, function(x) {
  x$entrez_id <- sample(x$entrez_id)
  return(x)
})

# curated reactome list
reactome_list <- readRDS("data/reactome_list.RDS")

# list of CV compounds
cv_risk_list <- readRDS("data/cv_risk_list.RDS")

# Pathway correlation function -----------------------------------------------------

pathway_correlations <- function(source_df,
                                 db_list,
                                 geneset,
                                 pval=0.05,
                                 min_gene=5) {
  # compare source set to each DB set
  result_list <- lapply(db_list, function(db) {
    
    # determine gene overlap
    # significant genes in each
    source_genes <- source_df[source_df$entrez_id %in% geneset &
                                source_df$PValue < pval,]$entrez_id
    db_genes <- db[db$entrez_id %in% geneset &
                     db$PValue < pval,]$entrez_id
    
    overlap_genes <- unique(c(source_genes, db_genes))
    # make sure they are expressed in both
    overlap_genes <- overlap_genes[overlap_genes %in% source_df$entrez_id &
                                     overlap_genes %in% db$entrez_id]
    
    # if not enough genes, return NA
    if (length(overlap_genes) < min_gene) {
      return(data.frame(cor_value=NA,
                        p_value=1,
                        deg_overlap=length(overlap_genes),
                        overlap_genes=""))
    }
    
    # determine correlation between them
    source_fc <- source_df[source_df$entrez_id %in% overlap_genes,]$logFC
    names(source_fc) <- source_df[source_df$entrez_id %in% overlap_genes,]$entrez_id
    
    db_fc <- db[db$entrez_id %in% overlap_genes,]$logFC
    names(db_fc) <- db[db$entrez_id %in% overlap_genes,]$entrez_id
    
    cor_result <- cor.test(source_fc[overlap_genes],
                           db_fc[overlap_genes], method = "spearman")
    
    # return results
    data.frame(cor_value=cor_result$estimate[[1]],
               p_value=cor_result$p.value[[1]],
               deg_overlap=length(overlap_genes),
               overlap_genes=paste0(overlap_genes, collapse=","))
    
  })
  
  # turn into dataframe
  result_df <- do.call("rbind", result_list)
  result_df$treatment_label <- rownames(result_df)
  
  # sort by correlation value
  result_df[order(result_df$cor_value, decreasing = T),]
  
}


# Calculate correlations -----------------------------------------------------

# compare each test treatment to HRDB treatment in EC
test2hrdb_correlations_ec <- lapply(test_compounds_deg_ec, function(source_df) {
  
  # for each pathway, calculate correlations
  results_list <- mclapply(reactome_list, function(geneset) {
    
    pathway_correlations(source_df,hrdb_deg_ec,geneset)
    
  }, mc.cores=30)
  
  # remove pathways with not enough genes for correlation
  keep <- unlist(lapply(results_list, function(x) {!all(is.na(x$cor_value))}))
  
  results_list[keep]
  
})

# compare each test treatment to HRDB treatment in SMC
test2hrdb_correlations_smc <- lapply(test_compounds_deg_smc, function(source_df) {
  
  # for each pathway, calculate correlations
  results_list <- mclapply(reactome_list, function(geneset) {
    
    pathway_correlations(source_df,hrdb_deg_smc,geneset)
    
  }, mc.cores=30)
  
  # remove pathways with not enough genes for correlation
  keep <- unlist(lapply(results_list, function(x) {!all(is.na(x$cor_value))}))
  
  results_list[keep]
  
})

# compare test treatments internally - EC
test2test_correlations_ec <- lapply(as.list(names(test_compounds_deg_ec)), function(treatment) {
  
  source_df <- test_compounds_deg_ec[[treatment]]
  db_list <- test_compounds_deg_ec[names(test_compounds_deg_ec) != treatment]
  
  # for each pathway, calculate correlations
  results_list <- mclapply(reactome_list, function(geneset) {
    
    pathway_correlations(source_df,db_list,geneset)
    
  }, mc.cores=30)
  
  # remove pathways with not enough genes for correlation
  keep <- unlist(lapply(results_list, function(x) {!all(is.na(x$cor_value))}))
  
  results_list[keep]
})
names(test2test_correlations_ec) <- names(test_compounds_deg_ec)

# compare BDG treatments internally - SMC
test2test_correlations_smc <- lapply(as.list(names(test_compounds_deg_smc)), function(treatment) {
  
  source_df <- test_compounds_deg_smc[[treatment]]
  db_list <- test_compounds_deg_smc[names(test_compounds_deg_smc) != treatment]
  
  # for each pathway, calculate correlations
  results_list <- mclapply(reactome_list, function(geneset) {
    
    pathway_correlations(source_df,db_list,geneset)
    
  }, mc.cores=30)
  
  # remove pathways with not enough genes for correlation
  keep <- unlist(lapply(results_list, function(x) {!all(is.na(x$cor_value))}))
  
  results_list[keep]
})
names(test2test_correlations_smc) <- names(test_compounds_deg_smc)

# genesettest correlation function -------------------------------------------

correlation_gst <- function(correlations_df,
                            treatment_list,
                            pval_threshold=0.05) {
  
  # remove NAs
  correlations_df <- correlations_df[!is.na(correlations_df$cor_value),]
  
  #treatment universe
  universe <- correlations_df$treatment_label %>% unique()
  
  # filter down to pvalue threshold
  correlations_df <- correlations_df[correlations_df$p_value <= pval_threshold,]
  
  # pull out all remaining treatments
  all_treatments <- correlations_df$treatment_label
  
  # perform gene set test
  pval <- geneSetTest(all_treatments %in% treatment_list,
                      correlations_df$cor_value, type="t", alternative = "up")
  
  # determine mean correlations
  avg_cor <- mean(correlations_df$cor_value)
  cv_risk_cor <- mean(correlations_df[correlations_df$treatment_label %in% treatment_list,]$cor_value)
  non_cv_risk_cor <- mean(correlations_df[!correlations_df$treatment_label %in% treatment_list,]$cor_value)
  
  if (!is.na(cv_risk_cor) & !is.na(non_cv_risk_cor)) {
    wilcox_result <- wilcox.test(correlations_df[correlations_df$treatment_label %in% treatment_list,]$cor_value,
                                 correlations_df[!correlations_df$treatment_label %in% treatment_list,]$cor_value,
                                 alternative = "greater")
    wilcox_pval <- wilcox_result$p.value
  } else {
    wilcox_pval <- 1
  }
  
  # fisher test
  a <- universe %in% dplyr::filter(correlations_df, cor_value > 0)$treatment_label
  b <- universe %in% treatment_list
  if (sum(a) != 0 & sum(b) != 0) {
    fisher_pval <- fisher.test(a, b, alternative = "greater")$p.value
  } else {
    fisher_pval <- 1
  }
  
  # pull out all CV treatments
  cv_treatments <- correlations_df$treatment_label[correlations_df$treatment_label %in% treatment_list]
  # NEED TO BE POSITIVELY CORRELATED
  pos_cv_treatments <- correlations_df[correlations_df$treatment_label %in% treatment_list &
                                         correlations_df$cor_value > 0,]$treatment_label
  
  # output results
  data.frame(treatment_count=length(all_treatments),
             cv_treatment_count=sum(all_treatments %in% treatment_list),
             pos_cv_treatment_count=length(pos_cv_treatments),
             GST.p_value=pval,
             Wilcox.p_value=wilcox_pval,
             fisher.p_value=fisher_pval,
             avg_cor=avg_cor,
             cv_risk_cor=cv_risk_cor,
             non_cv_risk_cor=non_cv_risk_cor,
             cv_treatments=paste0(pos_cv_treatments, collapse=","))
  
}

# gene set test results on correlations -------------------------------------

# EC
test2hrdb_correlations_ec_p <- lapply(test2hrdb_correlations_ec, function(correlation_list) {
  
  gst_results <- mclapply(correlation_list, function(correlations_df) {
    correlation_gst(correlations_df, cv_risk_list)
  }, mc.cores = 30)
  
  gst_results_df <- do.call("rbind", gst_results)
  gst_results_df$pathway <- rownames(gst_results_df)
  
  return(gst_results_df)
})

# SMC
test2hrdb_correlations_smc_p <- lapply(test2hrdb_correlations_smc, function(correlation_list) {
  
  gst_results <- mclapply(correlation_list, function(correlations_df) {
    correlation_gst(correlations_df, cv_risk_list)
  }, mc.cores=30)
  
  gst_results_df <- do.call("rbind", gst_results)
  gst_results_df$pathway <- rownames(gst_results_df)
  
  return(gst_results_df)
})

# combine gene set test results with test2test correlations -----------

# EC
test_correlation_results_ec <- lapply(as.list(names(test2hrdb_correlations_ec_p)), function(treatment) {
  
  gst_results_df <- test2hrdb_correlations_ec_p[[treatment]] # geneset test results
  test_correlations <- test2test_correlations_ec[[treatment]] # correlations with other test compounds
  
  # convert into dataframe
  test_correlations_df <- do.call("rbind", lapply(as.list(names(test_correlations)), function(pathway) {
    data <- test_correlations[[pathway]]
    data$pathway <- pathway
    
    return(data)
  }))
  
  # convert to wide form
  test_correlations_df <- dcast(test_correlations_df, pathway ~ treatment_label, value.var="cor_value")
  colnames(test_correlations_df)[-1] <- paste0(colnames(test_correlations_df)[-1], ".cor")
  
  # merge em up!
  gst_results_df <- merge(gst_results_df, test_correlations_df, by="pathway", all.x=T)
  
  gst_results_df[order(gst_results_df$Wilcox.p_value),]
})
names(test_correlation_results_ec) <- names(test2hrdb_correlations_ec_p)

# SMC
test_correlation_results_smc <- lapply(as.list(names(test2hrdb_correlations_smc_p)), function(treatment) {
  
  gst_results_df <- test2hrdb_correlations_smc_p[[treatment]] # geneset test results
  test_correlations <- test2test_correlations_smc[[treatment]] # correlations with other test compounds
  
  # convert into dataframe
  test_correlations_df <- do.call("rbind", lapply(as.list(names(test_correlations)), function(pathway) {
    data <- test_correlations[[pathway]]
    data$pathway <- pathway
    
    return(data)
  }))
  
  # convert to wide form
  test_correlations_df <- dcast(test_correlations_df, pathway ~ treatment_label, value.var="cor_value")
  colnames(test_correlations_df)[-1] <- paste0(colnames(test_correlations_df)[-1], ".cor")
  
  # merge em up!
  gst_results_df <- merge(gst_results_df, test_correlations_df, by="pathway", all.x=T)

  gst_results_df[order(gst_results_df$Wilcox.p_value),]
})
names(test_correlation_results_smc) <- names(test2hrdb_correlations_smc_p)

# output
saveRDS(test_correlation_results_ec, file=paste0(out_dir, "test_gst_correlations_ec.scrambled.RDS"))
saveRDS(test_correlation_results_smc, file=paste0(out_dir, "test_gst_correlations_smc.scrambled.RDS"))


