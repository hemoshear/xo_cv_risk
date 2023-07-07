library(parallel)

out_dir <- "results/"

dir.create(out_dir, showWarnings = F)

# Read in data ------------------------------------------------
hrdb_deg_ec <- readRDS("data/hrdb_deg_ec.RDS")
hrdb_deg_smc <- readRDS("data/hrdb_deg_smc.RDS")

test_compounds_deg_ec <- readRDS("data/test_compounds_deg_ec.RDS")
test_compounds_deg_smc <- readRDS("data/test_compounds_deg_smc.RDS")

reactome_list <- readRDS("data/reactome_list.RDS")

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

# compare each BDG treatment to HRDB treatment in EC
test2hrdb_correlations_ec <- lapply(test_compounds_deg_ec, function(source_df) {
  
  # for each pathway, calculate correlations
  results_list <- mclapply(reactome_list, function(geneset) {
    
    pathway_correlations(source_df,hrdb_deg_ec,geneset)
    
  }, mc.cores=30)
  
  # remove pathways with not enough genes for correlation
  keep <- unlist(lapply(results_list, function(x) {!all(is.na(x$cor_value))}))
  
  results_list[keep]
  
})

# compare each BDG treatment to HRDB treatment in SMC
bdg2hrdb_correlations_smc <- lapply(bdg_deg_smc, function(source_df) {
  
  # for each pathway, calculate correlations
  results_list <- mclapply(reactome.list, function(geneset) {
    
    pathway_correlations(source_df,eld_deg_smc,geneset)
    
  }, mc.cores=30)
  
  # remove pathways with not enough genes for correlation
  keep <- unlist(lapply(results_list, function(x) {!all(is.na(x$cor_value))}))
  
  results_list[keep]
  
})

