################################################################
visr.category("Differential expression analysis",
              info = "Characterize differentially expressed genes by comparing groups of cells",
              active.condition = "visr.param.output_dir != '' && visr.param.enable_de_analysis")
################################################################

DE_FORMULA_ALL_CLUSTERS = "across all cell clusters"
DE_FORMULA_ONE_CLUSTER = "by one cluster versus the others"
DE_FORMULA_PHENOTYPE = "by phenotype (cell attribute)"
visr.param("de_formula", label = "Perform DE", items = c(DE_FORMULA_ALL_CLUSTERS, DE_FORMULA_ONE_CLUSTER, DE_FORMULA_PHENOTYPE))

visr.param("de_cluster_id", label = "Perform DE on which Cluster id?", min = 1L, default = 1L,
           items = c(Inf), item.labels = c("All (slow)"),
           active.condition = sprintf("visr.param.de_formula == '%s'", DE_FORMULA_ONE_CLUSTER))

visr.param("de_phenotype_name", label = "Phenotype (cell attribute) name",
           active.condition = sprintf("visr.param.de_formula == '%s'", DE_FORMULA_PHENOTYPE))

visr.param("de_subset_by_marker_genes", label = "Perform DE on selected marker genes only",
           info = "You can select a specific set of genes that you know are important for your analysis. Otherwise will perform DE on only genes that pass the filtering process",
           default = F)

visr.param("marker_genes_list", label = "Marker genes (comma separated)", info="Optional comma separated list of short gene names to use as marker genes.",
           debugvalue = "MEF2C, MEF2D, MYF5, ANPEP, PDGFRA, MYOG, TPM1, TPM2, MYH2, MYH3, NCAM1, TNNT1, TNNT2, TNNC1, CDK1, CDK2, CCNB1, CCNB2, CCND1, CCNA1, ID1",
           active.condition = "visr.param.de_subset_by_marker_genes")

visr.param("num_plot_genes_jitter", label = "Draw level of expression for how many top genes?", default = 9L,
           info = "Plots the level of expression for each group of cells per gene,\nfor the specified number of most statistically significant genes.") # min = 0

################################################################
################################################################
################################################################

#'
#' detects de genes given model formula
#' @param cds CellDataSet       object
#' @param fullModelFormulaStr   formula string specifying the full model in differential expression test
#' @param output_de_genes_filename filename to output the de genes (optional)
detect_de_genes <- function(cds, fullModelFormulaStr, output_de_genes_filename = NULL) {
  stopifnot(is(cds, "CellDataSet"))

  de_genes <- monocle::differentialGeneTest(cds,
                                            fullModelFormulaStr = fullModelFormulaStr,
                                            reducedModelFormulaStr = "~1", # default
                                            relative_expr = TRUE, # default
                                            cores = visr.param.num_cores,
                                            verbose = TRUE)
  dim(de_genes)
  #de_genes_valid <- de_genes[which(de_genes$qval < visr.param.trajectory_max_qval),]
  #de_genes_valid <- de_genes_valid[order(de_genes_valid$qval),]
  de_genes_valid <- de_genes[order(de_genes$qval),]

  if (!is.null(output_de_genes_filename)) {
    visr.logProgress(paste("Writing DE genes to file", output_de_genes_filename))
    visr.writeDataTable(de_genes_valid, output_de_genes_filename)
  }

  # the most statistically significant genes
  visr.logProgress(paste("Plotting level of expression for top", visr.param.num_plot_genes_jitter, "most statistically significant gene(s)."))
  p <- monocle::plot_genes_jitter(cds[de_genes_valid$id[seq(visr.param.num_plot_genes_jitter)],],
                                  grouping = "Cluster", color_by = "Cluster",
                                  ncol = ceiling(sqrt(visr.param.num_plot_genes_jitter)), nrow = NULL)
  print(p)

  return (de_genes_valid)
}

#'
#' performs differential expression analysis
#'
perform_de_analysis <- function(monocle_app_object, output_dir) {
  my_cds <- monocle_app_object$cds
  stopifnot(is(my_cds, "CellDataSet"))

  if (visr.var.add_title_pages)
    plotTitlePage("Differential Expression Analysis", 2)

  if (visr.param.de_subset_by_marker_genes & visr.param.marker_genes_list != '') {
    marker_genes <- row.names(subset(fData(my_cds), gene_short_name %in% strsplit(visr.param.marker_genes_list, split = "[ ]*[,| ]+[ ]*")[[1]]))
    cds_de_subset <- my_cds[marker_genes,]
  } else {
    cds_de_subset <- my_cds[monocle_app_object$subset_gene_ids,]
  }

  if (visr.param.de_formula == DE_FORMULA_ALL_CLUSTERS) {
    visr.logProgress("Performing differential expression analysis across all clusters ...")
    de_genes_filename <- paste0(output_dir, "/de_genes_all_clusters.txt")
    de_genes <- detect_de_genes(cds_de_subset, '~Cluster', de_genes_filename)
  }
  else if (visr.param.de_formula == DE_FORMULA_ONE_CLUSTER) {
    # create vector of no's
    de_cluster_ids = if (visr.param.de_cluster_id == Inf) levels(pData(my_cds)$Cluster) else visr.param.de_cluster_id
    i <- 0
    for (de_cluster_id in de_cluster_ids) {
      i <- i + 1
      # change status to yes if the cell was in cluster 1
      my_vector <- rep('no', nrow(pData(my_cds)))
      my_vector[pData(my_cds)$Cluster == de_cluster_id] <- 'yes' #rep('yes', sum(pData(my_cds)$Cluster == de_cluster_id))

      # add vector to phenoData
      pData(my_cds)$test <- my_vector

      head(pData(my_cds))

      # TODO: perform DE based on a selected gene (replace ~test with the specified gene expresion column)
      # question: should we use normalized or unnormalized expression values
      visr.logProgress(paste("(", i, "of", length(de_cluster_ids), ")",
                             "Performing the differential expression analysis on cluster", de_cluster_id,
                             "\nfor the", nrow(unsup_clustering_genes)))

      de_genes_filename <- paste0(output_dir, "/de_genes_cluster", de_cluster_id, ".txt")

      de_genes <- detect_de_genes(cds_de_subset, '~test', de_genes_filename)
    }
  }
  else if (visr.param.de_formula == DE_FORMULA_PHENOTYPE) {
    de_genes_filename <- paste0(output_dir, "/de_genes_", visr.param.de_phenotype_name, ".txt")
    de_genes <- detect_de_genes(cds_de_subset, paste0('~', visr.param.de_phenotype_name), de_genes_filename)
  }

  monocle_app_object$cds <- my_cds
  monocle_app_object$de_genes <- de_genes
  return(monocle_app_object)
}
