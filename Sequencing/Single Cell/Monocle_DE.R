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
