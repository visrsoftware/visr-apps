################################################################
visr.category(label="Input")
################################################################

INPUT_TYPE_10X = "10X single cell dataset"
INPUT_TYPE_COUNT_MATRIX_TXT = "Count matrix (.txt)"
INPUT_TYPE_SPARSE_MATRIX = "Sparse matrix"
visr.param("input_type",
           info = "Specify the format of your input data. Whether it is a 10X cell ranger dataset or a normal or sparse an expression matrix",
           items = c(INPUT_TYPE_10X, INPUT_TYPE_COUNT_MATRIX_TXT, INPUT_TYPE_SPARSE_MATRIX), debugvalue = INPUT_TYPE_10X)

visr.param("data_dir_10x",
           label="10X dataset directory",
           info="10X Cell Ranger dataset directory. It should contain a sub-directory named 'outs'",
           type="filename", filename.mode = "dir",
           active.condition = sprintf("visr.param.input_type == '%s'", INPUT_TYPE_10X),
           debugvalue= "~/SFU/Datasets/SingleCell/pbmc3k/") # Peripheral Blood Mononuclear Cells (PBMCs) from a healthy donor

visr.param("expression_matrix", type = "filename",
           info = "Numeric matrix of expression values, where rows are genes, and columns are cells",
           active.condition = sprintf("(visr.param.input_type == '%s') || (visr.param.input_type == '%s')", INPUT_TYPE_COUNT_MATRIX_TXT, INPUT_TYPE_SPARSE_MATRIX))

visr.param("sample_sheet", type = "filename",
           info = "Pheno data where rows are cells, and columns are cell attributes (such as cell type, culture condition, day captured, etc.). It should have same number of rows as the number of columns in expression value matrix.",
           active.condition = sprintf("(visr.param.input_type == '%s') || (visr.param.input_type == '%s')", INPUT_TYPE_COUNT_MATRIX_TXT, INPUT_TYPE_SPARSE_MATRIX))

visr.param("gene_annotation", type = "filename",
           info = "Feature data where rows are features (e.g. genes), and columns are gene attributes (such as biotype, gc content, etc.). It should have same number of rows as the number of rows in expression value matrix.",
           active.condition = sprintf("(visr.param.input_type == '%s') || (visr.param.input_type == '%s')", INPUT_TYPE_COUNT_MATRIX_TXT, INPUT_TYPE_SPARSE_MATRIX))

DATA_TYPE_UMI = "UMIs, Raw transcript counts"
DATA_TYPE_FPKM = "FPKM, TPM"
DATA_TYPE_LOG_FPKM = "log FPKM/TPMs, Ct values from SC qPCR"
visr.param("data_type",
           label = "Data type",
           info = "Specify the appropriate data type. Monocle works well with both count-based measures (e.g. UMIs) and relative expression data.",
           items=c(DATA_TYPE_UMI, DATA_TYPE_FPKM, DATA_TYPE_LOG_FPKM),
           active.condition = "(visr.param.data_dir_10x != '') || (visr.param.expression_matrix != '' && visr.param.sample_sheet != '' && visr.param.gene_annotation != '')")

################################################################
visr.category("Output",
              active.condition = "(visr.param.data_dir_10x != '') || (visr.param.expression_matrix != '' && visr.param.sample_sheet != '' && visr.param.gene_annotation != '')")
################################################################

visr.param("output_dir",
           label="Output Directory to save the results",
           info="Output directory where the analysis results will be saved to",
           type="filename", filename.mode = "dir",
           debugvalue= "~/SFU/Datasets/SingleCell/output_monocle/")

visr.param("create_subdir", default = TRUE,
           label = "Create new sub-direcory",
           info = "Create a new sub directory with the name DATE_TIME (YYYYMMDD_hhmmss)")

################################################################
visr.category("Analysis steps",
              info = "Different analysis steps",
              active.condition = "visr.param.output_dir != '' && ((visr.param.data_dir_10x != '') || (visr.param.expression_matrix != '' && visr.param.sample_sheet != '' && visr.param.gene_annotation != ''))")
################################################################

visr.param("enable_filtering", label = "Filtering cells and subsetting genes", default = T,
           info = "Remove outlier cells and genes before further processing.")

visr.param("enable_dim_red", label = "Dimensionality reduction", default = T,
           info = "Reduce dimensionality of data from many genes to fewer number of components.")

visr.param("enable_clustering", label = "Clustering", default = T,
           active.condition = "visr.param.enable_filtering && visr.param.enable_dim_red",
           info = "Identify subtypes of cells using unsupervised clustering.")

visr.param("enable_de_analysis", label = "Differential expression analysis", default = T,
           active.condition = "visr.param.enable_filtering && visr.param.enable_dim_red && visr.param.enable_clustering",
           info = "Characterize differentially expressed genes by comparing groups of cells.")

visr.param("enable_trajectories", label = "Single-cell trajectories", default = T,
           active.condition = "visr.param.enable_filtering && visr.param.enable_dim_red && visr.param.enable_clustering && visr.param.enable_de_analysis",
           info = "Discover cells transition from one state to another.")

visr.param("find_pseudotime_genes", label = "Find pseudotime changing genes", default = T,
           info = "Find genes that change as a function of pseudotime",
           active.condition = "visr.param.enable_filtering && visr.param.enable_dim_red && visr.param.enable_clustering && visr.param.enable_de_analysis && visr.param.enable_trajectories")

visr.param("analyze_trajectory_branches", label = "Analyze branches in trajectories", default = T,
           info = "Analyze branches in single-cell trajectories to identify the genes that differ at a particular branch point",
           active.condition = "visr.param.enable_filtering && visr.param.enable_dim_red && visr.param.enable_clustering && visr.param.enable_de_analysis && visr.param.enable_trajectories")


