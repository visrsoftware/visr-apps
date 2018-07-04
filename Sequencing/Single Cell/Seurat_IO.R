
#workflow
visr.app.category("workflow")
visr.param("workflow", label = "Choose Analysis Workflow",
           items = c("single","integrated"),debugvalue = "single",
           item.labels = c("Analysis of one dataset","Integrated Analysis"))

#input
visr.app.category("Input")
visr.param("Import_method", label = "Choose Import Method",
           items = c("load_raw","load_seurat"),
           item.labels = c("Load CellRanger Output","Load Seurat Object"),default = "load_raw", 
           debugvalue = "load_seurat", active.condition = "visr.param.workflow == 'single'")
visr.param("Path_to_outs", type="filename",filename.mode = "dir", 
           info="Cell Ranger pipeline output directory. It should contain another directory named \"outs\"",
           debugvalue = "C:/Users/Yiwei Zhao/Desktop/BRC/tools/meta/single/",
           active.condition = "visr.param.Import_method == 'load_raw' && visr.param.workflow == 'single'")
visr.param("Path_to_seurat_object", type="filename",filename.mode = "load",
           debugvalue = "C:/Users/Yiwei Zhao/Desktop/pbmc.complete.Robj", info = "The object should have been",
           active.condition = "visr.param.Import_method == 'load_seurat' && visr.param.workflow == 'single'")

visr.param("Import_method2", label = "Choose Import Method",
           items = c("load_two","load_one"), default = "load_two", debugvalue = "load_one",
           item.labels = c("Merge two seurat objects for integrated analysis","Load merged seurat object"),
           active.condition = "visr.param.workflow == 'integrated'")
visr.param("seurat_integrated", label = "Merged Seurat Object", type="filename",filename.mode = "load",
           debugvalue = "C:/Users/Yiwei Zhao/Desktop/BRC/tools/SeuratR/Integrated/Aligned.Robj",
           info = "Path to seurat object. It should be a merged object from 2 datasets with different labels (labels should be stored in \"object@meta.data$group\".",
           active.condition = "visr.param.Import_method2 == 'load_one' && visr.param.workflow == 'integrated'")

visr.param("seurat_obj_1", label = "Seurat Object 1", type="filename",filename.mode = "load",
           debugvalue = "C:/Users/Yiwei Zhao/Desktop/BRC/tools/SeuratR/stim_2000.Robj",
           info = "Path to seurat object. The object should have been filtered and contain variable gene information.",
           active.condition = "visr.param.Import_method2 == 'load_two' && visr.param.workflow == 'integrated'")
visr.param("dataset1_name", label = "Label of dataset 1", default = "condition_A",
           active.condition = "visr.param.Import_method2 == 'load_two' && visr.param.workflow == 'integrated'")
visr.param("seurat_obj_2", label = "Seurat Object 2", type="filename",filename.mode = "load",
           debugvalue = "C:/Users/Yiwei Zhao/Desktop/BRC/tools/SeuratR/ctrl_2000.Robj",
           info = "Path to seurat object. The object should have been filtered and contain variable gene information.",
           active.condition = "visr.param.Import_method2 == 'load_two' && visr.param.workflow == 'integrated'")
visr.param("dataset2_name", label = "Label of dataset 2", default = "condition_B",
           active.condition = "visr.param.Import_method2 == 'load_two' && visr.param.workflow == 'integrated'")


load_raw_valid <- "(visr.param.Import_method == 'load_raw' && visr.param.Path_to_outs != '' && visr.param.workflow == 'single')"
load_seurat_valid <- "(visr.param.Import_method == 'load_seurat' && visr.param.Path_to_seurat_object != ''&& visr.param.workflow == 'single')"
load_data_valid <- sprintf("(%s || %s)", load_raw_valid, load_seurat_valid)

load_one_valid <- "(visr.param.workflow == 'integrated' && visr.param.Import_method2 == 'load_one' && visr.param.seurat_integrated != '')"
load_two_valid <- "(visr.param.workflow == 'integrated' && visr.param.Import_method2 == 'load_two' && visr.param.seurat_obj_1 != ''&& visr.param.seurat_obj_2 != '' && visr.param.dataset1_name != ''&& visr.param.dataset2_name != '')"
load_integrated_valid <-  sprintf("(%s || %s)", load_one_valid, load_two_valid)
input_valid <- sprintf("(%s || %s)", load_data_valid, load_integrated_valid)

#output
visr.app.category("Output", active.condition = input_valid)
visr.param("Output_directory", type = "filename", filename.mode="dir",
           label="Output directory to save the results",
           info="Output directory where the analysis results will be saved to",
           debugvalue="C:/Users/Yiwei Zhao/Desktop/temp/")
visr.param("create_subdir", default = T,
           label = "Create new sub-direcory",
           info = "Create a new sub directory with the name DATE_TIME (YYYYMMDD_hhmmss)")


