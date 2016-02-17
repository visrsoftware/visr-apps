#decent heatmaps in R: http://sebastianraschka.com/Articles/heatmaps_in_r.html
usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}

usePackage("gplots")

input_columns <- c("index","start","end")
input_sort_column<-"index"
input_sort_decreasing<-FALSE
input_normalize<-TRUE

input_plot_title<-"heatmap"
input_text_size<-20
input_color_low<-"blue"
input_color_mid<-"white"
input_color_high<-"red"

input_dendogram<-"none"
input_colRotate<-25
input_rowRotate<-45
input_hasKey<-TRUE
input_keySize<-1
input_rowLabels<-TRUE
input_colLabels<-TRUE
input_rowOrdered<-FALSE



visr.applyParameters()

# workaround to allow heatmap to be drawn when there is a single column
if (length(input_columns) == 1) {input_columns = c(input_columns, input_columns)}

heatmap_col<-colorRampPalette(c(input_color_low, input_color_mid, input_color_high))(n = 1000)

d<-subset(visr.input, select = input_columns)
#d<-scale(d)
if (input_sort_column != "") d<-d[order(visr.input[,input_sort_column], decreasing=input_sort_decreasing),]

m<-data.matrix(d)
fnorm<-function(x) {return (pmin( (x-min(x))/(2*median(x)-min(x)), 1.0))}
if (input_normalize) {m<-apply(m, 2, fnorm)}
#if (input_normalize) {m <- sweep(m, 2, apply(m, 2, min))}
#if (input_normalize) {m <- sweep(m, 2, apply(m, 2, max), "/")}

if (input_colLabels) input_colLabels<- input_columns else input_colLabels <- NULL
if (input_rowLabels) input_rowLabels<- rownames(visr.input) else input_rowLabels <- NULL


#m<-apply(d, 2, sd)
{{
heatmap.2(m, cexCol=input_text_size/10, col=heatmap_col, main=input_plot_title, 
          dendrogram=input_dendogram, Colv = "NA", Rowv = input_rowOrdered, 
          density.info="none", trace="none", 
          labRow = input_rowLabels, labCol=input_colLabels,
          key=input_hasKey, keysize = input_keySize, 
          srtCol=input_colRotate, srtRow=input_rowRotate,
          lmat=rbind(c(2),c(3),c(1),c(4)),lhei=c(1,1,10,1), lwid=c(1))
}}