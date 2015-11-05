usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org")}}

usePackage("ggplot2")
usePackage("reshape2")
usePackage("seriation")

input_columns <- c("index","start","end")
sort_column<-"index"
sort_decreasing<-FALSE
input_normalize<-TRUE

plot_title<-"heatmap"
text_size<-20
color_low<-"blue"
color_mid<-"white"
color_high<-"red"

visr.applyParameters()

heatmap_col2<-colorRampPalette(c(color_low, color_mid, color_high))(n = 1000)

d<-subset(visr.input, select = input_columns)
d<-d[order(visr.input[,sort_column], decreasing=sort_decreasing),]
m<-data.matrix(d)
if (input_normalize) {m <- sweep(m, 2, apply(m, 2, min))}
if (input_normalize) {m <- sweep(m, 2, apply(m, 2, max), "/")}

#m<-apply(m, 2, function(x) min((x-min(x))/(2*median(x)-min(x)),1.0))

#heatmap.2(m, cexCol=text_size/10, col=heatmap_col, main=plot_title, key=TRUE, dendrogram="none", Colv = "NA", Rowv = "NA", density.info="none", trace="none", labRow = FALSE, labCol=input_columns,lmat=rbind(c(2),c(3),c(1),c(4)),lhei=c(2,3,9,1), lwid=c(1))

# http://is-r.tumblr.com/post/32449990608/optimal-seriation-for-your-matrices
longData <- melt(m)
head(longData)
zp1 <- ggplot(longData, aes(x = Var2, y = Var1, fill = value))
zp1 <- zp1 + geom_tile()
zp1 <- zp1 + scale_fill_gradientn(colours = heatmap_col2)
zp1 <- zp1 + scale_x_discrete(expand = c(0, 0))
#zp1 <- zp1 + scale_y_discrete(expand = c(0, 0))
zp1 <- zp1 + theme(axis.text.x=element_text(angle=45, hjust = 1, size = text_size))
zp1 <- zp1 + theme(axis.text.y = element_blank())
zp1 <- zp1 + theme(axis.title.x = element_blank())
zp1 <- zp1 + theme(axis.title.y = element_blank())
zp1 <- zp1 + theme(axis.ticks=element_blank())
print(zp1)