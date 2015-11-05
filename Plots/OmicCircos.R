biocPackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {source("http://bioconductor.org/biocLite.R"); biocLite(pkg);require(pkg, character.only = TRUE)}}
biocPackage("GenomicRanges")

source(paste(getwd(),"/Graph/OmicCircosLib.R",sep = ""))

input_hlcolshade <- "#CCFF001A"
input_hlcolborder <- "#FF000080"

visr.applyParameters()

if (input_name == "") error_message<-"No chromosome column is specified.\n[ Parameters -> Segments -> Name ]"
if (input_start == "") error_message<-"No chromosome start position is specified\n[ Parameters -> Segments -> start ]"
if (input_end == "") error_message<-"No chromosome end position is specified.\n[ Parameters -> Segments -> end ]"
if (input_position == "") input_position <- input_start; error_message<-"Warning: chromosome start position is used as chromosome position."

if (input_cutoff == 0) input_cutoff = "n"

options(stringsAsFactors = FALSE)
set.seed(1234)

seg.f <- cbind(input_table[,c(input_name,input_start,input_end)],"NA","NA")

seg.num <- length(unique(seg.f[,1]))


seg.name <- unique(seg.f[,1])
db <- segAnglePo(seg.f, seg=seg.name)


colors <- rainbow(seg.num, alpha=input_alpha)

if(par()$fin[1] <= par()$fin[2]) {winx <- 800; winy <- 800*par()$fin[2]/par()$fin[1]}
if(par()$fin[1] > par()$fin[2]) {winx <- 800*par()$fin[1]/par()$fin[2]; winy <- 800}
par(mar=c(0, 0, 0, 0)) # set margins
plot(c(1,winx), c(1,winy), type="n", axes=FALSE, xlab="", ylab="", main="",cex=0.1)

# outmost anchor track
circos(R=input_R0, print.chr.lab=input_printchrlab, scale=input_scale0, type="chr", cir=db, col=colors, W=4, xc = 400,yc = 400,cex=0.1)

#values
seg.v <- input_table[,c(input_name,input_position,input_values)]

{{
  circos(mapping = seg.v,
         R = input_R0*input_R/100,
         W = input_R0*input_w/100,
         cir = db,
         type = input_type,
         col.v = 3,
         B = input_B,
         col.bar=input_colbar, 
         col.bar.po = input_colbarpo,
         cluster = input_cluster,
         scale = input_scale,
         cutoff = input_cutoff,
         #zoom="",
         lwd = input_lwd,
         cex = 1,
         col=colors[c(1,2)]) 
}}


#### highlight "hl"

{{
  if(input_highlight == TRUE) {
    highlight <- c(input_R0*input_hlinnerradius/100, input_R0*input_hlouterradius/100, input_hlstartsegment, input_hlstartposition, input_hlendsegment, input_hlendposition, input_hlcolshade, input_hlcolborder)
    circos(R=input_R0, cir=db, W=input_R0/2, mapping=highlight, type="hl", lwd=2)
  }
}}




