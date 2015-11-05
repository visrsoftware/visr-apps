pca_input <- c("index")
pca_cor<-TRUE

visr.applyParameters()

pca_out<-princomp(subset(visr.input, select = pca_input), cor = pca_cor)
plot(pca_out, main="Principal Components")
