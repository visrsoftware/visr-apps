source("visrutils.R")
visr.app.start("Independence",
               info = "Performs test of independence using different tests.",
               debugdata = mtcars)
visr.category("Model")
visr.param("vars1", label="Variables1", type="multi-column", debugvalue = c("gear", "mpg"))
visr.param("vars2", label="Variables2", type="multi-column", debugvalue = c("cyl", "disp"))
visr.category("Output")
visr.param("output", label="P-Value table", type = "output-table")
visr.app.end()

visr.applyParameters()

isCategorical <- function(x) {
  if (is.numeric(x)) {
    if (length(levels(as.factor(x))) > 10) {
      print(length(levels(as.factor(x))))
      return(F)
    }
  }
  return(T)
}

visr.param.output <- NULL

for (v1 in visr.param.vars1) {
  for (v2 in visr.param.vars2) {
    var1 <- visr.input[, v1]
    var2 <- visr.input[, v2]
    var1IsCat <- isCategorical(var1)
    var2IsCat <- isCategorical(var2)

    if (var1IsCat && var2IsCat) {
      # two categorical
      # chi2 <- chisq.test(table(visr.input[, c(v1, v2)]))
      chi2 <- chisq.test(table(data.frame(cbind(as.factor(var1)), as.factor(var2))))
      pvalue <- chi2$p.value
      info <- "chisq.test"
    } else if (!var1IsCat && !var2IsCat) {
      # two numeric
      ktest <- kruskal.test(var1 ~ var2)
      pvalue <- ktest$p.value
      info <- "kruskal.test"
    } else {
      # one categorical and one numerical
      aov1 <- aov(var1 ~ var2)
      pvalue <- summary(aov1)[[1]][["Pr(>F)"]][[1]]
      info <- "aov"
      # boxplot(var1 ~ var2)
    }
    visr.param.output <- rbind(visr.param.output, c(v1, v2, pvalue, info))
    print(paste(v1, v2, pvalue))
  }
}
colnames(visr.param.output) <- c("var1", "var2", "pvalue", "info")

