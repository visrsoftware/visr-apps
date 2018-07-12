source("visrutils.R")
visr.app.start("ANOVA",
               info = "One-way ANOVA test: calculates in-group variance and intra-group variance and compares them.",
               debugdata = iris)
visr.category("Model")
visr.param("val", label="response (value)", type="column-numerical", debugvalue = "Sepal.Width")
visr.param("fac", label="explanatory variable (factor)", type="column", debugvalue = "Species")
visr.category("Plot")
visr.param("plot", items = c("anova", "tukey"), item.labels = c("ANOVA-BoxPlot", "Tukey Honest Significant Differences"))
visr.app.end()

visr.applyParameters()

val = visr.input[,visr.param.val]
fac = as.factor(visr.input[,visr.param.fac])

aov1 = aov(val ~ fac)
print(summary(aov1))

tukey <- TukeyHSD(aov1)
print(tukey)

if (visr.param.plot == "anova") {
  boxplot(val ~ fac,
          xlab = visr.param.fac,
          ylab = visr.param.val,
          main = paste("ANOVA Pr(>F) =", format(summary(aov1)[[1]][["Pr(>F)"]][[1]], digits = 2)))
} else if (visr.param.plot == "tukey") {
  par(mar=c(5,10,5,1))
  plot(tukey, las=1)
}
