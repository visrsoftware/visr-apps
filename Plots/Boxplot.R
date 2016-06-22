source("visrutils.R")
visr.library("car")
visr.applyParameters()

input_table

if (input_ylab=="") {
  input_ylab<-input_var
}

if (input_group != "" & input_xlab == "") {
  input_xlab = "Groups"
}

colors<-gsub("\\s","", strsplit(input_col,",")[[1]])
colors<-gsub("0x","#", colors)


if (input_group == "") {
  mydata <- input_table[,input_var]

  if(input_sort) {
       mydata <- mydata[,order(apply(mydata, 2, FUN = median))]
  }

  boxplot(x = mydata,
          xlab = input_xlab,
          ylab = input_ylab,
          main = input_main,
          range = input_range,
          outline = input_outline,
          notch = input_notch,
          col = colors,
          log = input_log,
          horizontal = input_horizontal,
          varwidth = FALSE,
          las=2,
          pars=list(par(mar=c(input_cexMarginBottom,5,5,5))),
          cex.lab=input_cexAxislabel/10, cex.axis=input_cexAxis/10, cex.main=input_cexTitle/10
    )
}





if (input_group != "" && length(input_var) == 1) {{
  boxplot(input_table[,input_var]~input_table[,input_group],
          xlab = input_xlab,
          ylab = input_ylab,
          main = input_main,
          range = input_range,
          outline = input_outline,
          notch = input_notch,
          col = colors,
          log = input_log,
          horizontal = input_horizontal,
          varwidth = FALSE,
          las=2,
          pars=list(par(mar=c(input_cexMarginBottom,5,5,5))),
          cex.lab=input_cexAxislabel/10, cex.axis=input_cexAxis/10, cex.main=input_cexTitle/10
  )
}}

if(input_group != "" && length(input_var) != 1) {{
  par(mfrow=c(1,length(input_var)))
  for(i in 1:length(input_var)) {
    input_ylab <- input_var[i]
    Boxplot(y = input_table[,input_var[i]],
            g = input_table[,input_group],
            data = input_table,
            xlab = input_xlab,
            ylab = input_ylab,
            main = input_main,
            range = input_range,
            outline = input_outline,
            notch = input_notch,
            col = colors,
            log = input_log,
            horizontal = input_horizontal,
            varwidth = FALSE, id.n = 0,
            las=2,
            pars=list(par(mar=c(input_cexMarginBottom,5,5,5))),
            cex.lab=input_cexAxislabel/10, cex.axis=input_cexAxis/10, cex.main=input_cexTitle/10
    )
  }
}}
