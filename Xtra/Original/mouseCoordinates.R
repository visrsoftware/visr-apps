# write your R code below
plot(visr.input$start, main="Sample Plot")

x1<-grconvertX(visr.mouse.pressX, from = "device", to = "user")
y1<-grconvertY(visr.mouse.pressY, from = "device", to = "user")
x2<-grconvertX(visr.mouse.currX, from = "device", to = "user")
y2<-grconvertY(visr.mouse.currY, from = "device", to = "user")

rect(min(x1,x2), min(y1,y2), max(x1,x2), max(y1,y2),lty=2)
text(min(x1,x2), min(y1,y2), sprintf("(%4.2f, %4.2f)", min(x1,x2), min(y1,y2)), adj = c(1,1))

#if (x1 != x2 && y1 != y2) {plot(visr.input$start, main="Sample Plot", xlim=c(min(x1,x2),max(x1,x2)), ylim=c(min(y1,y2),max(y1,y2)))}
