source("visrutils.R")
visr.library("ggplot2")
visr.library("scales")
visr.library("grid")

visr.applyParameters()

input_legend_title <- ifelse(input_legend_title=="", input_color, input_legend_title)
{{
  if(input_color != "") {
    input_color_modified <- ifelse(input_color_factor, paste("as.factor(",input_color,")",sep=""), input_color)
    p <- ggplot(visr.input, aes_string(x=input_x, y=input_y, colour=input_color_modified)) #+ guides(colour=guide_legend(title="New Legend Title")) #+ scale_colour(name = input_color)
  } else {  
    p <- ggplot(visr.input, aes_string(x=input_x, y=input_y)) 
  }
}}


p <- p + geom_point(alpha = input_point_alpha, size = input_point_size) 
p <- p + theme(text = element_text(size=input_text_size))
p <- p + theme(legend.position = input_legend_position, legend.key = element_rect(colour = "black"))
p <- p + theme(legend.margin = unit(-0.5, "cm"))
{{
  if (input_color != "") {
    if (input_color_factor) {
      p <- p + scale_color_discrete(name=input_legend_title)
    } else {
      p <- p + scale_color_continuous(name=input_legend_title)
    }
  }
}}

{{
if (input_xscale == "log2") { 
  p <- p + scale_x_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) 
} else if (input_xscale == "log10") {  
  p <- p + scale_x_continuous(trans = log10_trans(), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) 
} else {  
  p <- p + scale_x_continuous(label=get(input_xscale)) 
}
}}

{{
if (input_yscale == "log2") {
  p <- p + scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) 
} else if (input_yscale == "log10") { 
  p <- p + scale_y_continuous(trans = log10_trans(), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) 
} else {  
  p <- p + scale_y_continuous(label=get(input_yscale))
}
}}

{{
if (input_label != "") {
  p <- p + geom_text(aes_string(input_x,input_y, label = input_label), hjust = 0.4, vjust = -0.5)
}
}}

p <- p + ggtitle(input_title)

print(p)