source("visrutils.R")

if (!visr.isGUI()) {
  visr.input <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"),
           header = T, sep = ";")
}

visr.app.start("UpSet",
               info = "A visualization technique to analyze set-based data. UpSet visualizes set intersections, their properties, and the items (elements) in the dataset.",
               instructions = " - Drag and drop the data table\n - Specify sets to look at\n - Click [ Run ]",
               debugdata = visr.input)

############################################################
visr.category("Input sets")
############################################################
visr.param("sets", "Sets to look at", type = "multi-column",
           debugvalue = c("Drama", "Comedy", "Action", "Thriller", "Western"))


############################################################
visr.category("Layout Options", collapsed = TRUE)
############################################################
visr.param("nintersects", "Number of intersections to plot", default = 40L, min = 1L)
visr.param("order_sets", "Order sets by size", info = "Order sets by size or keep them in the order of input", default = T)
visr.param("group_by", "Group intersections by", info = "How the intersections should be grouped", items = c("degree", "sets"))
visr.param("order_by", "Order intersections by", info = "How the intersections in the matrix should be ordered by",
           items = c("both", "freq", "degree"), item.labels = c("frequency & degree", "frequency", "degree"))
visr.param("decreasing_freq", "Order intersection by decreasing frequency",
           info = "How the intersections should be ordered.",
           default = T,
           active.condition = "visr.param.order_by != 'degree'")
visr.param("decreasing_degree", "Order intersections by decreasing degree",
           info = "How the intersections should be ordered.",
           default = F,
           active.condition = "visr.param.order_by != 'freq'")


############################################################
visr.category("Plot Options", collapsed = TRUE)
############################################################
visr.param("line_size", info = "Width of lines in matrix plot", default = 0.7)
visr.param("point_size", info = "Size of points in matrix plot", default = 2.2)

visr.param("scale_intersections", info = "The scale to be used for the intersection sizes",
           items = c("identity", "log10", "log2"))
visr.param("scale_sets", info = "The scale to be used for the set sizes",
           items = c("identity", "log10", "log2"))


############################################################
visr.category("Plot Labels", collapsed = TRUE)
############################################################
visr.param("main_bar_y_label", info = "The y-axis label of the, intersection size bar plot",
           default = "Intersection Size")
visr.param("sets_x_label", info = "The x-axis label of the set size bar plot",
           default = "Set Size")
visr.param("intersection_size_title", info = "Value to scale the text sizes" ,default = 1)
visr.param("intersection_size_tick_labels", info = "Value to scale the text sizes", default = 1)
visr.param("set_size_title", info = "Value to scale the text sizes", default = 1)
visr.param("set_size_tick_labels", info = "Value to scale the text sizes", default = 1)
visr.param("set_names", info = "Value to scale the text sizes", default = 1)
visr.param("numbers_above_bars", info = "Value to scale the text sizes", default = 1)

visr.param("set_size_angles", info = "Angle to rotate the set size plot x-axis text", default = 0L)
visr.param("show_numbers", info = "Show numbers of intersection sizes above bars", default = T)
visr.param("number_angles", info = "The angle of the numbers atop the intersection size bars", default = 0L,
           active.condition = "visr.param.show_numbers == true")

visr.app.end()

visr.applyParameters()


visr.library("UpSetR")

if (!visr.param.order_sets) {
  visr.param.sets = rev(visr.param.sets) # because plotting order is from bottom to top
}

p <- upset(visr.input,
           sets = visr.param.sets,
           nintersects = visr.param.nintersects,
           keep.order = !visr.param.order_sets,
           group.by = visr.param.group_by,
           order.by = ifelse(visr.param.order_by == "both", c("freq", "degree"), visr.param.order_by),
           decreasing = c(visr.param.decreasing_freq, visr.param.decreasing_degree),
           point.size = visr.param.point_size,
           line.size = visr.param.line_size,
           show.numbers = ifelse(visr.param.show_numbers, "yes", "no"),
           number.angles = visr.param.number_angles,
           set_size.angles = visr.param.set_size_angles,
           mainbar.y.label = visr.param.main_bar_y_label,
           sets.x.label = visr.param.sets_x_label,
           scale.intersections = visr.param.scale_intersections,
           scale.sets = visr.param.scale_sets,
           text.scale = c(visr.param.intersection_size_title,
                          visr.param.intersection_size_tick_labels,
                          visr.param.set_size_title,
                          visr.param.set_size_tick_labels,
                          visr.param.set_names,
                          visr.param.numbers_above_bars))


if (FALSE) {
  upset(nintersects = 40,
        set.metadata = NULL,
        matrix.color = "gray23",
        main.bar.color = "gray23",
        mainbar.y.max = NULL,
        sets.bar.color = "gray23",
        mb.ratio = c(0.7, 0.3),
        expression = NULL, att.pos = NULL,
        att.color = main.bar.color,
        cutoff = NULL, queries = NULL,
        query.legend = "none",
        shade.color = "gray88", shade.alpha = 0.25,
        matrix.dot.alpha = 0.5, empty.intersections = NULL, color.pal = 1,
        boxplot.summary = NULL, attribute.plots = NULL
        )
}
