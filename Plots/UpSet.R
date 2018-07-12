source("visrutils.R")

if (!visr.isGUI()) {# debug only
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
visr.param("input_type", "Input Type", items = c("binary", "id_set"), item.labels = c("binary encoding for sets", "(id, set) combination"))

visr.param("sets", "Sets to look at",
           info = "Select the columns representing the sets. All rows should either have a value of 1 or 0 for each selected columns \n indicating whether or not the row belongs to the corresponding set.",
           type = "multi-column",
           active.condition = "visr.param.input_type == 'binary'",
           debugvalue = c("Drama", "Comedy", "Action", "Thriller", "Western"))

visr.param("input_id_column", "id column",
           info = "The column with the data item IDs",
           type = "column", active.condition = "visr.param.input_type == 'id_set'")

visr.param("input_set_column", "set column",
           info = "The column with the data item set names",
           type = "column", active.condition = "visr.param.input_type == 'id_set'")

visr.param("output_binary", "Also output the binary encoded table", default = FALSE,
           active.condition = "visr.param.input_type == 'id_set'")

visr.param("output_binary_table", "Name for output table", type = "output-table",
           active.condition = "visr.param.input_type == 'id_set' && visr.param.output_binary")


############################################################
visr.category("Layout Options", collapsed = F)
############################################################
visr.param("nintersects", "Maximum number of intersections to show", type = 'integer',
           default = 40L, min = 1L, items = c('NA'), item.labels = c('All'))
visr.param("empty.intersections", "Display empty sets",
           info = "Additionally display empty sets up to specified maximum number of intersections",
           default = F, debugvalue = T)
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
visr.category("Axis Options", collapsed = T)
############################################################

visr.param("mainbar.y.max", "Y-axis max for intersection sizes", info = "The maximum y value of the intersection size bar plot scale. May be useful when aligning multiple UpSet plots horizontally.",
           type = "double", default = "NULL", items = c('NULL'), item.labels = c('Auto'), debugvalue = NULL)

visr.param("scale_intersections", "Y-axis scaling for intersection sizes", info = "The scale to be used for the intersection sizes",
           items = c("identity", "log10", "log2"))

visr.param("main_bar_y_label", 'Y-axis label', info = "The y-axis label of the, intersection size bar plot",
           default = "Intersection Size")

visr.param("scale_sets", "X-axis scaling for set sizes", info = "The scale to be used for the set sizes",
           items = c("identity", "log10", "log2"))

visr.param("sets_x_label", 'X-axis label', info = "The x-axis label of the set size bar plot",
           default = "Set Size")

############################################################
visr.category("Size Options", collapsed = T)
############################################################
visr.param("line_size", info = "Width of lines in matrix plot", default = 0.7)

visr.param("point_size", info = "Size of points in matrix plot", default = 2.2)

visr.param("matrix.ratio", "matrix plot size ratio", info = 'Ratio between matrix plot (intersections) and main bar plot', default = 0.7, min = 0.01, max = 0.99)

visr.param("intersection_size_title", info = "Value to scale the text sizes" ,default = 1)
visr.param("intersection_size_tick_labels", info = "Value to scale the text sizes", default = 1)
visr.param("set_size_title", info = "Value to scale the text sizes", default = 1)
visr.param("set_size_tick_labels", info = "Value to scale the text sizes", default = 1)
visr.param("set_names", info = "Value to scale the text sizes", default = 1)
visr.param("numbers_above_bars", info = "Value to scale the text sizes", default = 1)

visr.param("set_size_angles", info = "Angle to rotate the set size plot x-axis text", default = 0L)
visr.param("show_numbers", info = "Show numbers of intersection sizes above bars", default = T)
visr.param("number_angles", info = "The angle of the numbers atop the intersection size bars", default = 0L,
           active.condition = "visr.param.show_numbers")

############################################################
visr.category("Color Options", collapsed = TRUE)
############################################################
visr.param("main.bar.color", "Color of the main bar plot", type="color", default = "gray23")
visr.param("sets.bar.color", "Color of set size bar plot", type="color", default = "gray23")

visr.param("matrix.color", "Color of intersection points", type="color", default = "gray23")
visr.param("matrix.dot.alpha", "Transparency of empty intersection points", default = 0.5, min = 0, max = 1)

visr.param("shade.color", "Color of row shading in matrix", type="color", default = "gray88")
visr.param("shade.alpha", "Transparency of shading in matrix", default = 0.25, min = 0, max = 1)

############################################################
visr.category("Box plots", info = "Boxplots representing the distribution of a selected attribute for each intersection.", collapsed = TRUE)
############################################################
visr.param("boxplot1", type = "column-numerical",
           info = "Column used for boxplots representing the distribution of the selected column for each intersection.",
           debugvalue = "AvgRating")

visr.param("boxplot2", type = "column-numerical",
           info = "Column used for boxplots representing the distribution of the selected column for each intersection.",
           debugvalue = "ReleaseDate")

visr.app.end()

visr.applyParameters()


visr.library("UpSetR")

if (!visr.param.order_sets) {
  visr.param.sets = rev(visr.param.sets) # because plotting order is from bottom to top
}

boxplot.summary = c (if (visr.param.boxplot1 != '') visr.param.boxplot1 else NULL,
                     if (visr.param.boxplot2 != '') visr.param.boxplot2 else NULL)

if (!is.null(boxplot.summary))
  visr.param.empty.intersections = F # currently there is a bug when both are specified

p <- upset(visr.input,
           sets = visr.param.sets,
           nintersects = visr.param.nintersects,
           empty.intersections = if (visr.param.empty.intersections) 'on' else NULL,
           keep.order = !visr.param.order_sets,
           group.by = visr.param.group_by,
           order.by = if (visr.param.order_by == "both") c("freq", "degree") else visr.param.order_by,
           decreasing = c(visr.param.decreasing_freq, visr.param.decreasing_degree),
           point.size = visr.param.point_size,
           line.size = visr.param.line_size,
           show.numbers = if (visr.param.show_numbers) "yes" else "no",
           number.angles = visr.param.number_angles,
           set_size.angles = visr.param.set_size_angles,
           mainbar.y.max = visr.param.mainbar.y.max,
           mainbar.y.label = visr.param.main_bar_y_label,
           sets.x.label = visr.param.sets_x_label,
           scale.intersections = visr.param.scale_intersections,
           scale.sets = visr.param.scale_sets,
           text.scale = c(visr.param.intersection_size_title,
                          visr.param.intersection_size_tick_labels,
                          visr.param.set_size_title,
                          visr.param.set_size_tick_labels,
                          visr.param.set_names,
                          visr.param.numbers_above_bars),
           matrix.color = visr.param.matrix.color,
           matrix.dot.alpha = visr.param.matrix.dot.alpha,
           main.bar.color = visr.param.main.bar.color,
           sets.bar.color = visr.param.sets.bar.color,
           shade.color = visr.param.shade.color,
           shade.alpha = visr.param.shade.alpha,
           # queries = list(list(query = intersects, params = list("Action", "Drama"), active = T)),
           # attribute.plots = list(gridrows = 50, plots = list(list(plot = histogram, x = "ReleaseDate", queries = F))),
           boxplot.summary = boxplot.summary,
           mb.ratio = c(visr.param.matrix.ratio, 1.0 - visr.param.matrix.ratio)
)


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
