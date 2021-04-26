#' zchoo Monday, Apr 26, 2021 02:24:31 PM
## This is a script to generate .png files for SV gallery

#' @name ridge.plot
#' @title ridge.plot
#'
#' @description
#'
#' Generates ridge plot of event burden against background distribution of event counts
#'
#' @param complex.fname (character) output from complex event caller
#' @param background.fname (character) text file with event burdens (default sv.burden.txt)
#' @param out.fname (character) ridge plot output file name (default ridgeplot.png)
ridge.plot = function(complex.fname = NULL,
                      background.fname = "/data/sv.burden.txt",
                      out.fname = "ridgeplot.png") {
    if (!file.exists(complex.fname)) {
        stop("complex.fname does not exist")
    }
    if (!file.exists(background.fname)) {
        stop("background.fname does not exist")
    }

    ## read complex output
    this.complex = readRDS(complex.fname)

    ## read background distribution of events and extract event types
    sv.burden.dt = fread(background.fname)
    ev.types = unique(sv.burden.dt$type)

    if (is.null(this.complex$meta$events) || nrow(this.complex$meta$events) == 0) {
        warning("No detected SVs! Plotting background distribution only.")
        this.burden.dt = data.table(type = ev.types, burden = 0)
    } else {
        ev.counts = table(this.complex$meta$events$type)
        this.burden.dt = data.table(type = ev.types,
                                    burden = ifelse(ev.types %in% names(ev.counts), ev.counts[ev.types], 0))
    }

    ## synchronize factor levels
    sv.burden.dt[, type := factor(type, levels = ev.types)]
    this.burden.dt[, type := factor(type, levels = ev.types)]

    ## prepare plot
    pt = ggplot(sv.burden.dt[type != "unclassified",], aes(x = burden, y = type, fill = type)) +
        geom_density_ridges(bandwidth = 0.1, alpha = 0.5, scale = 0.9,
                            rel_min_height = 0.01, color = "white",
                            jittered_points = TRUE,
                            position = position_points_jitter(width = 0.01, height = 0),
                            point_shape = '|', point_size = 3, point_alpha = 0.1, point_colour = "black") +
        geom_segment(data = this.burden.dt[ type != "unclassified"],
                     aes(x = burden, xend = burden, y = as.numeric(type), yend = as.numeric(type) + 0.9),
                     color = "red") +
        scale_x_continuous(trans = "log1p", breaks = c(0, 1, 10, 100)) +
        labs(x = "Event Burden", y = "") +
        theme_ridges(center = TRUE) +
        theme(legend.position = "none",
              axis.title = element_text(size = 20, family = "sans"),
              axis.text.x = element_text(size = 15, family = "sans"),
              axis.text.y = element_text(size = 20, family = "sans"))

    ppng(print(pt), filename = out.fname)
}
