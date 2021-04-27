#' zchoo Monday, Apr 26, 2021 02:24:31 PM
## This is a script to generate .png files for SV gallery

#' @name gallery.wrapper
#' @title gallery.wrapper
#'
#' @description
#'
#' wrapper function to prep input for SlickR
#' @param complex.fname (character)
#' @param cvgt.fname
#' @param background.fname
#' @param server
#' @param pair
#' @param ev.types
#' @param pad
#' @param height
#' @param width
#' @param outdir
#' 
#' @return data.table with columns from complex and additionally plot.fname and plot.link
gallery.wrapper = function(complex.fname = NULL,
                           background.fname = "/data/sv.burden.txt",
                           cvgt.fname = "./coverage.gtrack.rds",
                           server = "",
                           pair = "",
                           ev.types = c("qrp", "tic", "qpdup", "qrdel",
                                        "bfb", "dm", "chromoplexy", "chromothripsis",
                                        "tyfonas", "rigma", "pyrgo"),
                           pad = 5e5,
                           height = 1000,
                           width = 1000,
                           outdir = "./") {
    ## generate ridge plot
    ridgeplot.fname = ridge.plot(complex.fname = complex.fname,
                                 background.fname = background.fname,
                                 ev.types = ev.types,
                                 height = height,
                                 width = width,
                                 outdir = outdir)

    svplot.dt = sv.plot(complex.fname = complex.fname,
                        cvgt.fname = cvgt.fname,
                        server = server,
                        pair = pair,
                        ev.types = ev.types,
                        pad = pad,
                        height = height,
                        width = width,
                        outdir = outdir)

    out = rbind(data.table(plot.fname = ridgeplot.fname),
                svplot.dt,
                fill = TRUE,
                use.names = TRUE)
    
    return (out)
}
                           
#' @name sv.plot
#' @title sv.plot
#'
#' @description
#'
#' Creates .png files corresponding to gTrack plots for each complex SV pattern
#'
#' @param complex.fname (character) output from complex event caller
#' @param cvgt.fname (character) coverage gTrack
#' @param server (character) server url
#' @param pair (character) pair id
#' @param ev.types (character) complex event types
#' @param pad (numeric) window padding in bp default 5e5
#' @param height (numeric) plot height default 500
#' @param width (numeric) plot width default 500
#' @param outdir (character) where to save plots?
#'
#' @return data.table with columns id, plot.fname, and plot.link
sv.plot = function(complex.fname = NULL,
                   cvgt.fname = "./coverage.gtrack.rds",
                   server = "",
                   pair = "",
                   ev.types = c("qrp", "tic", "qpdup", "qrdel",
                                "bfb", "dm", "chromoplexy", "chromothripsis",
                                "tyfonas", "rigma", "pyrgo"),
                   pad = 5e5,
                   height = 500,
                   width = 500,
                   outdir = "./") {
    if (!file.exists(complex.fname)) {
        stop("complex.fname does not exist")
    }
    if (!file.exists(cvgt.fname)) {
        stop("coverage gTrack does not exist")
    }

    ## grab complex and coverage gTracks
    this.complex = readRDS(complex.fname)
    this.complex.gt = this.complex$gt
    cvgt = readRDS(cvgt.fname)

    ## extract just complex events
    complex.ev = this.complex$meta$events[type %in% ev.types,]

    ## grab plot titles and file names
    complex.ev[, plot.fname := file.path(outdir, paste("event", ev.id, "png", sep = "."))]
    complex.ev[, plot.title := paste(type, "|", "event id:", ev.id)]

    ## make gGnome.js url query parameters
    complex.ev[, ev.js.range := {
        gr = gr.reduce(gr.stripstrand(parse.gr(footprint)) + 1e4)
        paste(gr.string(gr), collapse = "%20|%20")
    }, by = ev.id]
    complex.ev[, plot.link := paste0(server, "index.html?file=", pair, ".json&location=", ev.js.range, "&view=")]

    ## save plots
    pts = lapply(1:nrow(complex.ev),
                 function(ix) {
                     ppng(plot(c(cvgt, this.complex.gt), parse.grl(complex.ev$footprint[ix]) + pad),
                          title = complex.ev$plot.title[ix],
                          filename = complex.ev$plot.fname[ix],
                          height = height,
                          width = width)
                 })

    return(complex.ev[, .(plot.fname, plot.link)])
}
                   

#' @name ridge.plot
#' @title ridge.plot
#'
#' @description
#'
#' Generates ridge plot of event burden against background distribution of event counts
#'
#' @param complex.fname (character) output from complex event caller
#' @param background.fname (character) text file with event burdens (default sv.burden.txt)
#' @param ev.types (character)
#' @param height (numeric) height of png default 1e3
#' @param width (numeric) width of png default 1e3
#' @param outdir (character) ridge plot output directory 
ridge.plot = function(complex.fname = NULL,
                      background.fname = "/data/sv.burden.txt",
                      ev.types = c("qrp", "tic", "qpdup", "qrdel",
                                   "bfb", "dm", "chromoplexy", "chromothripsis",
                                   "tyfonas", "rigma", "pyrgo"),
                      height = 1000,
                      width = 1000,
                      outdir = "./") {
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
    ev.types = intersect(unique(sv.burden.dt$type), ev.types)

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

    ## save plot
    out.fname = file.path(outdir, "ridgeplot.png")
    ppng(print(pt), filename = out.fname, height = height, width = width)
    return(out.fname)
}
