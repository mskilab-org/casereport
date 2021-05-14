#' zchoo Monday, Apr 26, 2021 02:24:31 PM
## This is a script to generate .png files for SV gallery

#' @name cgc.gtrack
#'
#' @param cgc.fname (character)
#' @param gencode (character)
#' @param outdir (character)
#' @param verbose (logical)
#'
#' @return file name of gencode gTrack
cgc.gtrack = function(cgc.fname = "./data/cgc.tsv",
                      gencode.fname = NULL,
                      outdir = "./") {
    
    if (!file.exists(cgc.fname)) {
        stop("cgc.fname does not exist")
    }
    if (is.null(gencode.fname) || !file.exists(gencode.fname)) {
        stop("gencode does not exist")
    }

    cgc.gene.symbols = fread(cgc.fname)[["Gene Symbol"]]
    gff = skidb::read_gencode(fn = gencode.fname)
    cgc.gt = track.gencode(gencode = gff, genes = cgc.gene.symbols)

    fn = file.path(outdir, "cgc.gtrack.rds")
    saveRDS(cgc.gt, fn)
    return(fn)
}

#' @name gallery.wrapper
#' @title gallery.wrapper
#'
#' @description
#'
#' wrapper function to prep input for SlickR
#' @param complex.fname (character)
#' @param cvgt.fname
#' @param gngt.fname
#' @param cgcgt.fname CGC gene gTrack
#' @param agt.fname allele gTrack path
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
                           gngt.fname = "./data/gt.ge.hg19.rds",
                           cgcgt.fname = NULL,
                           agt.fname = NULL,
                           server = "",
                           pair = "",
                           ev.types = c("qrp", "tic", "qpdup", "qrdel",
                                        "bfb", "dm", "chromoplexy", "chromothripsis",
                                        "tyfonas", "rigma", "pyrgo", "cpxdm"),
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
                        gngt.fname = gngt.fname,
                        cgcgt.fname = cgcgt.fname,
                        agt.fname = agt.fname,
                        server = server,
                        pair = pair,
                        ev.types = ev.types,
                        pad = pad,
                        height = height,
                        width = width,
                        outdir = outdir)

    ## out = rbind(data.table(plot.fname = ridgeplot.fname),
    ##             svplot.dt,
    ##             fill = TRUE,
    ##             use.names = TRUE)
    
    ## return (out)
    return(svplot.dt)
}

#' @name cn.plot
#' @title cn.plot
#'
#' @description
#'
#' Creates .png files corresponding to gTrack plots for each amplified/deleted gene
#'
#' @param drivers.fname text file with drivers + CN alterations + genomic locations
#' @param complex.fname (character) output from complex event caller
#' @param cvgt.fname (character) coverage gTrack file path
#' @param gngt.fname (character) gencode gTrack file path
#' @param cgcgt.fname (character) CGC gTrack file path
#' @param agt.fname (character) allele gTrack file name
#' @param server (character) server url
#' @param pair (character) pair id
#' @param ev.types (character) complex event types
#' @param pad (numeric) window padding in bp default 5e5
#' @param height (numeric) plot height default 500
#' @param width (numeric) plot width default 500
#' @param outdir (character) where to save plots?
#'
#' @return data.table with columns id, plot.fname, and plot.link
cn.plot = function(drivers.fname = NULL,
                   complex.fname = NULL,
                   cvgt.fname = "./coverage.gtrack.rds",
                   gngt.fname = "./data/gt.ge.hg19.rds",
                   cgcgt.fname = NULL,
                   agt.fname = NULL,
                   server = "",
                   pair = "",
                   ev.types = c("qrp", "tic", "qpdup", "qrdel",
                                "bfb", "dm", "chromoplexy", "chromothripsis",
                                "tyfonas", "rigma", "pyrgo", "cpxdm"),
                   pad = 5e5,
                   height = 1000,
                   width = 1000,
                   outdir = "./") {
    if (!file.exists(drivers.fname)) {
        stop("drivers.fname does not exist")
    }
    if (!file.exists(complex.fname)) {
        stop("complex.fname does not exist")
    }
    if (!file.exists(cvgt.fname)) {
        stop("coverage gTrack does not exist")
    }

    ## grab drivers and convert to GRanges
    if (grepl("rds", drivers.fname)) {
        drivers.dt = readRDS(drivers.fname) %>% as.data.table
    } else {
        drivers.dt = fread(drivers.fname) %>% as.data.table
    }

    ## set gene names
    if (!("gene_name" %in% colnames(drivers.dt)) & ("gene" %in% colnames(drivers.dt))) {
        drivers.dt[, gene_name := gene]
    }

    message(nrow(drivers.dt))

    ## if empty data table, return
    if (nrow(drivers.dt) > 0) {

        ## grab complex and coverage gTracks
        this.complex = readRDS(complex.fname)
        this.complex.gt = this.complex$gt
        cvgt = readRDS(cvgt.fname) ## coverage
        gngt = readRDS(gngt.fname) ## gencode gTrack


        if ("seqnames" %in% colnames(drivers.dt)) {
            drivers.gr = dt2gr(drivers.dt)
        } else {

            ## otherwise grab gr from gencode data and intersect with things
            ge.stacked = gngt@data[[1]] %>% range %>% stack
            drivers.ge.dt = gr2dt(ge.stacked %Q% (name %in% drivers.dt$gene_name))
            drivers.ge.dt[, gene_name := name]
            drivers.gr = dt2gr(drivers.ge.dt)
            drivers.dt = gr2dt(drivers.gr)
            
        }

        ## grab plot titles and file names
        drivers.dt[, plot.fname := file.path(outdir, paste(gene_name, "png", sep = "."))]

        ## make gGnome.js url query parameters
        drivers.dt[, ev.js.range := gr.string(drivers.gr)]
        drivers.dt[, plot.link := paste0(server, "index.html?file=", pair, ".json&location=", ev.js.range, "&view=")]

        ## read CGC gTrack if provided
        if (!is.null(cgcgt.fname)) {
            if (file.exists(cgcgt.fname)) {
                cgc.gt = readRDS(cgcgt.fname)
            } else {
                cgc.gt = NULL
            }
        } else {
            cgc.gt = NULL
        }

        ## read allele gTrack if provided
        if (!is.null(agt.fname)) {
            if (file.exists(agt.fname)) {
                agt = readRDS(agt.fname)
            } else {
                agt = NULL
            }
        } else {
            agt = NULL
        }

        ## gTrack formatting
        cvgt$ylab = "CN"
        cvgt$name = "cov"
        cvgt$yaxis.pretty = 3
        cvgt$xaxis.chronly = TRUE

        this.complex.gt$ylab = "CN"
        this.complex.gt$name = "JaBbA"
        this.complex.gt$yaxis.pretty = 3
        this.complex.gt$xaxis.chronly = TRUE

        gngt$cex.label = 0.01 ## no labels for full gencode
        gngt$xaxis.chronly = TRUE
        gngt$name = "genes"
        gngt$ywid = 0.1
        gngt$height = 2
        gngt$yaxis.cex = 0.8

        if (!is.null(cgc.gt)) {
            cgc.gt$cex.label = 0.5
            cgc.gt$xaxis.chronly = TRUE
            cgc.gt$name = "CGC"
            cgc.gt$ywid = 0.1
            cgc.gt$height = 5
            cgc.gt$stack.gap = 0.5
            gngt$yaxis.cex = 0.8

            ## concatenate final gTracks
            gt = c(gngt, cgc.gt, cvgt, this.complex.gt)
        } else {
            gt = c(gngt, cvgt, this.complex.gt)
        }

        if (!is.null(agt)) {
            agt$ylab = "CN"
            agt$yaxis.pretty = 3
            agt$xaxis.chronly = TRUE

            gt = c(agt, gt)
        }

        ## make one plot per range in drivers gr
        pts = lapply(1:length(drivers.gr),
                     function(ix) {
                         ## prepare window
                         win = drivers.gr[ix]
                         if (pad > 0 & pad <= 1) {
                             adjust = pmax(5e5, pad * width(win))
                             win = GenomicRanges::trim(win + adjust)
                         } else {
                             win = GenomicRanges::trim(win + pad)
                         }
                         ppng(plot(gt, win, legend.params = list(plot = FALSE)),
                              title = drivers.gr$gene_name[ix], ## title is the gene name
                              filename = drivers.dt$plot.fname[ix],
                              height = height,
                              width = width)
                     })

        return(drivers.dt[, .(gene_name, plot.fname, plot.link)])
    } else {
        warning("No drivers with CN change!")
        return(data.table(gene_name = character(), plot.fname = character(), plot.link = character()))
    }
}

#' @name sv.plot
#' @title sv.plot
#'
#' @description
#'
#' Creates .png files corresponding to gTrack plots for each complex SV pattern
#'
#' @param complex.fname (character) output from complex event caller
#' @param cvgt.fname (character) coverage gTrack file path
#' @param gngt.fname (character) gencode gTrack file path
#' @param cgcgt.fname (character) CGC gTrack file path
#' @param agt.fname (character) allele gTrack file name
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
                   gngt.fname = "./data/gt.ge.hg19.rds",
                   cgcgt.fname = NULL,
                   agt.fname = NULL,
                   server = "",
                   pair = "",
                   ev.types = c("qrp", "tic", "qpdup", "qrdel",
                                "bfb", "dm", "chromoplexy", "chromothripsis",
                                "tyfonas", "rigma", "pyrgo", "cpxdm"),
                   pad = 5e5,
                   height = 1000,
                   width = 1000,
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
    cvgt = readRDS(cvgt.fname) ## coverage
    gngt = readRDS(gngt.fname) ## gencode gTrack

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

    ## read CGC gTrack if provided
    if (!is.null(cgcgt.fname)) {
        if (file.exists(cgcgt.fname)) {
            cgc.gt = readRDS(cgcgt.fname)
        } else {
            cgc.gt = NULL
        }
    } else {
        cgc.gt = NULL
    }

    ## read allele gTrack if provided
    if (!is.null(agt.fname)) {
        if (file.exists(agt.fname)) {
            agt = readRDS(agt.fname)
        } else {
            agt = NULL
        }
    } else {
        agt = NULL
    }

    ## gTrack formatting
    cvgt$ylab = "CN"
    cvgt$name = "cov"
    cvgt$yaxis.pretty = 3
    cvgt$xaxis.chronly = TRUE

    this.complex.gt$ylab = "CN"
    this.complex.gt$name = "JaBbA"
    this.complex.gt$yaxis.pretty = 3
    this.complex.gt$xaxis.chronly = TRUE

    gngt$cex.label = 0.01 ## no labels for full gencode
    gngt$xaxis.chronly = TRUE
    gngt$name = "genes"
    gngt$ywid = 0.1
    gngt$height = 2
    gngt$yaxis.cex = 0.8

    if (!is.null(cgc.gt)) {
        cgc.gt$cex.label = 0.5
        cgc.gt$xaxis.chronly = TRUE
        cgc.gt$name = "CGC"
        cgc.gt$ywid = 0.1
        cgc.gt$height = 5
        cgc.gt$stack.gap = 0.5
        gngt$yaxis.cex = 0.8

        ## concatenate final gTracks
        gt = c(gngt, cgc.gt, cvgt, this.complex.gt)
    } else {
        gt = c(gngt, cvgt, this.complex.gt)
    }

    if (!is.null(agt)) {
        agt$ylab = "CN"
        agt$yaxis.pretty = 3
        agt$xaxis.chronly = TRUE

        gt = c(agt, gt)
    }

    ## save plots
    pts = lapply(1:nrow(complex.ev),
                 function(ix) {
                     ## prepare window
                     win = parse.grl(complex.ev$footprint[[ix]]) %>% unlist
                     if (pad > 0 & pad <= 1) {
                         adjust = pmax(1e5, pad * width(win))
                         win = GenomicRanges::trim(win + adjust)
                     } else {
                         win = GenomicRanges::trim(win + pad)
                     }
                     ppng(plot(gt, win, legend.params = list(plot = FALSE)),
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
#' @param color (character) color of sample burden line default red
#' @param lwd (numeric) size (width) of sample burden line default 2
#' @param outdir (character) ridge plot output directory 
ridge.plot = function(complex.fname = NULL,
                      background.fname = "/data/sv.burden.txt",
                      ev.types = c("qrp", "tic", "qpdup", "qrdel",
                                   "bfb", "dm", "chromoplexy", "chromothripsis",
                                   "tyfonas", "rigma", "pyrgo"),
                      height = 1000,
                      width = 1000,
                      color = "red",
                      lwd = 2,
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
        this.burden.dt = data.table(type = ev.types, burden = ifelse(ev.types %in% names(ev.counts), ev.counts[ev.types], 0))
    }

    ## annotate each with a percentile
    percentile.fn = lapply(setNames(this.burden.dt$type, this.burden.dt$type),
                           function (t) {
                               return(ecdf(sv.burden.dt[type == t, burden]))
                           })

    nt = sapply(1:nrow(this.burden.dt), function(x) {
        percentile.fn[[this.burden.dt$type[x]]](this.burden.dt$burden[x])
    })
    
    this.burden.dt[, ntile := nt]
    this.burden.dt[, ntile.label := paste("Percentile:", format(ntile * 100, digits = 4), "%")]
    this.burden.dt[, count.label := paste("Count:", burden)]
    this.burden.dt[, final.label := ifelse(burden > 0,
                                           paste(count.label, ntile.label, sep = "\n"),
                                           count.label)]

    ## synchronizeg factor levels
    sv.burden.dt[, type := factor(type, levels = ev.types)]
    this.burden.dt[, type := factor(type, levels = ev.types)]

    ## prepare plot
    pt = ggplot(sv.burden.dt[type %in% ev.types,], aes(x = burden, y = type, fill = type)) +
        geom_density_ridges(bandwidth = 0.1, alpha = 0.5, scale = 0.9,
                            rel_min_height = 0.01, color = "white",
                            jittered_points = TRUE,
                            position = position_points_jitter(width = 0.01, height = 0),
                            point_shape = '|', point_size = 3, point_alpha = 0.1, point_colour = "black") +
        geom_segment(data = this.burden.dt[ type %in% ev.types],
                     aes(x = burden, xend = burden, y = as.numeric(type), yend = as.numeric(type) + 0.9),
                     color = color, size = lwd) +
        geom_label(data = this.burden.dt[type %in% ev.types],
                   aes(x = burden, y = as.numeric(type) + 0.5, label = final.label),
                   nudge_x = 0.2,
                   hjust = "left",
                   color = "black",
                   label.size = 0,
                   alpha = 0.5) + 
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
