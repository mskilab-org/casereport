#' zchoo Monday, Apr 26, 2021 02:24:31 PM
## This is a script to generate .png files for SV gallery

#' @name create_genes_gtrack
#' @title create_genes_gtrack
#'
#' @param gencode.fname (character) path to GENCODE GFF
#' @param genes (character vector) genes to mark using "label" GRangesList metadata column (all other genes will have label == NA)
#' @param coding.only (logical) only include entries with gene_type == 'protein_coding'
#' @param level.thresh (numeric) only include entries with level <= level.thresh
#' @param cached.dir (character) path to directory in which to save the output gTrack
#' @param verbose (logical)
#'
#' @return file name of gencode gTrack
#' @export
create_genes_gtrack = function(genes = c(),
                               gencode.fname = NULL,
                               coding.only = TRUE,
                               level.thresh = 2,
                               cached.dir = Sys.getenv('GENCODE_DIR')) {
    if (is.null(gencode.fname) || (!file.exists(gencode.fname) && !RCurl::url.exists(gencode.fname))) {
        stop('gencode does not exist at: "', gencode.fname,' "')
    }

    message('Generating genes gTrack')
    gencode = rtracklayer::import(gencode.fname)
    # adding gene_label column (see https://github.com/mskilab/gTrack/pull/18 for details)
    gene_names = gencode$gene_name
    gene_names[which(!(gene_names %in% genes))] = NA
    gencode$gene_label = gene_names
    if (isTRUE(coding.only)){
        genode = gencode %Q% (gene_type == 'protein_coding')
        if ('level' %in% names(mcols(gencode))){
            gencode = gencode %Q% (level <= level.thresh)
        }
    }
    genes.gt = gTrack::track.gencode(gencode = gencode,
                                     cached.dir = cached.dir)
    return(genes.gt)
}

#' @name gallery.wrapper
#' @title gallery.wrapper
#'
#' @description
#'
#' wrapper function to prep input for SlickR
#' @param complex.fname (character)
#' @param cvgt.fname (character) coverage gTrack file path
#' @param gngt.fname (character) gencode gTrack file path
#' @param cgcgt.fname (character) CGC gTrack file path
#' @param agt.fname (character) allele gTrack file name
#' @param background.fname (character) text file with event burdens (default sv.burden.txt)
#' @param server (character) path to gGnome.js server url
#' @param pair (character) pair id
#' @param ev.types (character) complex event types
#' @param pad (numeric) window padding in bp default 5e5
#' @param height (numeric) plot height default 1000
#' @param width (numeric) plot width default 1000
#' @param outdir (character) where to save plots
#' @return data.table with columns from complex and additionally plot.fname and plot.link
gallery.wrapper = function(complex.fname = NULL,
                           background.fname = "/data/sv.burden.txt",
                           cvgt.fname = "./coverage.gtrack.rds",
                           gngt.fname = "./data/gt.ge.hg19.rds",
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

    # generate ridge plot
    ridgeplot.fname = ridge.plot(complex.fname = complex.fname,
                                 background.fname = background.fname,
                                 ev.types = ev.types,
                                 height = height,
                                 width = width,
                                 outdir = outdir)

    svplot.dt = sv.plot(complex.fname = complex.fname,
                        cvgt.fname = cvgt.fname,
                        gngt.fname = gngt.fname,
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

#' @name grab.window
#' @title grab.window
#'
#' @description
#'
#' helper function for adjusting the window around range of interest
#'
#' if range falls within a complex event, then return the event footprint
#' if gene is amplified but not in a complex event but is in a complex amplicon:
#' - return amplicon
#' if gene is not in a complex event OR in a complex amplicon:
#' - return the range of the node(s) containing the gene
#'
#' @param gr (GRanges)
#' @param complex.fname (character)
#' @param ploidy (numeric)
#' @param amp.thresh (numeric)
#' @param ev.types (character) list of acceptable event types (try to exclude simple events)
#' @param return.type (character) one of "data.table" or "GRanges"
#'
#' @return data.table with node.fp, amp.fp, and ev.fp for node, event, and amplicon footprints
grab.window = function(gr, complex.fname,
                       amp.thresh = 4,
                       ploidy = 2,
                       ev.types = c(),
                       return.type = "data.table")
{

    if (is.null(complex.fname) || !file.exists(complex.fname)) {
        stop("complex.fname is not valid")
    }
    if (!(return.type %in% c("data.table", "GRanges"))) {
        stop("invalid return type")
    }

    ## copy to prevent mutation
    gr = copy(gr)

    ## check if zero-length input given
    if (length(gr)== 0) {
        gr$ev.fp = c()
        gr$amp.fp = c()
        gr$node.fp = c()
        gr$win = c()
        if (return.type == "data.table") {
            return (gr2dt(gr))
        }
        return(gr)
    }

    ## read complex

    gg = readRDS(complex.fname)

    ## grab events as GRanges
    if (!is.null(gg$meta$events) && gg$meta$events[, .N]){
        evs = gg$meta$events[type %in% ev.types]
        if (nrow(evs) > 0) {
            ev.grl = parse.grl(evs$footprint)
            values(ev.grl) = evs
            ev.gr = stack(ev.grl)

            ## warning: for genes overlapping multiple events, will pull ALL footprints (potentially huge :
            ev.ov = gr.findoverlaps(gr, ev.gr, scol = c("footprint"), return.type = "data.table")
            if (nrow(ev.ov) > 0) {
                tmp = ev.ov[, .(footprint = paste(unique(footprint), collapse = ",")),
                            by = query.id]
                gr$ev.fp = tmp[match(1:length(gr), tmp$query.id), footprint]
            } else {
                gr$ev.fp = NA_character_ ## keep consistent type
            }
        } else {
            gr$ev.fp = NA_character_
        }
    } else {
        gr$ev.fp = NA_character_
    }


    ## grab amplicons as GRanges
    keep = (gg$nodes$dt$cn/ploidy) > amp.thresh
    if (any(keep, na.rm = TRUE)){
        gg$clusters(keep)
        ## grab nodes with non-NA cluster
        amp.gr = (gg$nodes$gr %Q% (!is.na(cluster))) %>% gr.stripstrand

        if (length(amp.gr) > 0) {

            ## get footprint of the whole cluster of easier querying later on
            amp.dt = gr2dt(amp.gr)
            amp.dt$grstr = gr.string(amp.gr)
            amp.dt[, footprint := paste(unique(grstr), collapse = ","), by = "cluster"]
            amp.gr$footprint = amp.dt$footprint

            ## get overlaps with genes
            amp.ov = gr.findoverlaps(gr, amp.gr, scol = c("cluster", "footprint"), return.type = "data.table")
            if (nrow(amp.ov) > 0) {
                tmp.amp = amp.ov[, .(footprint = paste(unique(footprint), collapse = ",")), by = query.id]
                gr$amp.fp = tmp.amp[match(1:length(gr), tmp.amp$query.id), footprint]
            } else {
                gr$amp.fp = NA_character_
            }
        } else {
            gr$amp.fp = NA_character_
        }
    } else {
        gr$amp.fp = NA_character_
    }

  
  ## grab node footprint
    node.gr = gg$nodes$gr[, c()]
    node.gr$footprint = gr.string(node.gr)
    node.ov = gr.findoverlaps(gr, node.gr, scol = c("footprint"), return.type = "data.table")
    if (nrow(node.ov) > 0) {
        tmp.node = node.ov[, .(footprint = paste(unique(footprint)), collapse = ","), by = query.id]
        gr$node.fp = tmp.node[match(1:length(gr), tmp.node$query.id), footprint]
    } else {
        gr$node.fp = NA_character_
    }

    ## get final window
    gr$win = ifelse(!is.na(gr$ev.fp),
                    gr$ev.fp,
             ifelse(!is.na(gr$amp.fp),
                    gr$amp.fp,
                    gr$node.fp))

    ## get event id, amplicon id, and node id for each range
    if (return.type == "data.table") {
        return(gr2dt(gr))
    }
    return(gr)
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
#' @param agt.fname (character) allele gTrack file name
#' @param server (character) path to gGnome.js server url
#' @param pair (character) pair id
#' @param ev.types (character) complex event types
#' @param amp.thresh (default 2) amplicon CN threshold for plotting ## possibly compute and save separately?
#' @param ploidy (default 2)
#' @param pad (numeric) window padding in bp default 5e5
#' @param height (numeric) plot height default 500
#' @param width (numeric) plot width default 500
#' @param outdir (character) where to save plots?
#' @param verbose (logical) default TRUE
#'
#' @return data.table with columns id, plot.fname, and plot.link
cn.plot = function(drivers.fname = NULL,
                   complex.fname = NULL,
                   cvgt.fname = "./coverage.gtrack.rds",
                   gngt.fname = "./data/gt.ge.hg19.rds",
                   agt.fname = NULL,
                   server = "",
                   pair = "",
                   gg.rds="",
                   ev.types = c("qrp", "tic", "qpdup", "qrdel",
                                "bfb", "dm", "chromoplexy", "chromothripsis",
                                "tyfonas", "rigma", "pyrgo", "cpxdm"),
                   amp.thresh = 2,
                   ploidy = 2,
                   pad = 0.5,
                   height = 1000,
                   width = 1000,
                   outdir = "./",
                   overwrite  = TRUE,
                   verbose = TRUE) {
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
        drivers.dt[, plot.link := paste0('//', server, "index.html?file=", pair, ".json&location=", ev.js.range, "&view=")]

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
        cvgt$y0 = 0

        this.complex.gt$ylab = "CN"
        this.complex.gt$name = "JaBbA"
        this.complex.gt$yaxis.pretty = 3
        this.complex.gt$xaxis.chronly = TRUE
        this.complex.gt$y0 = 0

        gngt$xaxis.chronly = TRUE
        gngt$labels.suppress.grl = TRUE
        gngt$name = "genes"
        gngt$ywid = 0.1
        gngt$height = 4
        gngt$yaxis.cex = 0.8
        gngt$grl.labelfield = ''

        ## form gencode gTrack for just over/underexpressed genes
        drivers.gt.data = gngt@data[[1]][names(gngt@data[[1]]) %in% drivers.dt[, gene_name]]
        drivers.gt = gngt
        drivers.gt@data[[1]] = drivers.gt.data
        drivers.gt$labels.suppress = TRUE
        drivers.gt$name = "drivers"
        drivers.gt$labels.suppress.gr = TRUE
        drivers.gt$labels.suppress.grl = FALSE
        drivers.gt$height = 5
        drivers.gt$xaxis.chronly = TRUE
        drivers.gt$ywid = 0.1
        drivers.gt$yaxis.cex = 0.8

        if (!is.null(agt)) {
            agt$ylab = "CN"
            agt$yaxis.pretty = 3
            agt$xaxis.chronly = TRUE
            agt$y0 = 0
        }

        ## grab windows
        if (verbose) {
            message("Grabbing windows for plotting")
        }
        win.gr.other = grab.window(gr = drivers.gr, complex.fname = complex.fname,
                             return.type = "GRanges", amp.thresh = amp.thresh,
                             ploidy = ploidy, ev.types = ev.types)
        
        gg=gG(jabba = gg.rds)
        gg2=gg$copy$subgraph(trim(gg$nodes[gg$nodes$dt$cn<gg$meta$ploidy]$gr+pad))
        gg2$clusters()
        win.gr=GRangesList()
        for(gene in drivers.dt$gene_name){
            sqname=drivers.dt[gene_name==gene]$seqnames
            sstart=drivers.dt[gene_name==gene]$start
            send=drivers.dt[gene_name==gene]$end
            cl = (gg2$nodes %&% GRanges(seqnames = sqname, ranges = IRanges(start = sstart, end = send)))$dt$cluster
            thisWin = gg2$nodes[gg2$nodes$dt$cluster == cl]$footprint
            win.gr[[gene]]=thisWin
        }
        drivers.gr$win = win.gr.other$win ## copy over plot window possibly return vector in the future
        ## make one plot per range in drivers gr
        pts = lapply(1:length(drivers.gr),
                     function(ix) {

                         if (!file.exists(drivers.dt$plot.fname[ix]) | overwrite){
                             ## prepare window
                             if(grepl("homdel",drivers.dt$cnv[ix]) | grepl("hetdel",drivers.dt$cnv[ix])){
                                win=win.gr[[ix]]
                            }else{
                                if (is.null(drivers.gr$win) || is.na(drivers.gr$win[[ix]])) {
                                    win = drivers.gr[ix]
                                } else {
                                     win = parse.grl(drivers.gr$win[[ix]]) %>% stack %>% gr.stripstrand
                                }

                                if (pad > 0 & pad <= 1) {
                                 adjust = pmax(5e5, pad * width(win))
                                 message(gr.string(win))
                                 message("adjust size: ", adjust)
                                 win = GenomicRanges::trim(win + adjust)
                                 message(gr.string(win))
                                } else {
                                 win = GenomicRanges::trim(win + pad)
                                }
                            }

                             # assigning greater cex value to the driver gene in order to highlight it 
                             drivers.df = values(drivers.gt@data[[1]])
                             win_width = sum(width(win))
                             cex.val = 0.3 # if the window is very wide then we will only show the highlight driver
                             cex.vals = rep(cex.val, length(values(drivers.gt@data[[1]])$id))
                             cex.val = ifelse(win_width > 100e6, 0.4, 0.3) # if the window is very wide then we will only show the highlight driver
                             if (win_width > 20e9){
                                 # if the window is really wide then we remove the labels from other drivers
                                 # we could do this in a smarter way where only in cases in which labels overlap then we remove them, but it doesn't seem worth the effort
                                 driver.labels = rep('', length(values(drivers.gt@data[[1]])$id))
                                 names(driver.labels) = values(drivers.gt@data[[1]])$id
                                 driver.labels[as.character(drivers.gr$gene_name[ix])] = as.character(drivers.gr$gene_name[ix])
                                 drivers.df$driver.label = driver.labels
                                 values(drivers.gt@data[[1]]) = drivers.df
                                 drivers.gt$grl.labelfield = 'driver.label'
                             } else {
                                 drivers.gt$grl.labelfield = 'id'
                             }
                             names(cex.vals) = values(drivers.gt@data[[1]])$id
                             cex.vals[as.character(drivers.gr$gene_name[ix])] = 0.6
                             values(drivers.gt@data[[1]])$cex.label = cex.vals
                             drivers.gt$grl.cexfield = 'cex.label'
                             
                             # TODO: print only the specific driver gene if the window larger than win.thresh
                             gt = c(gngt, drivers.gt, cvgt, this.complex.gt)
                             if (!is.null(agt)){
                                gt = c(agt, gt)
                             }
                             gt$stack.gap = win_width / 10
                             ppng(plot(gt, win, legend.params = list(plot = FALSE)),
                                  title = drivers.gr$gene_name[ix], ## title is the gene name
                                  filename = drivers.dt$plot.fname[ix],
                                  height = height,
                                  width = width)
                         }
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
#' @param genes.gt.fname (character) CGC gTrack file path
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
                   genes.gt.fname = NULL,
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

    if (is.null(this.complex$meta$events) || !nrow(this.complex$meta$events)){
        return(list())
    }
    ## extract just complex events
    complex.ev = this.complex$meta$events[type %in% ev.types,]


    if (nrow(complex.ev)>0) {

        ## grab plot titles and file names

        complex.ev[, plot.fname := file.path(outdir, paste("event", ev.id, "png", sep = "."))]
        complex.ev[, plot.title := paste(type, "|", "event id:", ev.id)]

        ## make gGnome.js url query parameters
        complex.ev[, ev.js.range := {
            gr = gr.reduce(gr.stripstrand(parse.gr(footprint)) + 1e4)
            paste(gr.string(gr), collapse = "%20|%20")
        }, by = ev.id]
        complex.ev[, plot.link := paste0('//', server, "index.html?file=", pair, ".json&location=", ev.js.range, "&view=")]

        ## read CGC gTrack if provided
        if (!is.null(genes.gt.fname)) {
            if (file.exists(genes.gt.fname)) {
                genes.gt = readRDS(genes.gt.fname)
            } else {
                genes.gt = NULL
            }
        } else {
            genes.gt = NULL
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

        if (!is.null(genes.gt)) {
            genes.gt$cex.label = 0.5
            genes.gt$labels.suppress = FALSE
            genes.gt$xaxis.chronly = TRUE
            genes.gt$name = "Genes"
            genes.gt$ywid = 0.1
            genes.gt$height = 5
            genes.gt$stack.gap = 1e6
            gngt$yaxis.cex = 0.8

            ## concatenate final gTracks
            gt = c(genes.gt, cvgt, this.complex.gt)
        } else {
            gt = c(cvgt, this.complex.gt)
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
                             adjust = pmax(1e6, pad * width(win))
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
        return(complex.ev[, .(type, plot.fname, plot.link)])
    } else {
        return(data.table(plot.fname = character(0), plot.link = character(0)))
    }


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
