#' zchoo Tuesday, Apr 27, 2021 10:49:13 AM
#' this is to generate data tables and plots for fusions

#' @name fusion.wrapper
#' @title fusion.wrapper
#'
#' @description
#'
#' Wrapper for fusions
#'
#' @param fusions.fname (character)
#' @param complex.fname (character)
#' @param cvgt.fname (character)
#' @param gngt.fname (character)
#' @param cgc.fname (character)
#' @param ev.types (character)
#' @param pad (numeric)
#' @param height (numeric)
#' @param width (numeric)
#' @param outdir (character)
#'
#' @return data.table
fusion.wrapper = function(fusions.fname = NULL,
                          complex.fname = NULL,
                          cvgt.fname = NULL,
                          gngt.fname = NULL,
                          cgc.fname = "/data/cgc.tsv",
                          ev.types = c("qrp", "qpdup", "qrdel",
                                       "tic", "bfb", "dm", "chromoplexy",
                                       "chromothripsis", "tyfonas", "rigma", "pyrgo"),
                          height = 1000,
                          width = 1000,
                          pad = 1e5,
                          outdir = "./") {

    filtered.fusions = fusion.table(fusions.fname = fusions.fname,
                                    complex.fname = complex.fname,
                                    cgc.fname = cgc.fname,
                                    ev.types = ev.types,
                                    outdir = outdir)

    filtered.fusions = fusion.plot(fs = filtered.fusions,
                                   complex.fname = complex.fname,
                                   cvgt.fname = cvgt.fname,
                                   gngt.fname = gngt.fname,
                                   pad = pad,
                                   height = height,
                                   width = width,
                                   outdir = outdir)

    return (filtered.fusions$dt)
}

#' @name fusion.table
#' @title fusion.table
#'
#' @description
#'
#' Prepare table of fusion events for displaying
#'
#' @param fusions.fname (character) file name containing fusion gWalks
#' @param complex.fname (character) file name to events
#' @param cgc.fname (character) cancer gene census file name
#' @param ev.types (character) event types
#' @param outdir (character) output directory
#'
#' @return gWalk (filtered) with metadata columns:
#' - walk.id
#' - name (e.g. involved genes)
#' - num.aa (total number of aminos)
#' - gene.pc (protein coordinates)
#' - driver (logical, is a gene a driver)
#' - ev.id (complex event IDs overlapping with walk)
#' - ev.type (complex event type overlapping with walk)
fusion.table = function(fusions.fname = NULL,
                        complex.fname = NULL,
                        cgc.fname = "/data/cgc.tsv",
                        ev.types = c("qrp", "qpdup", "qrdel",
                                     "tic", "bfb", "dm", "chromoplexy",
                                     "chromothripsis", "tyfonas", "rigma", "pyrgo"),
                        outdir = "./") {

    if (!file.exists(fusions.fname)) {
        stop("fusions.fname does not exist")
    }
    if (!file.exists(complex.fname)) {
        stop("complex.fname does not exist")
    }
    if (!file.exists(cgc.fname)) {
        stop("cgc.fname does not exist")
    }

    ## filter to include only in-frame non-silent
    this.fusions = readRDS(fusions.fname) ## gWalk object
    filtered.fusions = this.fusions[in.frame == TRUE & silent == FALSE & numgenes > 1]

    ## compute total number of amino acids and mark
    n.aa = sapply(filtered.fusions$dt$gene.pc, function(pc) { sum(width(parse.grl(pc))) })
    filtered.fusions$set(total.aa = n.aa)

    ## filter so that gene pairs are unique (if multiple choose the one with most AA's)
    filtered.fusions = filtered.fusions[order(n.aa, decreasing = TRUE)]
    filtered.fusions = filtered.fusions[which(!duplicated(filtered.fusions$dt$name))]

    ## Cancer Gene Census genes
    cgc.gene.symbols = fread(cgc.fname)[["Gene Symbol"]]

    ## annotate genes if they are in cgc
    cgc.gene = sapply(filtered.fusions$dt$name, function (gns) {any(unlist(strsplit(gns, ",")) %in% cgc.gene.symbols)})

    filtered.fusions$set(driver = cgc.gene)

    ## grab GRanges for each walk
    fs.grl = filtered.fusions$grl
    values(fs.grl) = filtered.fusions$dt[, .(walk.id)]
    fs.gr = stack(fs.grl)
    ## fs.dt = as.data.table(filtered.fusions$grl) ## adds column group
    ## fs.dt[, ":="(walk.id = filtered.fusions$dt$walk.id[group])]

    ## overlap with complex events
    this.ev = readRDS(complex.fname)$meta$events[type %in% ev.types,]
    ev.grl = parse.grl(this.ev$footprint)
    values(ev.grl) = this.ev
    ev.gr = stack(ev.grl)
    ## ev.footprints = as.data.table(parse.grl(this.ev$footprint))
    ## ev.footprints[, ":="(ev.id = this.ev$ev.id[group], type = this.ev$type[group])]

    ov = gr.findoverlaps(fs.gr, ev.gr,
                         qcol = c("walk.id"),
                         scol = c("ev.id", "type"),
                         return.type = "data.table")
    
    ov = ov[, .(ev.id = paste(unique(ev.id), sep = ","), type = paste(unique(type), sep = ",")), by = walk.id]

    dt = merge(filtered.fusions$dt, ov, by = "walk.id", all.x = TRUE)

    ## add ev.id and type to  metadata 
    filtered.fusions$set(ev.id = dt$ev.id)
    filtered.fusions$set(ev.type = dt$type)
    
    return(filtered.fusions)
}

#' @name fusion.plot
#' @title fusion.plot
#'
#' @description
#'
#' Create .png files for each fusion
#'
#' @param fs (gWalk) gWalk object containing filtered fusions
#' @param complex.fname (character)
#' @param cvgt.fname (character) coverage gTrack file name
#' @param gngt.fname (character) gencode gTrack file name
#' @param pad (numeric) gWalk pad for plotting default 1e5
#' @param height (numeric) plot height default 1e3
#' @param width (numeric) plot width default 1e3
#' @param outdir (character) output directory
#'
#' @return gWalk with additional columns plot.fname (input to slickR)
fusion.plot = function(fs = NULL,
                       complex.fname = NULL,
                       cvgt.fname = NULL,
                       gngt.fname = "/data/gt.ge.hg19.rds",
                       pad = 1e5,
                       height = 1e3,
                       width = 1e3,
                       outdir = "./") {

    if (!file.exists(cvgt.fname)) {
        stop("cvgt.fname does not exist")
    }
    if (!file.exists(gngt.fname)) {
        stop("gngt.fname does not exist")
    }
    if (!file.exists(complex.fname)) {
        stop("complex.fname does not exist")
    }

    fs = fs$copy

    ## read gTracks
    cvgt = readRDS(cvgt.fname)
    gngt = readRDS(gngt.fname)
    this.complex.gt = readRDS(complex.fname)$gt
    
    plot.fnames = sapply(seq_along(fs),
                         function (ix) {
                             fn = file.path(outdir, "fusions", paste0("walk", fs$dt$walk.id[ix], ".png"))
                             fs.gt = fs[ix]$gt
                             fs.gt$name = "gWalk"
                             ppng(plot(c(gngt, cvgt, this.complex.gt, fs.gt),
                                       fs[ix]$footprint + pad,
                                       legend.params = list(plot = FALSE)),
                                  title = paste(fs$dt$name[ix], "|", "walk", fs$dt$walk.id[ix]),
                                  filename = fn,
                                  height = height,
                                  width = width)
                             return(fn)
                         })
    fs$set(plot.fname = plot.fnames)
    return(fs)
}
