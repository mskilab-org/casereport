#' zchoo Tuesday, Apr 27, 2021 10:49:13 AM
#' this is to generate data tables and plots for fusions

#' @name star2grl
#' @title star2grl
#'
#' Quick utility function to convert star-fusion breakpoints to GRangesList
#'
#' @param fname (character) name of star-fusion output file
#' @param mc.cores (numeric) default 16
star2grl = function(fname, mc.cores = 16) {
    dt = fread(fname)
    grl = mclapply(1:nrow(dt),
                   function(ix) {
                       left.bp.str = strsplit(dt[ix, LeftBreakpoint], ":")[[1]]
                       right.bp.str = strsplit(dt[ix, RightBreakpoint], ":")[[1]]
                       gr = GRanges(seqnames = c(left.bp.str[1], right.bp.str[1]),
                                    ranges = IRanges(start = as.numeric(c(left.bp.str[2],
                                                                          right.bp.str[2])),
                                                     width = 1),
                                    ## reverse strands to match our strand designation for junctions
                                    strand = c(ifelse(left.bp.str[3] == "+", "-", "+"),
                                               ifelse(right.bp.str[3] == "+", "-", "+"))
                                    )
                       return (gr)
                   }, mc.cores = mc.cores) %>% GRangesList
    values(grl) = dt
    return(grl)
}

#' @name wgs.circos
#' @title wgs.circos
#'
#' Quick utility function for circos plot with read depth, junctions, and segments
#' (copied from skitools)
#' 
#' @param junctions Junction object with optional metadata field  $col to specify color
#' @param cov GRanges of scatter points with optional fields $col
#' @param segs GRanges of segments with optional fields $col and $border
#' @param win GRanges window to limit plot to
#' @param cytoband GRanges of cytoband
#' @param y.field field in cov that specifies the y axis to draw
#' @param cex.points cex for cov points
#' @param max.ranges max ranges for cov points (1e4)
#' @param ylim ylim on cov (default automatically computed)
#' @param cytoband.path path to UCSC style cytoband path
#' @param y.quantile quantile normalization
#' @param chr.sum whether to chr.sub everything 
#' @author Marcin Imielinski
#' @export
wgs.circos = function(junctions = jJ(),
                      cov = NULL,
                      segs = NULL,
                      win = NULL,
                      field = 'ratio',
                      cytoband = NULL,
                      y.field = field,
                      ylim = NA,
                      cytoband.path = '~/DB/UCSC/hg19.cytoband.txt',
                      cex.points = 1,
                      ideogram.outer = TRUE,
                      scatter = TRUE,
                      bar = FALSE,
                      line = FALSE,
                      gap.after = 1,
                      labels.cex = 1,
                      y.quantile = 0.9999,
                      chr.sub = TRUE,
                      max.ranges = 1e4,
                      axis.frac = 0.02,
                      palette = 'BrBg', ...)
{

    if (!file.exists(cytoband.path))
        stop('cytoband not file, must be UCSC style tsv')

    if (is.null(cytoband))
        cytoband = circlize::read.cytoband(cytoband.path)$df

    cytoband = as.data.table(cytoband)
    setnames(cytoband, c('seqnames', 'start', 'end', 'band', 'stain'))

    if (chr.sub)
        cytoband[, seqnames := gsub('chr', '', seqnames)]
    
    if (!is.null(win))
    {
        if (is.character(win) | is.integer(win) | is.numeric(win) | is.factor(win))
            win = parse.gr(as.character(win))

        if (inherits(win, 'data.frame'))
            win = dt2gr(win)

        cytoband  = as.data.table(dt2gr(cytoband) %*% win)[, .(seqnames, start, end, band, stain)]
    }

    total.width = cytoband[, sum(as.numeric(end-start))]
    if (!is.na(axis.frac) && axis.frac>0)
    {
        axis.width = ceiling(axis.frac*total.width)
        cytoband = rbind(cytoband, data.table(seqnames = 'axis', start = 0, end = axis.width, band = '', stain = ''), fill = TRUE)
    }

    if (chr.sub)
    {
        ix = ((junctions$left %>% gr.sub('chr', ''))  %^% dt2gr(cytoband)) &
            ((junctions$right %>% gr.sub('chr', '')) %^% dt2gr(cytoband))
        junctions = junctions[ix]
    }
    else
    {
        ix = junctions$left %^% dt2gr(cytoband) & junctions$right %^% dt2gr(cytoband)
        junctions = junctions[ix]
    }

    cytoband[, seqnames := as.character(seqnames)]
    args  = list(...)
    ## some important pars
    labels.cex = ifelse(is.null(args$labels.cex), 1, args$labels.cex)
    bands.height = ifelse(is.null(args$bands.height), 0.1, args$bands.height)
    cn.height = ifelse(is.null(args$cn.height), 0.3, args$cn.height)
    link.h.ratio = ifelse(is.null(args$link.h.ratio), 0.75, args$link.h.ratio)

    ## mark with colors by class
    col.dt = data.table(class = c("INV-like", "TRA-like", "DUP-like", "DEL-like"),
                        col = c(alpha("purple", 0.5),
                                alpha("green", 0.5),
                                alpha("red", 0.5),
                                alpha("blue", 0.5)))
    bpdt = junctions$dt
    bpdt[, col := col.dt$col[match(bpdt$class, col.dt$class)]]
    
    bp1 = junctions$left %>% gr2dt
    bp2 = junctions$right%>% gr2dt
    circlize::circos.clear()
    circlize::circos.par(start.degree = 90, gap.after = gap.after*1)
    circlize::circos.genomicInitialize(cytoband, sector.names = unique(cytoband$seqnames), plotType = NULL, 
                                       track.height = bands.height,
                                       labels.cex = labels.cex)

    circlize::circos.genomicTrackPlotRegion(cytoband, stack = TRUE,
                                            panel.fun = function(region, value, ...) {
                                                xlim = circlize::get.cell.meta.data("xlim")
                                                ylim = circlize::get.cell.meta.data("ylim")
                                                chr = circlize::get.cell.meta.data("sector.index") %>% gsub('chr', '', .)
                                                if (circlize::get.cell.meta.data("sector.index") != 'axis')
                                                {
                                                    circlize::circos.text(mean(xlim), 0.9, chr, cex = 1.5, facing = "clockwise", adj = c(0,1),
                                                                          niceFacing = TRUE)
                                                }
                                            }, track.height = 0.1, bg.border = NA)

    ## inner ideogram
    if (ideogram.outer)
    {
        circlize::circos.genomicTrackPlotRegion(cytoband, stack = TRUE,
                                                panel.fun = function(region, value, ...) {
                                                    xlim = circlize::get.cell.meta.data("xlim")
                                                    ylim = circlize::get.cell.meta.data("ylim")
                                                    chr = circlize::get.cell.meta.data("sector.index")
                                                    if (circlize::get.cell.meta.data("sector.index") != 'axis')
                                                    {
                                                        at = pretty(xlim, n = 3)
                                                        circlize::circos.axis(direction = "outside", labels.facing = "outside", major.at = at, minor.ticks = 10, labels = (at/1e6) %>% as.integer, labels.cex = labels.cex*0.3)
                                                        circlize::circos.genomicRect(region, value, col =  circlize::cytoband.col(value[[2]]), border = NA)
                                                        circlize::circos.rect(xlim[1], ylim[1], xlim[2], ylim[2], border = "black")
                                                    }
                                                }, track.height = 0.05, bg.border = NA)
    }
    
    ## coverage scatter plot
    if (!is.null(cov))
    {
        if (inherits(cov, 'data.frame'))
            cov = dt2gr(cov)

        cov = cov[!is.na(values(cov)[[y.field]])]
        cov = cov[!is.infinite(values(cov)[[y.field]])]

        if (is.na(ylim))
            ylim = c(0, quantile(values(cov)[[y.field]], y.quantile, na.rm = TRUE))
        
        cov$y = values(cov)[[y.field]] %>% as.numeric
        cov$y = cov$y %>% pmin(ylim[2]) %>% pmax(ylim[1])

        if (is.null(cov$col))
            cov$col = 'black'

        cov = cov[sample(length(cov), pmin(length(cov), max.ranges))]
        uchr = unique(cytoband$seqnames)
        cov = cov %&% dt2gr(cytoband)
        covdt = gr2dt(cov)[, seqnames := factor(seqnames, uchr)]
        circlize::circos.genomicTrackPlotRegion(covdt[, .(seqnames, start, end, y, as.character(col), ytop = y)],
                                                ylim = ylim,
                                                track.height = cn.height,
                                                bg.border = ifelse(uchr == 'axis', NA, alpha('black', 0.2)),
                                                panel.fun = function(region, value, ...) {
                                                    if (circlize::get.cell.meta.data("sector.index") != 'axis')
                                                    {
                                                        if (circlize::get.cell.meta.data("sector.index") == uchr[1])
                                                            circlize::circos.yaxis(side = 'left')                                    
                                                        if (scatter)
                                                            circlize::circos.genomicPoints(region, value, numeric.column = 1, col = value[[2]], pch = 16, cex = cex.points, ...)
                                                        if (bar)
                                                            circlize::circos.genomicRect(region, value[[1]], ytop.column = 1, border = value[[2]], col = value[[2]], pch = 16, cex = cex.points, ...)
                                                        if (line)
                                                            circlize::circos.genomicLines(region, value[[1]], col = value[[2]], pch = 16, cex = cex.points, ...)
                                                    }
                                                })
    }
    circlize::circos.par(cell.padding = c(0, 0, 0, 0))

    if (!is.null(segs))
    {
        if (inherits(segs, 'data.frame'))
            segs = dt2gr(segs)

        if (chr.sub)
            segs = segs %>% gr.sub('chr', '')

        segs = segs[segs %^% dt2gr(cytoband), ]

        segs = as.data.table(segs)
        if (is.null(segs$col))
            segs$col = 'gray'

        if (is.null(segs$border))
            segs$border = segs$col

        if (chr.sub)
            segs[, seqnames := gsub('chr', '', seqnames)]

        circlize::circos.genomicTrackPlotRegion(segs[, .(seqnames, start, end, col, border)], stack = TRUE,
                                                panel.fun = function(region, value, ...) {
                                                    circlize::circos.genomicRect(region, value, col = value[[1]], border = value[[2]])
                                                    xlim = circlize::get.cell.meta.data("xlim")
                                                    ylim = circlize::get.cell.meta.data("ylim")
                                                    chr = circlize::get.cell.meta.data("sector.index")
                                        #                                    circlize::circos.rect(xlim[1], ylim[1], xlim[2], ylim[2], border = "black")
                                                }, track.height = 0.05, bg.border = NA)
    }

    circlize::circos.par(cell.padding = c(0, 0, 0, 0))


    ## inner ideogram
    if (!ideogram.outer)
    {
        circlize::circos.genomicTrackPlotRegion(cytoband, stack = TRUE,
                                                panel.fun = function(region, value, ...) {
                                                    xlim = circlize::get.cell.meta.data("xlim")
                                                    ylim = circlize::get.cell.meta.data("ylim")
                                                    chr = circlize::get.cell.meta.data("sector.index")
                                                    if (circlize::get.cell.meta.data("sector.index") != 'axis')
                                                    {
                                                        at = pretty(xlim, n = 3)
                                                        circlize::circos.axis(direction = "outside", labels.facing = "outside", major.at = at, minor.ticks = 10, labels = (at/1e6) %>% as.integer, labels.cex = labels.cex*0.3)
                                                        circlize::circos.genomicRect(region, value, col = circlize::cytoband.col(value[[2]]), border = NA)
                                                        circlize::circos.rect(xlim[1], ylim[1], xlim[2], ylim[2], border = "black")
                                                    }
                                                }, track.height = 0.05, bg.border = NA)
    }

    if (nrow(bpdt))
    {

        if (is.null(bpdt$lwd))
            bpdt$lwd = NA_integer_

        bpdt[is.na(lwd), lwd := 2]

        if (is.null(bpdt$col))
            bpdt$col = NA_character_

        bpdt[is.na(col), col := 'red']

        if (is.null(bpdt$lty))
            bpdt$lty = NA_integer_

        bpdt[is.na(lty), lty := 1]

        if (nrow(bpdt))
            bpdt$span  = cut(junctions$span, c(0, 1e6, 3e8, Inf))

        spmap = structure(c(0.05, 0.2, 1), names = levels(bpdt$span))
        ixs = split(1:nrow(bpdt), bpdt$span)
        lapply(names(ixs), function(i)
            circlize::circos.genomicLink(
                bp1[ixs[[i]], .(seqnames, start, end)],
                bp2[ixs[[i]], .(seqnames, start, end)],
                h = spmap[i],
                                        #                  rou = circlize:::get_most_inside_radius()*c(0.1, 0.5, 1)[bpdt$span[ixs[[i]]] %>% as.integer],
                col = bpdt[ixs[[i]], ]$col,
                lwd =  bpdt[ixs[[i]], ]$lwd,
                lty =  bpdt[ixs[[i]], ]$lty,
                h.ratio = link.h.ratio,
                border=NA)
            )
    }

    ## create legend for junctions
    lgd_links = ComplexHeatmap::Legend(
        at = col.dt$class,
        type = "lines",
        legend_gp = gpar(col = col.dt$col, lwd = 5),
        title_position = "topleft",
        title = "Junctions",
        labels_gp = gpar(col = "black", fontsize = 15),
        title_gp = gpar(col = "black", fontsize = 25),
        grid_width = unit(10, "mm"))

    
    draw(lgd_links,
         x = unit(1, "npc") - unit(4, "mm"),
         y = unit(4, "mm"),
         just = c("right", "bottom"))

    
    circlize::circos.clear()
}

#' @name fusion.wrapper
#' @title fusion.wrapper
#'
#' @description
#'
#' Wrapper for fusions
#'
#' @param fusions.fname (character)
#' @param complex.fname (character)
#' @param cvgt (character or gTrack)
#' @param agt.fname (character)
#' @param gngt gTrack object or character describing path to RDS file containing a gTrack object
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
                          cvgt = NULL,
                          gngt = NULL,
                          agt.fname = NULL,
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
                                    ev.types = ev.types)

    filtered.fusions = fusion.plot(fs = filtered.fusions,
                                   complex.fname = complex.fname,
                                   cvgt = cvgt,
                                   gngt = gngt,
                                   agt.fname = agt.fname,
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
                                     "chromothripsis", "tyfonas", "rigma", "pyrgo"))
{
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
    cgc.dt = rbindlist(
        lapply(1:length(filtered.fusions),
               function(ix) {
                   ## xtYao: name seems to be integers in new gGnome??
                   ## let's use "genes" moving forward
                   ## gns = unlist(strsplit(filtered.fusions$dt$name[ix], ","))
                   gns = unlist(strsplit(filtered.fusions$dt$genes[ix], ","))
                   gene.in.cgc = any(gns %in% cgc.gene.symbols)
                   gns.filtered = gns[which(gns %in% cgc.gene.symbols)]
                   cgc.names = paste(gns.filtered, collapse = ", ")
                   return(data.table(driver = gene.in.cgc, driver.name = cgc.names))
               }),
        fill = TRUE)

    filtered.fusions$set(driver = cgc.dt$driver, driver.name = cgc.dt$driver.name)

    ## get chromosomes touched by each
    fs.chroms = sapply(1:length(filtered.fusions),
                       function(ix) {
                           ix.seqnames = seqnames(filtered.fusions[ix]$grl[[1]]) %>% as.character
                           ix.split.seqnames = split(ix.seqnames,
                                                     cumsum(c(1, abs(diff(as.numeric(as.factor(ix.seqnames)))))))
                           out = lapply(1:length(ix.split.seqnames),
                                        function(j) {unique(ix.split.seqnames[[j]])})
                           return(paste(out, collapse = "->"))
                       })

    filtered.fusions$set(chroms = fs.chroms)

    ## grab GRanges for each walk
    fs.grl = filtered.fusions$grl
    values(fs.grl) = filtered.fusions$dt[, .(walk.id)]
    fs.gr = stack(fs.grl)

    ## overlap with complex events
    this.ev = readRDS(complex.fname)$meta$events[type %in% ev.types,]
    ev.grl = parse.grl(this.ev$footprint)
    values(ev.grl) = this.ev
    ev.gr = stack(ev.grl)

    ov = gr.findoverlaps(fs.gr, ev.gr,
                         qcol = c("walk.id"),
                         scol = c("ev.id", "type"),
                         return.type = "data.table")


    if (ov[,.N] > 0){
        ov = ov[, .(ev.id = paste(unique(ev.id), sep = ","), type = paste(unique(type), sep = ",")), by = walk.id]
        dt = merge(filtered.fusions$dt, ov, by = "walk.id", all.x = TRUE)
    } else {
        dt = filtered.fusions$dt
    }

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
#' @param cvgt gTrack object or character describing path to RDS file containing the coverage gTrack object
#' @param gngt gTrack object or character describing path to RDS file containing a gTrack object
#' @param agt.fname (character) allele gTrack file name
#' @param pad (numeric) gWalk pad for plotting default 1e5
#' @param height (numeric) plot height default 1e3
#' @param width (numeric) plot width default 1e3
#' @param outdir (character) output directory
#' @param prefix (character) a charachter prefix to use for file names
#' @param show.walk (logical) plot a track for the walk
#' @param yaxis.pretty vector or scalar positive integer specifying how many ticks to optimally draw on yaxis (formatting)(logical) plot a track for the walk. When provided this will modify the value for agt cvgt and the complex gTrack. The default is 3. If you wish to maintain the unique value from each original gTrack then set yaxis.pretty = NULL.
#' @param zoom_on_junctions (logical) instead of taking the footprint of the walk, take a window around the coordinates of the ALT junctions (window width is determined by the pad argument)
#' @param show.junction.support (logical) compute and plot junction.support (junction.support is done with realign = FALSE for speed). You must also supply a tumor.bam, ref_bwa, and ref_seq in order to run junction.support
#' @param ref_bwa (character or RSeqLib::BWA) a path to an RDS containing the RSeqLib::BWA object for the reference genome
#' @param ref_seq (character or DNAStringSet) DNAStringSet or a path to an RDS containing the Biostrings::DNAStringSet object for the reference genome
#' @param tumor.bam (character) path to the tumor BAM file. This is only used when show.junction.support == TRUE.
#'
#' @return gWalk with additional columns plot.fname (input to slickR)
fusion.plot = function(fs = NULL,
                       complex.fname = NULL,
                       cvgt = NULL,
                       agt.fname = NULL,
                       gngt = "/data/gt.ge.hg19.rds",
                       pad = 1e5,
                       height = NULL,
                       width = NULL,
                       outdir = "./",
                       prefix = '',
                       show.walk = TRUE,
                       yaxis.pretty = 3,
                       zoom_on_junctions = FALSE,
                       show.junction.support = FALSE,
                       ref_bwa = NULL,
                       ref_seq = NULL,
                       tumor.bam = NULL,
                       plot.type = 'png') {

    if (!inherits(cvgt, 'gTrack') & !is.null(cvgt)){
        if (is.character(cvgt)){
            if (!file.exists(cvgt)) {
               stop('cvgt input not valid: ', cvgt, ' does not exist')
            }
       } else {
            stop("cvgt not valid. cvgt must be either a gTrack object or character describing path to RDS file containing a gTrack object, but an objcet of class \"", paste(class(cvgt), collapse = ' '), "\" was provided.")
       }
        cvgt = readRDS(cvgt)
    }

    if (!inherits(gngt, 'gTrack')){
        if (is.character(gngt)){
            if (!file.exists(gngt)) {
               stop('gngt input not valid: ', gngt, ' does not exist')
            }
       } else {
            stop("gngt not valid. gngt must be either a gTrack object or character describing path to RDS file containing a gTrack object, but an objcet of class \"", paste(class(gngt), collapse = ' '), "\" was provided.")
       }
        gngt = readRDS(gngt)
    }

    if (!file.exists(complex.fname)) {
        stop("complex.fname does not exist")
    }

    if (show.junction.support & is.character(tumor.bam) & !is.null(ref_bwa) & !is.null(ref_seq)){
        # if junction.support was requested then let's load the ref_bwa and ref_seq
        if (is.character(ref_bwa)){
            if (!file.exists(ref_bwa)){
                stop('Invalid ref_bwa. ', ref_bwa, ' does not exist.')
            }
            message('Loading reference BWA object.')
            bwa_object = readRDS(ref_bwa)
        } else {
            if (inherits(ref_bwa, 'BWA')){
                bwa_object = ref_bwa
            } else {
                stop('ref_bwa must be a RSeqLib::BWA object or a path to an RDS containing such object, but \"', paste(class(ref_bwa), collapse = ' '), '\" was provided.')
            }
        }
        if (is.character(ref_seq)){
            if (!file.exists(ref_seq)){
                stop('Invalid ref_seq. ', ref_seq, ' does not exist.')
            }
            message('Loading reference DNAStringSet object.')
            ref_seq = readRDS(ref_seq)
        } else {
            if (!inherits(ref_seq, 'DNAStringSet')){
                stop('ref_seq must be a DNAStringSet object or a path to an RDS containing such object, but \"', paste(class(ref_seq), collapse = ' '), '\" was provided.')
            }
        }
        if (!file.exists(tumor.bam)){
            stop('Invalid tumor.bam. ', tumor.bam, ' does not exist.')
        }
    } else {
        if (show.junction.support){
            message('Setting show.junction.support to FALSE since not all required info was provided. See help menu.')
            show.junction.support = FALSE
        }
    }

    gngt$xaxis.chronly = TRUE
    gngt$name = "genes"

    base.gt = gngt # we always plot the genes track

    ## read allele gTrack if provided
    if (!is.null(agt.fname)) {
        if (file.exists(agt.fname)) {
            agt = readRDS(agt.fname)
            agt$ylab = "CN"
            if (!is.null(yaxis.pretty)){
                agt$yaxis.pretty = yaxis.pretty
            }
            agt$xaxis.chronly = TRUE
            base.gt = c(base.gt, agt)
        } else {
            agt = NULL
        }
    } else {
        agt = NULL
    }

    fs = fs$copy

    if (!is.null(cvgt)){
        ## format gTrack
        cvgt$ylab = "CN"
        cvgt$name = "cov"
        if (!is.null(yaxis.pretty)){
            cvgt$yaxis.pretty = yaxis.pretty
        }
        cvgt$xaxis.chronly = TRUE
        # add to the baseline gTrack
        base.gt = c(base.gt, cvgt)
    }
    if (!is.null(complex.fname)){
        this.complex.gt = readRDS(complex.fname)$gt
        ## format gTrack
        this.complex.gt$ylab = "CN"
        this.complex.gt$name = "JaBbA"
        if (!is.null(yaxis.pretty)){
            this.complex.gt$yaxis.pretty = yaxis.pretty
        }
        this.complex.gt$chronly = TRUE
        # add to the baseline gTrack
        base.gt = c(base.gt, this.complex.gt)
    }

    plot.fnames = sapply(seq_along(fs),
                         function (ix) {
                            fn = file.path(dcreate(outdir), "fusions", paste0(prefix, "walk", fs$dt$walk.id[ix], ".", plot.type))
                            if (show.walk == TRUE){
                                fs.gt = fs[ix]$gt

                                ## formatting
                                fs.gt$name = "gWalk"
                                fs.gt$labels.suppress = TRUE
                                fs.gt$labels.suppress.grl = TRUE
                                fs.gt$xaxis.chronly = TRUE

                                gt = c(base.gt, fs.gt)
                            }

                            gt$xaxis.chronly = TRUE

                            if (show.junction.support){
                                message('Computing junction support')
                                j = fs$edges[type == 'ALT']$junctions
                                reads = read.bam(tumor.bam, j$breakpoints + 1e4) %>% grl.unlist
                                reads = reads %Q% which(!is.na(seq))

                                jsupp = junction.support(
                                    reads,
                                    junctions = j,
                                    ref = ref_seq,
                                    pad = 1e3,
                                    realign = FALSE,
                                    bwa = bwa_object)

                               jsupp.gt = gTrack(jsupp %>% split(., .$grl.ix) %>% unname, name = "reads")
                               gt = c(gt, jsupp.gt)
                            }

                            ## format window
                            if (zoom_on_junctions == TRUE){
                                win = fs$edges[type == 'ALT']$junctions$grl %>% grl.unlist
                            } else {
                                win = fs[ix]$footprint
                            }

                            if (pad > 0 & pad <= 1) {
                                adjust = pmax(1e5, pad * width(win))
                                win = GenomicRanges::trim(win + adjust)
                            } else {
                                win = GenomicRanges::trim(win + pad)
                            }

                            if (plot.type == 'png'){
                                plot.fun = skitools::ppng
                                if (is.null(height)) height = 1e3
                                if (is.null(width)) width = 1e3
                            } else {
                                plot.fun = skitools::ppdf
                                if (is.null(height)) height = 10
                                if (is.null(width)) width = 10
                            }

                                plot.fun(plot(gt,
                                     win,
                                     legend.params = list(plot = FALSE)),
                                     title = paste(fs$dt$genes[ix], "|", "walk", fs$dt$walk.id[ix]),
                                     filename = fn,
                                     height = height,
                                     width = width)

                            return(fn)
                         })

    fs$set(plot.fname = plot.fnames)
    return(fs)
}

