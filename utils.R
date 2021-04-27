message("Loading Functions")
##DLRS for coverage quality
dlrs <- function(x) {
    nx <- length(x)
    if (nx<3) {
        stop("Vector length>2 needed for computation")
    }
    tmp <- embed(x,2)
    diffs <- tmp[,2]-tmp[,1]
    dlrs <- IQR(diffs,na.rm=TRUE)/(sqrt(2)*1.34)
    return(dlrs)
}

covcbs = function(x, field = "ratio", name = "sample", max.ranges = 1e4,
                  lwd.border = 0.2, purity = NULL, ploidy = NULL, rebin = NULL,
                  ...){
    x = readRDS(x)
    if (!is.null(purity) & !is.null(ploidy)){
        x$cn = rel2abs(x, purity = purity, ploidy = ploidy, field = field)
        field = "cn"
    }
    if (!is.null(rebin)){
        x = rebin(x, binwidth = rebin, field = field,
                            FUN = median, na.rm = TRUE)
    }
    gTrack(x,
           y.field = field,
           name = name,
           max.ranges = max.ranges,
           circles = TRUE,
           lwd.border = lwd.border,
           ...)
}


fready = function (..., pattern = "\\W+", sub = "_") 
{
    tab = fread(...)
    nms = dedup(gsub(pattern, sub, names(tab), perl = TRUE), 
        suffix = ".") %>% gsub("^[^A-Za-z]", "", ., perl = TRUE)
    setnames(tab, nms)
    return(tab)
}


dedup = function (x, suffix = ".") 
{
    dup = duplicated(x)
    udup = setdiff(unique(x[dup]), NA)
    udup.ix = lapply(udup, function(y) which(x == y))
    udup.suffices = lapply(udup.ix, function(y) c("", paste(suffix, 
        2:length(y), sep = "")))
    out = x
    out[unlist(udup.ix)] = paste(out[unlist(udup.ix)], unlist(udup.suffices), 
        sep = "")
    return(out)
}

rel2abs = function (gr, purity = NA, ploidy = NA, gamma = NA, beta = NA, 
    field = "ratio", field.ncn = "ncn") 
{
    mu = values(gr)[, field]
    mu[is.infinite(mu)] = NA
    w = as.numeric(width(gr))
    w[is.na(mu)] = NA
    sw = sum(w, na.rm = T)
    mutl = sum(mu * w, na.rm = T)
    ncn = rep(2, length(mu))
    if (!is.null(field.ncn)) 
        if (field.ncn %in% names(values(gr))) 
            ncn = values(gr)[, field.ncn]
    ploidy_normal = sum(w * ncn, na.rm = T)/sw
    if (is.na(gamma)) 
        gamma = 2 * (1 - purity)/purity
    if (is.na(beta)) 
        beta = ((1 - purity) * ploidy_normal + purity * ploidy) * 
            sw/(purity * mutl)
    return(beta * mu - ncn * gamma/2)
}


rebin = function (cov, binwidth, field = names(values(cov))[1], FUN = median, na.rm = TRUE) 
{
    tmp = as.data.table(cov[, c()])
    tmp$value = values(cov)[[field]]
    outdt = tmp[, FUN(value, na.rm = na.rm), by = .(seqnames, 
        start = floor(start/binwidth) * binwidth + 1)]
    outdt[, `:=`(end, start + binwidth - 1)]
    out = dt2gr(outdt)
    names(values(out)) = field
    return(out)
}

ppng = function (expr, filename = "plot.png", height = 1000, width = 1000, 
          dim = NULL, cex = 1, title = NULL, cex.pointsize = min(cex), 
          cex.title = 1, ...) 
{
    if (length(cex) == 1) 
        cex = rep(cex, 2)
    height = cex[1] * height
    width = cex[2] * width
    DEFAULT.OUTDIR = Sys.getenv("PPNG.DIR")
    if (nchar(DEFAULT.OUTDIR) == 0) 
        DEFAULT.OUTDIR = normalizePath("~/public_html/")
    if (!grepl("^[~/]", filename)) 
        filename = paste(DEFAULT.OUTDIR, filename, sep = "/")
    if (!file.exists(dirname(filename))) 
        system(paste("mkdir -p", dirname(filename)))
    cat("rendering to", filename, "\n")
    png(filename, height = height, width = width, pointsize = 24 * 
                                                      cex.pointsize, ...)
    if (!is.null(dim)) {
        if (length(dim) == 1) 
            dim = rep(dim, 2)
        dim = dim[1:2]
        layout(matrix(1:prod(dim), nrow = dim[1], ncol = dim[2], 
                      byrow = TRUE))
    }
    eval(expr)
    if (!is.null(title)) 
        title(title, cex.main = cex.title * max(cex))
    dev.off()
}



grok_vcf = function (x, snpeff.ontology = NULL, label = NA, keep.modifier = TRUE, long = FALSE, 
    oneliner = FALSE, verbose = FALSE) 
{
    fn = c("allele", "annotation", "impact", "gene", "gene_id", 
        "feature_type", "feature_id", "transcript_type", "rank", 
        "variant.c", "variant.p", "cdna_pos", "cds_pos", "protein_pos", 
        "distance")
    if (is.character(x)) {
        out = suppressWarnings(skidb::read_vcf(x))
        if (is.na(label)) 
            label = x
    }
    else out = x
    if (is.na(label)) 
        label = ""
    if (verbose) 
        message("Grokking vcf ", label)
    if (!long) {
        vcf = out
        if (length(vcf) > 0) {
            if (!is.null(vcf$ANN)) {
                vcf$eff = unstrsplit(vcf$ANN)
                vcf$modifier = !grepl("(HIGH)|(LOW)|(MODERATE)", 
                  vcf$eff)
                if (!keep.modifier) 
                  vcf = vcf[!vcf$modifier]
            }
            vcf$ref = as.character(vcf$REF)
            vcf$alt = as.character(unstrsplit(vcf$ALT))
            vcf = vcf[, sapply(values(vcf), class) %in% c("factor", 
                "numeric", "integer", "logical", "character")]
            vcf$var.id = 1:length(vcf)
            vcf$type = ifelse(nchar(vcf$ref) == nchar(vcf$alt), 
                "SNV", ifelse(nchar(vcf$ref) < nchar(vcf$alt), 
                  "INS", "DEL"))
            vcf$label = label
        }
        return(vcf)
    }
    else if (length(out) > 0) {
        out$REF = as.character(out$REF)
        out$ALT = as.character(unstrsplit(out$ALT))
        out$vartype = ifelse(nchar(out$REF) == nchar(out$ALT), 
            "SNV", ifelse(nchar(out$REF) < nchar(out$ALT), "INS", 
                "DEL"))
        if (is.null(out$ANN)) 
            stop("no $ANN column, check to see if annotated VCF is formatted in the SnpEff style")
        else out$eff = unstrsplit(out$ANN)
        out$modifier = !grepl("(HIGH)|(LOW)|(MODERATE)", out$eff)
        if (!keep.modifier) 
            out = out[!out$modifier]
        if (inherits(out$ANN, "character")) 
            annlist = strsplit(out$ANN, ",")
        else annlist = out$ANN %>% as.list
        tmp = lapply(annlist, function(y) do.call(rbind, lapply(strsplit(y, 
            "\\|"), "[", 1:15)))
        tmpix = rep(1:length(out), elementNROWS(tmp))
        meta = as.data.frame(do.call(rbind, tmp))
        colnames(meta) = fn
        meta$varid = tmpix
        meta$file = label
        out2 = out[tmpix]
        rownames(meta) = NULL
        values(out2) = cbind(values(out2), meta)
        names(out2) = NULL
        out2$ANN = NULL
        precedence = c("trunc", "cnadel", "cnadup", "complexsv", 
            "splice", "inframe_indel", "fusion", "missense", 
            "promoter", "regulatory", "noncoding", "inv", "synonymous", 
            "")
        ## eff = readRDS(system.file("extdata", "snpeff_ontology.rds", 
        ##     package = "skitools"))[, `:=`(short, factor(short, 
        ##                                                 precedence))][!is.na(short), ]
        if (!is.null(snpeff.ontology)){
            eff = readRDS(snpeff.ontology)
        }
        .short = function(vcf) {
            tmp = strsplit(as.character(vcf$annotation), "\\&")
            dtl = data.table(eff = unlist(tmp), id = rep(1:length(tmp), 
                lengths(tmp))) %>% merge(eff, by = "eff", allow.cartesian = TRUE) %>% 
                unique(by = "id")
            setkey(dtl, id)
            vcf$short = dtl[.(1:length(vcf)), short]
            return(vcf)
        }
        out2 = .short(out2)
        if (oneliner) 
            out2$oneliner = paste(ifelse(!is.na(out2$gene), as.character(out2$gene), 
                as.character(out2$annotation)), ifelse(nchar(as.character(out2$variant.p)) > 
                0, as.character(out2$variant.p), as.character(out2$variant.c)))
    }
    return(out2)
}


## convert a GRanges to the string for gGnome.js browser
js.range = function(gr){
    if (is.character(gr)){
        gr = parse.gr(gr)
    }
    gr = gr.stripstrand(gr.reduce(gr[, c()]))
    if (any(width(gr)==1)){
        gr = gr + 1
    }
    paste(gr.string(gr), collapse = "%20|%20")
}


#' @name circos
#' @title circos
#'
#' Quick utility function for circos plot with read depth, junctions, and segments
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
circos = function(junctions = jJ(), cov = NULL, ncov = NULL, segs = NULL, win = NULL, field = 'ratio', cytoband = NULL, y.field = field, ylim = NA, cytoband.path = '~/DB/UCSC/hg19.cytoband.txt', cex.points = 1, ideogram.outer = TRUE, scatter = TRUE, bar = FALSE, line = FALSE, gap.after = 1, labels.cex = 1, y.quantile = 0.9999, chr.sub = TRUE, max.ranges = 1e4, axis.frac = 0.02, palette = 'BrBg', ...)
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

        sl = seqlengths(win)
        sl = sl[which(is.element(names(sl), unique(as.character(seqnames(win)))))]
        cytoband = cytoband[is.element(seqnames, names(sl))] ## ignore the chr not in win
        ## Xiaotong fixes the reordering of chromosome problem
        cytoband  = gr2dt(
            dt2gr(cytoband, seqlengths = sl) %*% win
        )[, .(seqnames = factor(as.character(seqnames), levels = names(sl)), start, end, band, stain)]
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

    ## cytoband[, seqnames := as.character(seqnames)]
    cytoband[, seqnames := as.factor(seqnames)]
    args  = list(...)
    ## some important pars
    labels.cex = ifelse(is.null(args$labels.cex), 1, args$labels.cex)
    bands.height = ifelse(is.null(args$bands.height), 0.1, args$bands.height)
    cn.height = ifelse(is.null(args$cn.height), 0.3, args$cn.height)
    link.h.ratio = ifelse(is.null(args$link.h.ratio), 0.75, args$link.h.ratio)
    bpdt = junctions$dt
    bp1 = junctions$left %>% gr2dt
    bp2 = junctions$right%>% gr2dt

    circlize::circos.clear()
    circlize::circos.par(start.degree = 90, gap.after = gap.after*1)
    circlize::circos.genomicInitialize(cytoband,
                                       sector.names = levels(cytoband$seqnames),
                                       plotType = NULL, 
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

    ## normal coverage
    ## coverage scatter plot
    if (!is.null(ncov))
    {
        if (inherits(ncov, 'data.frame'))
            ncov = dt2gr(ncov)

        ncov = ncov[!is.na(values(ncov)[[y.field]])]
        ncov = ncov[!is.infinite(values(ncov)[[y.field]])]

        if (is.na(ylim))
            ylim = c(0, quantile(values(ncov)[[y.field]], y.quantile, na.rm = TRUE))
        
        ncov$y = values(ncov)[[y.field]] %>% as.numeric
        ncov$y = ncov$y %>% pmin(ylim[2]) %>% pmax(ylim[1])

        if (is.null(ncov$col))
            ncov$col = 'black'

        ncov = ncov[sample(length(ncov), pmin(length(ncov), max.ranges))]
        uchr = unique(cytoband$seqnames)
        ncov = ncov %&% dt2gr(cytoband)
        ncovdt = gr2dt(ncov)[, seqnames := factor(seqnames, uchr)]
        circlize::circos.genomicTrackPlotRegion(
            ncovdt[, .(seqnames, start, end, y, as.character(col), ytop = y)],
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

        bpdt[is.na(lwd), lwd := 1]

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
                ## rou = circlize:::get_most_inside_radius()*c(0.1, 0.5, 1)[bpdt$span[ixs[[i]]] %>% as.integer],
                col = bpdt[ixs[[i]], ]$col,
                lwd = bpdt[ixs[[i]], ]$lwd,
                lty = bpdt[ixs[[i]], ]$lty,
                h.ratio = link.h.ratio,
                border=NA)
            )
    }
    circlize::circos.clear()
}


get_oncogenes_with_amp = function(oncotable){
    amplified_oncogenes = oncotable[grepl('ONC', role)][type == 'amp']
    strout = 'There are no oncogenes with copy number amplifications.'
    if (amplified_oncogenes[, .N] > 0){
        # TODO: we can later add the CN for each of these genes
        strout = paste0(paste(unique(amplified_oncogenes$gene), collapse = ', '), '.')
    }
    return(strout)
}

get_TSG_with_homdels = function(oncotable){
    homdel_tsgs = oncotable[grepl('TSG', role)][type == 'homdel']
    strout = 'There are no tumor suppressor genes with homozygous deletions.'
    if (homdel_tsgs[, .N] > 0){
        # TODO: we can later add the CN for each of these genes
        strout = paste0(paste(unique(homdel_tsgs$gene), collapse = ', '), '.')
    }
    return(strout)
}

check_file = function(fn, overwrite = FALSE, verbose = TRUE){
    if (file.exists(fn) & file.size(fn) > 0 & !overwrite){
        if (verbose){
            message('Found ', fn, ' so reading it. If you wish to regenerate, please repeat with "overwrite = TRUE".')
        }
        return(TRUE)
    }
        return(FALSE)
}

