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

########## load oncoKB API functions #############

#' @name get_oncokb_gene_entry_url
#' @title get_oncokb_gene_entry_url
#' @description
#'
#' Get the url for the OncoKB entry for an actionable alteration
#'
#' @param oncokb.response the JSON response, which is an output of the httr::GET query
#' @return The URL if exists, otherwise NA
#' @export
#' @author Alon Shaiber
get_oncokb_gene_entry_url = function(oncokb.response){
    if (!inherits(oncokb.response, 'list')){
        oncokb.response = list(oncokb.response)
    }
    urls = sapply(oncokb.response, function(response){
        baseurl = 'https://www.oncokb.org/gene'
        if (!inherits(response, 'response')){
            stop('You must provide an response object or a list of response objects.')
        }
        oncokb.content = httr::content(response)
        if (oncokb.content$geneExist == FALSE){
            return(NA)
        }
        gene = strsplit(oncokb.content$geneSummary, ',')[[1]][1]
        if (length(gene) == 0 | !is.character(gene)){
            return(NA)
        }
        if (!(length(oncokb.content$treatments) == 0)){
        # Notice: I am using the treatments object to get the AA mutation
        # in the future we can switch to reading it from our VCF file
            alterations = oncokb.content$treatments[[1]]$alterations
            if (is.null(alterations)){
                return(NA)
            }
            if (!inherits(alterations, 'list')){
                return(NA)
            }
            alteration = alterations[[1]]
            if (is.character(alteration) & length(alteration) > 0){
                library(RCurl)
                # if a URL exists for the specific alteration then let's go there
                url_alteration = paste0(baseurl, '/', gene, '/', alteration, '/')
                if (url.exists(url_alteration)){
                    return(url_alteration)
                }
            }
        } else {
            # If there was no URL for the specific alteration then we will try to go to the gene
            url_gene = paste0(baseurl, '/', gene, '/')
            if (url.exists(url_gene)){
                return(url_gene)
            }
        }
        return(NA)
    })
    return(urls)
}


#' @name get_oncokb_response
#' @title get_oncokb_response
#' @description
#'
#' Using httr to query genomic changes on OncoKB in accordance with the web API specification in https://www.oncokb.org/swagger-ui/index.html .
#'
#' @param variants.dt data.table with the variants to query. As a minimum, must contain the following columns: seqnames, start, end, REF, ALT.
#' @param oncokb.token a token to use when querying OncoKB. If you don't have a token, you must first obtain one from: https://www.oncokb.org/apiAccess
#' @return list of response object. Response objects can be parsed using httr (for example: httr::content).
#' @export
#' @author Alon Shaiber
get_oncokb_response = function(variants.dt, oncokb.token,
                                  oncokb.url = 'https://www.oncokb.org/api/v1/annotate',
                                  reference = 'GRCh37'){
    refdict = list('hg19' = 'GRCh37', 'hg38' = 'GRCh38')
    if (reference %in% names(refdict)){
        reference = refdict[[reference]]
    }
    if (!(reference %in% refdict)){
        stop('An invalid reference was provided. The only valid options are: ', paste(unique(c(names(refdict), refdict)), collapse = ', '))
    }
    all.res = lapply(
        1:nrow(variants.dt),
        function(i){
            res = GET(
                url = paste0(
                    oncokb.url, "/mutations/byGenomicChange?",
                    "genomicLocation=", variants.dt[i, asc(seqnames)], "%2C",
                    variants.dt[i, start], "%2C",
                    variants.dt[i, end], "%2C",
                    variants.dt[i, REF], "%2C", variants.dt[i, ALT],
                    "&referenceGenome=", reference),
                add_headers(authorization = paste("Bearer", oncokb.token)),
                add_headers(accept = "application/json")
            )
        })
}

#' @name get_oncokb_annotations
#' @title get_oncokb_annotations
#' @description
#'
#' Get a data.table containing specific fields from a list of OncoKB responses
#'
#' @param oncokb.response either a single response object or a list of response objects
#' @param fields the fields to extract from the OncoKB response
#' @return data.table with the values of the fields in each response
#' @export
#' @author Alon Shaiber
get_oncokb_annotations = function(oncokb.response,
                                  fields = c('highestSensitiveLevel',
                                             'highestResistanceLevel',
                                             'highestDiagnosticImplicationLevel',
                                             'highestPrognosticImplicationLevel')){
    if (!inherits(oncokb.response, 'list')){
        oncokb.response = list(oncokb.response)
    }
    annotations = lapply(oncokb.response, function(response){
        if (!inherits(response, 'response')){
            stop('You must provide an response object or a list of response objects.')
        }
        oncokb.content = httr::content(response)
        invalid.fields = setdiff(fields, names(oncokb.content))
        if (length(invalid.fields) > 0){
            stop('The following fields are not valid and are missing from the oncoKB response: ', invalid.fields)
        }
        vdt = setNames(data.table(matrix(ncol = length(fields), nrow = 1)), fields)
        for (field in fields){
            value = oncokb.content[[field]]
            if (length(value) > 1){
                stop('The following oncoKB field is not valid: ', field, ' since it contains more than one item.')
            }
            if (!is.null(value)){
                vdt[, (field) := value]
            }
        }
        return(vdt)
        })
    annotations = rbindlist(annotations, fill = TRUE)
    return(annotations)
}

#' @name test_oncokb_api
#' @title test_oncokb_api
#' @description
#'
#' A function to test our internal API with OncoKB
#'
#' @author Alon Shaiber
test_oncokb_api = function(oncokb.token){
    # this is an example actionable genomic alteration
    example1 = data.table(seqnames = '7', start = '140453136', end = '140453136', REF = 'A', ALT = 'T')
    # this is a gene that exists in the database, but a variant that does not
    example2 = data.table(seqnames = '2', start = '204732714', end = '204732714', REF = 'A', ALT = 'G')
    # this is a gene that is not in the database at all
    example3 = data.table(seqnames = '1', start = '69511', end = '69511', REF = 'A', ALT = 'G')
    som = rbind(example1, example2, example3)
    oncokb = get_oncokb_response(som, oncokb.token = oncokb.token)
    oncokb_annotations = get_oncokb_annotations(oncokb)
    onco_kb_entry_url = get_oncokb_gene_entry_url(oncokb)
    if (is.na(oncokb_annotations[1, highestSensitiveLevel])) stop('Something went wrong. The query for a known actionable genomic alteration did not return the expected value.')
    message('OncoKB hg19 query was finished with success.')

    # example using GRCh38
    example_GRCh38 = data.table(seqnames = 'chr7', start = '140753336', end = '140753337', REF = 'A', ALT = 'T')
    oncokb38 = get_oncokb_response(example_GRCh38, reference = 'GRCh38', oncokb.token = oncokb.token)
    oncokb_annotations = get_oncokb_annotations(oncokb38)
    onco_kb_entry_url = get_oncokb_gene_entry_url(oncokb38)
    if (is.na(oncokb_annotations[1, highestSensitiveLevel])) stop('Something went wrong. The query for a known actionable genomic alteration did not return the expected value.')
    message('OncoKB hg38 query was finished with success.')
}
########## done loading oncoKB API functions #############

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


file.good = function(f){
    file.exists(f) & (file.size(f)>0)
}



theme_pub = function(base_size=14, base_family="Helvetica") {
    (ggthemes::theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(
                    face = "plain",
                    size = rel(1.2),
                    hjust = 0.5),
                text = element_text(face = "plain"),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = "black", size = 0.5),
                panel.spacing = unit(0.1, "inch"),
                axis.title = element_text(face = "plain",size = rel(1)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(),
                axis.line = element_line(colour="black", size = 0.01),
                axis.ticks = element_line(size = 0.05),
                axis.ticks.length = grid::unit(0.02, "inch"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.position = "bottom",
                legend.direction = "horizontal",
                legend.key.size= grid::unit(0.1, "inch"),
                legend.background = element_blank(),
                plot.margin = grid::unit(c(0.2,0.2,0.2,0.2),"inch"),
                strip.background=element_blank()
                ))
}

get_gene_ampdel_annotations = function(genes_cn, amp.thresh, del.thresh){
    #' zchoo Monday, May 03, 2021 02:42:51 PM
    ## check that genes_cn is actually a data table
    if (!is.data.table(genes_cn)) {
        genes_cn = as.data.table(genes_cn)
    }
    genes_cn[, cnv := '']
    genes_cn[min_normalized_cn >= amp.thresh, cnv := 'amp']
    genes_cn[min_cn > 1 & min_normalized_cn <= del.thresh, cnv := 'del']
    genes_cn[min_cn == 1 & min_cn < ncn, cnv := 'hetdel']
    genes_cn[min_cn == 0, cnv := 'homdel']
    return(genes_cn)
}
    

#' @title check_GRanges_compatibility
#' @description
#'
#' Copmated the seqlevels of two GRanges objects. Returns TRUE/FALSE if the seqlevels are identical/non-identical and reports differences (in printed messages).
#'
#' @param gr1 GRanges object
#' @param gr2 GRanges object
#' @param name1 name of to describe gr1 (by default: 'first')
#' @param name2 name of to describe gr2 (by default: 'second')
#' @return TRUE if seqlevels are identical, otherwise, FALSE
#' @author Alon Shaiber
check_GRanges_compatibility = function(gr1, gr2, name1 = 'first', name2 = 'second'){
      # check which seqnames overlap and which don't 
      non_overlapping_seqnames1 = setdiff(seqlevels(gr1), seqlevels(gr2))
      non_overlapping_seqnames2 = setdiff(seqlevels(gr2), seqlevels(gr1))
      overlap = intersect(seqlevels(gr1), seqlevels(gr2))
      message('The following seqnames are only in the ', name1, ' GRanges, but not in the ', name2, ' GRanges: ', paste(non_overlapping_seqnames1, collapse = ', '))
      message('The following seqnames are only in the ', name2, ' GRanges, but not in the ', name1, ' GRanges: ', paste(non_overlapping_seqnames2, collapse = ', '))
      message('The follosing seqnames are in both GRanges objects: ', paste(overlap, collapse = ', '))
      if (length(non_overlapping_seqnames1) > 0 | length(non_overlapping_seqnames2) > 0){
          return(FALSE)
      }
      return(TRUE)
}
    


#' @title get_gene_copy_numbers
#' @description
#'
#' Takes a jabba_rds output and returns a GRanges with the genes that have either amplifications or deletions
#'
#' @param gg either path to rds of a gGraph or an object containing the gGraph with JaBbA output
#' @param gene_ranges GRanges of genes (must contain field "gene_name"). Alternatively a path to a file that could be parsed by rtracklayer::import (such as gtf) is acceptable.
#' @param nseg GRanges with field "ncn" - the normal copy number (if not provided then ncn = 2 is used)
#' @param nseg GRanges with field "ncn" - the normal copy number (if not provided then ncn = 2 is used)
#' @param gene_id_col the name of the column to be used in order to identify genes (must be unique for each gene, so usually "gene_name" is not the right choice).
#' @param simplify_seqnames when set to TRUE, then gr.sub is ran on the seqnames of the gGraph segments and the genes GRanges
#' @param ploidy tumor ploidy default 2
#' @param mfields the metadata fields that the output should inherit from the genes GRanges
#' @param output_type either GRanges or data.table
#' @return GRanges or data.table with genes CN
#' @author Alon Shaiber
#' @export 
get_gene_copy_numbers = function(gg, gene_ranges, nseg = NULL, gene_id_col = 'gene_id',
                                      simplify_seqnames = FALSE,
                                      mfields = c("gene_name", "source", "gene_id", "gene_type", "level", "hgnc_id", "havana_gene"),
                                      output_type = 'data.table',
                                 ploidy = 2){
    if (is.character(gg)){
      gg = readRDS(gg)
    }
    if (!inherits(gene_ranges, 'GRanges')){
        # try to import with rtracklayer
        gene_ranges = rtracklayer::import(gene_ranges)
    }

    if (!(output_type %in% c('GRanges', 'data.table'))){
        stop('Invalid output_type: ', output_type, '. outputtype must be either "GRanges" or "data.table".')
    }
    ngr = gg$nodes$gr
    if (simplify_seqnames){
        ngr = gr.sub(ngr)
        gene_ranges = gr.sub(gene_ranges)
    }
    GRanges_are_compatible = check_GRanges_compatibility(ngr, gene_ranges, 'gGraph segments', 'genes')

    if (!is.null(nseg)){
        ngr = ngr %$% nseg[, c('ncn')]
    } else {
        # if there is no nseg then assume ncn = 2
        message('No normal copy number segmentation was provided so assuming CN = 2 for all seqnames.')
        ngr$ncn = 2
    }
    ndt = gr2dt(ngr)

    seq_widths = as.numeric(width(ngr))
    # since we are comparing to CN data which is integer then we will also round the normal ploidy to the nearest integer.
    normal_ploidy = round(sum(seq_widths * ngr$ncn, na.rm = T) / sum(seq_widths, na.rm = T))

    # normalize the CN by ploidy and by local normal copy number
    ## ndt[, normalized_cn := cn * normal_ploidy / (jab$ploidy * ncn)] ## error because jab does not exist!
    ndt[, normalized_cn := cn * normal_ploidy / (ploidy * ncn)]

    # overlapping copy number segments with gene ranges
    gene_cn_segments = dt2gr(ndt, seqlengths = seqlengths(gg)) %*% gene_ranges %>% gr2dt
    # let's find genes that overlap with multiple copy number segments 
    # we would want to report the minimum and maximum CN for these genes as well as the number of CN segments overlapping the gene
    # we could do the same computation for all genes, but it is much more efficient to do it separately since the split_genes are a minority
    split_genes = gene_cn_segments[duplicated(get(gene_id_col)), get(gene_id_col)]

    gene_cn_non_split_genes = gene_cn_segments[!(get(gene_id_col) %in% split_genes)]
    gene_cn_non_split_genes[, `:=`(max_normalized_cn = normalized_cn,
                                   min_normalized_cn = normalized_cn,
                                   max_cn = cn,
                                   min_cn = cn,
                                   number_of_cn_segments = 1,
                                   cn = NULL,
                                   normalized_cn = NULL)]

    gene_cn_split_genes_min = gene_cn_segments[get(gene_id_col) %in% split_genes, .SD[which.min(cn)], by = gene_id_col]
    gene_cn_split_genes_min[, `:=`(min_normalized_cn = normalized_cn,
                                   min_cn = cn,
                                   cn = NULL,
                                   normalized_cn = NULL)]


    gene_cn_split_genes_max = gene_cn_segments[get(gene_id_col) %in% split_genes,
                                           .SD[which.max(cn)], by = gene_id_col][, .(get(gene_id_col),
                                                                max_normalized_cn = normalized_cn,
                                                                max_cn = cn)]
    setnames(gene_cn_split_genes_max, 'V1', gene_id_col)
    
    number_of_segments_per_split_gene = gene_cn_segments[get(gene_id_col) %in% split_genes, .(number_of_cn_segments = .N), by = gene_id_col]

    gene_cn_split_genes = merge(gene_cn_split_genes_min, gene_cn_split_genes_max, by = gene_id_col)
    gene_cn_split_genes = merge(gene_cn_split_genes, number_of_segments_per_split_gene, by = gene_id_col)

    gene_cn_table = rbind(gene_cn_split_genes, gene_cn_non_split_genes)

    if (output_type == 'data.table'){
        return(gene_cn_table)
    }
    return(dt2gr(gene_cn_table, seqlengths = seqlengths(gene_ranges)))
}
    
# This is under construction.
# a helper function to reduce the gencode GRanges to a a reduced version
reduce_gencode = function(gencode, gene_id_col = 'gene_id'){
    if (!inherits(gencode, 'GRanges')){
        stop('Invalid input. Input must be of class GRanges.')
    }

    genes_dt = gr2dt(gencode)
    # split to GRL according to gene_id
    genes_grl = split(gencode, gencode[, get(gene_id_col)])
    # reduce each gene to a single range. The 1e3 pad is just in order to overlap segments of the gene
    genes_grl_reduced = reduce(genes_grl + 1e3) - 1e3
    genes_gr = unlist(genes_grl_reduced)

    genes_gr[, gene_id_col] = names(genes_gr)

    new_genes_dt = gr2dt(genes_gr)

    new_genes_dt_with_metadata = merge(new_genes_dt, genes_dt[!duplicated(gene_id), ..mfields], by = 'gene_id')

    gencode_reduced = dt2gr(new_genes_dt_with_metadata)

    gencode_reduced$gene_id = as.character(gencode_reduced$gene_id)
    gencode_reduced$gene_type = as.character(gencode_reduced$gene_type)
    gencode_reduced$gene_name = as.character(gencode_reduced$gene_name)
    gencode_reduced$level = as.numeric(gencode_reduced$level)
    gencode_reduced$source = as.character(gencode_reduced$source)
    gencode_reduced$hgnc_id = as.character(gencode_reduced$hgnc_id)
    gencode_reduced$havana_gene = as.character(gencode_reduced$havana_gene)

}

#' zchoo Wednesday, May 05, 2021 11:22:17 AM

#' @name counts.to.allele.cn
#' @title counts.to.allele.cn
#'
#' @description
#'
#' allele rel2abs
#'
#' @param counts (numeric) numeric vector of counts
#' @param purity (numeric) purity
#' @param ploidy (numeric)
#'
#' @return cn (numeric) counts transformed to CN via purity + ploidy
counts.to.allele.cn = function(counts, purity, ploidy) {
    y = counts
    y.bar = 2 * mean(y, na.rm = TRUE)

    ## purity and ploidy
    alpha = purity
    tau = ploidy

    ## linear equation
    denom = alpha * tau + 2 * (1 - alpha)
    beta = (y.bar * alpha) / denom
    gamma =(y.bar * (1 - alpha)) / denom

    cn = (y - gamma) / beta
    return(cn)
}

#' @name grab.agtrack
#' @title grab.agtrack
#'
#' @description
#'
#' returns allele gtrack given sites.txt from het pileup
#'
#' @param agt.fname (character) path to sites.txt
#' @param min.frac (numeric) between 0 and 1, min frequency in normal to count as het site
#' @param max.frac (numeric) between 0 and 1, max frequency in normal to count as het site
#' @param purity (numeric)
#' @param ploidy (numeric)
#' @param major.col (character) major allele color
#' @param minor.col (character) minor allele color
#' @param max.ranges (numeric) gTrack param
#' @param max.ranges (numeric) gTrack param
#' @param lwd.border (numeric)
#' @param name (character) gTrack name
#' @param ... additional args to be passed to gTrack
#' 
#' @return allele gTrack
grab.agtrack = function(agt.fname = NULL,
                        min.frac = 0.2,
                        max.frac = 0.8,
                        purity = NULL,
                        ploidy = NULL,
                        major.col = "red",
                        minor.col = "blue",
                        max.ranges = 1e4,
                        lwd.border = 0.2,
                        name = "SNP",
                        ...) {
    if (is.null(agt.fname) || !file.exists(agt.fname)) {
        stop("agt.fname does not exist")
    }
    if (is.null(purity) | is.null(ploidy)) {
        stop("purity or ploidy not supplied")
    }
    ## prepare and filter
    agt.dt = fread(agt.fname)[alt.frac.n > min.frac & alt.frac.n < max.frac,]
    ## add major and minor
    agt.dt[, which.major := ifelse(alt.count.t > ref.count.t, "alt", "ref")]
    agt.dt[, major.count := ifelse(which.major == "alt", alt.count.t, ref.count.t)]
    agt.dt[, minor.count := ifelse(which.major == "alt", ref.count.t, alt.count.t)]

    ## melt the data frame
    agt.melted = rbind(agt.dt[, .(seqnames, start, end, count = major.count, allele = "major", col = major.col)],
                       agt.dt[, .(seqnames, start, end, count = minor.count, allele = "minor", col = minor.col)]
                       )

    ## rel2abs for alleles
    agt.melted[, cn := counts.to.allele.cn(count, purity, ploidy)]

    ## make GRanges
    agt.gr = dt2gr(agt.melted[, .(seqnames, start, end, cn, allele, col)])

    ## make gTrack
    agt = gTrack(agt.gr, y.field = "cn", name = name, 
                 ylab = "CN",
                 y.cap = FALSE, labels.suppress = TRUE,
                 circles = TRUE, lwd.border = lwd.border,
                 ...)

    return (agt)
}
