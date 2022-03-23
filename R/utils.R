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

#' @name covcbs
#' @title covcbs
#'
#' @param x (character) path to coverage GRanges
#' @param field (character) field containing t/n ratio
#' @param name (character) name of gTrack
#' @param max.ranges (numeric) maximum ranges to plot, default 1e4
#' @param lwd.border (numeric) lwd of circles in gTrack, default 0.2
#' @param purity (numeric) sample purity estimate
#' @param ploidy (numeric) sample ploidy estimate
#' @param rebin (numeric) bin size (bp) if rebinning coverage
#' @param strip.chr (logical) strip chromosome prefix of coverage file
covcbs = function(x, field = "ratio", name = "sample", max.ranges = 1e4,
                  lwd.border = 0.2, purity = NULL, ploidy = NULL, rebin = NULL,
                  strip.chr = FALSE,
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
    if (strip.chr) {
        sl = seqlengths(x)
        names(sl) = gsub(pattern = "chr", replacement = "", names(sl))
        sn = gsub(pattern = "chr", replacement = "", as.character(seqnames(x)))
        tmp = GRanges(seqnames = sn, ranges = IRanges(start = start(x), end = end(x)))
        values(tmp)[[field]] = values(x)[[field]]
        x = tmp
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

rel2abs = function(gr, purity = NA, ploidy = NA, gamma = NA, beta = NA, field = 'ratio', field.ncn = 'ncn', data_mean = NA, ncn.gr = NA, allele = FALSE, return.params = FALSE)
{
  mu = values(gr)[, field]
  mu[is.infinite(mu)] = NA
  w = as.numeric(width(gr))
  w[is.na(mu)] = NA
  sw = sum(w, na.rm = T)
  if (is.na(data_mean)){
      data_mean = sum(mu * w, na.rm = T) / sw
  }

  ncn = NA
  if (!is.null(field.ncn))
    if (field.ncn %in% names(values(gr)))
      ncn = values(gr)[, field.ncn]

  if (is.na(ncn)){
      if (!is.na(ncn.gr)){
          if (!inherits(ncn.gr, 'GRanges')){
              stop('ncn.gr must be of class GRanges, but ', class(GRanges), ' was provided.')
          }
          ncn = values(gr %$% ncn.gr[, field.ncn])[, field.ncn]
      } else {
      ncn = rep(2, length(mu))
      }
  }


  ploidy_normal = sum(w * ncn, na.rm = T) / sw  ## this will be = 2 if ncn is trivially 2

  if (allele) {
      y.bar = ploidy_normal * data_mean
      denom = purity * ploidy + ploidy_normal * (1 - purity)
      if (is.na(beta)) {
          beta = y.bar * purity / denom
      }
      if (is.na(gamma)) {
          gamma = (y.bar * (1 - purity)) / denom
      }

      if (return.params) {
          out = c(slope = 1/beta,
                  intercept = -gamma/beta)
          return(out)
      }
      
      return ((mu - gamma) / beta)
  }


  if (is.na(gamma))
    gamma = 2*(1-purity)/purity

  if (is.na(beta))
    beta = ((1-purity)*ploidy_normal + purity*ploidy) / (purity * data_mean)

  if (return.params) {
      out = c(slope = ((1-purity)*2 + purity * ploidy) / (purity * data_mean),
              intercept = -gamma)
      return (out)
  }

  return(beta * mu - ncn * gamma / 2)
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
                 cex.title = 1, mar = c(0.1,3,0.1,1), xaxs = "r", yaxs = "r", line = -5, adj = 0.5, ...) 
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
    par(mar = mar, xaxs = xaxs, yaxs = yaxs)
    if (!is.null(dim)) {
        if (length(dim) == 1) 
            dim = rep(dim, 2)
        dim = dim[1:2]
        layout(matrix(1:prod(dim), nrow = dim[1], ncol = dim[2], 
                      byrow = TRUE))
    }
    eval(expr)
    if (!is.null(title)) 
        title(title, cex.main = cex.title * max(cex), line = line, adj = adj)
    dev.off()
}



#' @name grok_vcf
#' @title modded grok_vcf
#'
#'
#' @param x path to vcf
#' @param label
#' @param keep.modifier
#' @param long
#' @param oneliner
#' @param verbose
#' @param geno
#' @param tmp.dir
#' @param gr
#' @return GRanges
#' @author Marcin Imielinski
#' @export
grok_vcf = function(x, label = NA, keep.modifier = TRUE, long = FALSE,
                    oneliner = FALSE, verbose = FALSE, geno = NULL,
                    tmp.dir = tempdir(), gr = NULL)
{
  fn = c('allele', 'annotation', 'impact', 'gene', 'gene_id', 'feature_type', 'feature_id', 'transcript_type', 'rank', 'variant.c', 'variant.p', 'cdna_pos', 'cds_pos', 'protein_pos', 'distance')

  if (is.character(x))
    {
        out = suppressWarnings(read_vcf(x, tmp.dir = tmp.dir, geno = geno, gr = gr))
        if (length(out) == 0) {
            return(out)
        }
      if (is.na(label))
        label = x
    }
  else
    out = x

  if (is.na(label))
    label = ''

  if (verbose)
    message('Grokking vcf ', label)

  if (!long)
  {
        vcf = out
        if (length(vcf)>0)
        {
          if (!is.null(vcf$ANN))
          {
            vcf$eff = unstrsplit(vcf$ANN)
            vcf$modifier = !grepl('(HIGH)|(LOW)|(MODERATE)', vcf$eff)
            if (!keep.modifier)
              vcf = vcf[!vcf$modifier]
          }
          vcf$ref = as.character(vcf$REF)
          vcf$alt = as.character(unstrsplit(vcf$ALT))
          vcf = vcf[, sapply(values(vcf), class) %in% c('factor', 'numeric', 'integer', 'logical', 'character')]
          vcf$var.id = 1:length(vcf)
          vcf$type = ifelse(nchar(vcf$ref)==nchar(vcf$alt), 'SNV',
                     ifelse(nchar(vcf$ref)<nchar(vcf$alt),
                            'INS', 'DEL'))
          vcf$label = label
        }
        return(vcf)
  }
  else if (length(out)>0)
    {
        out$REF = as.character(out$REF)
        out$ALT = as.character(unstrsplit(out$ALT))
        out$vartype = ifelse(nchar(out$REF) == nchar(out$ALT), 'SNV',
                      ifelse(nchar(out$REF) < nchar(out$ALT), 'INS', 'DEL'))
        tmp = lapply(out$ANN, function(y) do.call(rbind, lapply(strsplit(y, '\\|'), '[', 1:15)))
        tmpix = rep(1:length(out), sapply(tmp, NROW))
        meta = as.data.frame(do.call(rbind, tmp))
        colnames(meta) = fn
        meta$varid = tmpix
        meta$file = x
        out2 = out[tmpix]
        rownames(meta) = NULL
        values(out2) = cbind(values(out2), meta)
        names(out2) = NULL
        out2$ANN = NULL
        if (oneliner)
          out2$oneliner = paste(
            ifelse(!is.na(out2$gene),
                   as.character(out2$gene),
                   as.character(out2$annotation)),
            ifelse(nchar(as.character(out2$variant.p))>0,
                   as.character(out2$variant.p),
                   as.character(out2$variant.c)))
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
                if (!requireNamespace("RCurl", quietly = TRUE)) {
                      stop("Package RCurl needed.")
                }
                # if a URL exists for the specific alteration then let's go there
                url_alteration = paste0(baseurl, '/', gene, '/', alteration, '/')
                if (RCurl::url.exists(url_alteration)){
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
#' @param jabba_rds Junction object with optional metadata field  $col to specify color
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

#' @name check_file
#' @title check_file
#'
#' @description
#'
#' check whether file needs to be generated/regenerated
#' @param fn (character) path
#' @param overwrite (logical) overwrite path even if existing, default FALSE
#' @param verbose (logical)
#'
#' @return TRUE if file exists and does not need to be regenerated, otherwise FALSE
check_file = function(fn, overwrite = FALSE, verbose = TRUE){
    if (file.good(fn) & !overwrite) {
        if (verbose){
            message('Found ', fn, ' so reading it. If you wish to regenerate, please repeat with "overwrite = TRUE".')
        }
        return(TRUE)
    }
    return(FALSE)
}


file.good = function(f){
    if (is.null(f) || !is.character(f)){
        return(FALSE)
    }
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
    if ('cn.low' %in% names(genes_cn) && 'cn.high' %in% names(genes_cn)){
        genes_cn[, loh := '']
        genes_cn[is.na(cn.low), loh := NA]
        genes_cn[cn.low == 0 & ncn > 1 & min_cn > 0, loh := 'loh']
    }
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
#' @param complex.fname (character) path to complex event calls
#' @param nseg GRanges with field "ncn" - the normal copy number (if not provided then ncn = 2 is used)
#' @param nseg GRanges with field "ncn" - the normal copy number (if not provided then ncn = 2 is used)
#' @param gene_id_col the name of the column to be used in order to identify genes (must be unique for each gene, so usually "gene_name" is not the right choice).
#' @param simplify_seqnames when set to TRUE, then gr.sub is ran on the seqnames of the gGraph segments and the genes GRanges
#' @param ploidy tumor ploidy default 2
#' @param mfields the metadata fields that the output should inherit from the genes GRanges
#' @param ev.types complex event types (for now excluding simple dups/dels...)
#' @param output_type either GRanges or data.table
#' @return GRanges or data.table with genes CN
#' @author Alon Shaiber
#' @export 
get_gene_copy_numbers = function(gg, gene_ranges,
                                 complex.fname = NULL,
                                 nseg = NULL,
                                 gene_id_col = 'gene_id',
                                 simplify_seqnames = FALSE,
                                 mfields = c("gene_name", "source", "gene_id",
                                             "gene_type", "level", "hgnc_id", "havana_gene"),
                                 ev.types = c("qrp", "qpdup", "qrdel",
                                              "tic", "bfb", "dm", "chromoplexy",
                                              "chromothripsis", "tyfonas", "rigma", "pyrgo", "cpxdm"),
                                 output_type = 'data.table',
                                 ploidy = 2,
                                 verbose = TRUE){
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
    if (verbose)
        GRanges_are_compatible = check_GRanges_compatibility(ngr, gene_ranges, 'gGraph segments', 'genes')

    if (is.character(nseg)){
        if (file.exists(nseg) & endsWith(nseg, '.rds')){
            nseg = readRDS(nseg)
        }
        if (!inherits(nseg, 'GRanges')){
            nseg = NULL
        } else {
            if (verbose)
                GRanges_are_compatible = check_GRanges_compatibility(ngr, nseg, 'gGraph segments', 'nseg')
        }
    }
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
    ndt[, normalized_cn := ifelse(ncn == 0, 0, cn * normal_ploidy / (ploidy * ncn))]

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

    if ('cn.low' %in% names(ndt) && 'cn.high' %in% names(ndt)){
        # we will simply take the cn.low and cn.high from the segment with minimal CN
        gene_cn_split_genes_min = gene_cn_segments[get(gene_id_col) %in% split_genes, .SD[which.min(cn)], by = gene_id_col][,.(get(gene_id_col), cn, ncn, normalized_cn, cn.low, cn.high)]
    } else {
        gene_cn_split_genes_min = gene_cn_segments[get(gene_id_col) %in% split_genes, .SD[which.min(cn)], by = gene_id_col][,.(get(gene_id_col), cn, ncn, normalized_cn)]
    }
    gene_cn_split_genes_min[, `:=`(min_normalized_cn = normalized_cn,
                                   min_cn = cn)]
    setnames(gene_cn_split_genes_min, 'V1', gene_id_col)

    # set the ranges properly using the gene_ranges object (to guarantee that gene ranges in the output are identical to the ranges defined in gene_ranges and not some (quite arbitrary) subset of the range due to segmentation in the genome graph.
    mfields = intersect(names(values(gene_ranges)), mfields)
    if (length(mfields) == 0){
        warning('There is no overlap between the fields in the gene_ranges object and the fields provided in "mfields". We will just default to use all fields from the gene_ranges object: ', names(gene_ranges_dt))
        mfields = names(values(gene_ranges))
    }
    gene_ranges_dt = gr2dt(gene_ranges[,mfields])
    gene_cn_split_genes_min = merge.data.table(gene_ranges_dt, gene_cn_split_genes_min, by = gene_id_col) 

    gene_cn_split_genes_max = gene_cn_segments[get(gene_id_col) %in% split_genes,
                                               .SD[which.max(cn)], by = gene_id_col][, .(get(gene_id_col),
                                                                                         max_normalized_cn = normalized_cn,
                                                                                         max_cn = cn)]
    setnames(gene_cn_split_genes_max, 'V1', gene_id_col)
    
    number_of_segments_per_split_gene = gene_cn_segments[get(gene_id_col) %in% split_genes, .(number_of_cn_segments = .N), by = gene_id_col]

    gene_cn_split_genes = merge.data.table(gene_cn_split_genes_min, gene_cn_split_genes_max, by = gene_id_col)
    gene_cn_split_genes = merge.data.table(gene_cn_split_genes, number_of_segments_per_split_gene, by = gene_id_col)

    keep.fields = c('ncn', 'min_normalized_cn', 'min_cn', 'max_normalized_cn', 'max_cn', 'number_of_cn_segments')
    if ('cn.low' %in% names(ndt) && 'cn.high' %in% names(ndt)){
        keep.fields = c(keep.fields, 'cn.low', 'cn.high')
    }
    keep.fields = c(keep.fields, mfields, 'seqnames', 'start', 'end', 'strand')
    gene_cn_table = rbind(gene_cn_split_genes[, ..keep.fields], gene_cn_non_split_genes[, ..keep.fields])

    if (!is.null(gg$meta$events) && nrow(gg$meta$events)){
        this.ev = gg$meta$events[type %in% ev.types,]
        if (nrow(this.ev)>0){
            ev.grl = parse.grl(this.ev$footprint)
        values(ev.grl) = this.ev
        ev.gr = stack(ev.grl)
        } else {
            ev.gr = GRanges()
        }
    } else {
        ev.gr = GRanges()
    }
    
    ## get genes as granges
    gene_cn_gr = dt2gr(gene_cn_table, seqlengths = seqlengths(gene_ranges))

    if (length(ev.gr)){
        ov = gr.findoverlaps(gene_cn_gr, ev.gr,
                             qcol = c("gene_name"),
                             scol = c("ev.id", "type"),
                             return.type = "data.table")
    } else {
        ov = data.table()
    }
    

    if (ov[,.N] > 0){
        ov[, uid := paste0(gene_name, '_', type, '_', ev.id)]
        ov = ov[!duplicated(uid)]
        ov = ov[, .(ev.id = paste(ev.id, collapse = ","), ev.type = paste(type, collapse = ",")), by = gene_name]
        gene_cn_table = merge.data.table(gene_cn_table, ov, by = "gene_name", all.x = TRUE)
    } else {
        gene_cn_table[, ":="(ev.id = NA, type = NA)]
    }

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

    new_genes_dt_with_metadata = merge.data.table(new_genes_dt, genes_dt[!duplicated(gene_id), ..mfields], by = 'gene_id')

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
                 ylab = "CN", y0 = 0,
                 y.cap = FALSE, labels.suppress = TRUE,
                 circles = TRUE, lwd.border = lwd.border,
                 ...)

    return (agt)
}

#' @name rna_quantile
#' @title rna_quantile
#'
#' @description
#'
#' compute quantile of each transcript
#'
#' @param tpm.cohort file path data.table with id column gene and other columns representing pairs
#' @param pair character pair of interest. must be a column in tpm.cohort if tpm.pair is NULL
#' @param tpm.pair file path of data.table with column gene, pair
#'
#' @return data table with column gene and quantile representing expression quantile for the pair of interest relative to rest of cohort
rna_quantile = function(tpm.cohort, pair, tpm.pair = NULL) {
    ## this is now a required argument
    ## if (is.null(pair)) {
    ##     stop("must supply pair")
    ## }

    if (is.character(tpm.cohort) && file.good(tpm.cohort)){
        tpm.cohort = data.table::fread(tpm.cohort, header = TRUE)
    } else if (inherits(tpm.cohort, "data.table")){
        stopifnot(all(is.element("gene", colnames(tpm.cohort))))
    }

    if (!is.null(tpm.pair)) {
        if (is.character(tpm.pair) && file.good(tpm.pair)){
            tpm.pair = data.table::fread(tpm.pair, header = TRUE)
        }

        if (!inherits(tpm.pair, "data.frame")){
            stop("'tpm.pair' must be a data.frame")
        }
        tpm.pair = data.table(tpm.pair)

        ## we enforce that the sample id must be a column in the individual data matrix
        if (!is.element(pair, colnames(tpm.pair))){
            stop("The TPM value must be in the column named after the sample id.")
        }
        
        ## if sample already exists in cohort, replace with the input
        if (pair %in% colnames(tpm.cohort)) {
            tpm.cohort[[pair]] = NULL ## overwite pair if it exists
            tpm.cohort = merge.data.table(tpm.cohort, tpm.pair[, c("gene", pair), with = FALSE], by = "gene", all.x = TRUE)
        } else {
            tpm.cohort = merge.data.table(tpm.cohort, tpm.pair[, c("gene", pair), with = FALSE], by = "gene", all.x = TRUE)
        }
    }
    
    if (!(pair %in% colnames(tpm.cohort))){
        stop("pair must be in tpm.cohort columns")
    }

    id.cols = grep("gene|target_id|tumor_type|Transcript", colnames(tpm.cohort), value = TRUE)
    data.cols = setdiff(colnames(tpm.cohort), id.cols)
    melted.tpm.cohort = data.table::melt(tpm.cohort,
                                         id.vars = id.cols,
                                         measure.vars = data.cols,
                                         variable.name = "pair")

    melted.tpm.cohort[, qt := rank(as.double(.SD$value))/.N, by = gene]
    return(melted.tpm.cohort[!is.na(value)])
}

#' @name deconstructsigs_histogram
#' @title deconstructsigs_histogram
#'
#' @description
#'
#' TODO
#'
#' create histogram plot for deconstructSigs
#' 
#' @param sigs.fn (character) deconstructSigs results path
#' @param sigs.cohort.fn (character) path to cohort
#' @param pair (character)
#' @param cohort.type (character) e.g. supplied, Cell, tumor type ( actually optional )
#' @param sigMet (data.table) data table with description of mutaitonal signatures (sig.metadata.txt in the casereport data folder)
#' @param ... additional params passed ppng
#'
#' @return histogram which you can then ppng etc.
deconstructsigs_histogram = function(sigs.fn = NULL,
                                     sigs.cohort.fn = NULL,
                                     id = "",
                                     cohort.type = "",
				     outdir = "~",
			  	     sigMet = NULL,
                                     ...) {

    allsig = data.table::fread(sigs.cohort.fn)
    allsig = data.table::melt(allsig, id = "pair", variable.name = "Signature", value.name = "sig_count")

    ## counts of variants for just this pair
    var = fread(sigs.fn)[nchar(max.post)>0]
    sct = var[, .(sig_count = .N), by = .(Signature = max.post)]
    sct[, pair := id]

    ## recast as character to avoid releveling for now
    allsig[, Signature := as.character(Signature)]
    sct[, Signature := as.character(Signature)]

    if (id %in% allsig[, pair]){
        ## remove existing record from the cohort first
        allsig = allsig[pair != id]
    }

    allsig = rbind(allsig, sct)
    
    ## calculate percentile
    allsig[, perc := rank(.SD$sig_count)/.N, by = Signature]

    ## highlight the tracks where the signature is non-zero in this sample
    allsig = allsig[!is.na(Signature) & (Signature %in% sct$Signature)]

    ## only keep seqlevels with nonzero variants in this pair
    keep.slevels = allsig[pair == id, Signature] %>% as.character
    allsig = allsig[!is.na(Signature) & Signature %in% keep.slevels,]
    allsig[, Signature := gsub("Signature.", "", Signature)]

    ## reorder seqlevels by count
    new.slevels = allsig[pair == id,][order(sig_count), Signature]
    allsig[, Signature := factor(Signature, levels = new.slevels)]
	
    fwrite(allsig[pair== id,],file.path(outdir,"Sig.csv"))
    
    thisMet=sigMet[sigMet$Signature %in% allsig$Signature,]
    thisMet=thisMet[, Signature := factor(Signature, levels = new.slevels)]
    
    
    allsig=merge(allsig,thisMet,by='Signature')
    allsig$Signature_Description=paste0(allsig$MP.Summary,"\n (",as.character(allsig$Signature),")")
    allsig$Signature_Description=str_replace(allsig$Signature_Description,"@","\n")
    allsig$Signature_Description=str_replace(allsig$Signature_Description,"@temozolomide","\ntemozolomide")
    
    new.sdlevels = allsig[pair == id,][order(sig_count), Signature_Description]
    allsig[, Signature_Description := factor(Signature_Description, levels = new.sdlevels)]
	
    sigbar = ggplot(allsig, aes(y = Signature_Description, x = sig_count, fill = Signature)) +
        geom_density_ridges(bandwidth = 0.1,
                            alpha = 0.5,
                            scale = 0.9,
                            rel_min_height = 0.01,
                            color = NA,
                            jittered_points = TRUE,
                            position = position_points_jitter(width = 0.01, height = 0),
                            point_shape = '|',
                            point_size = 3,
                            point_alpha = 0.3,
                            point_colour = "black") +
        geom_segment(data = allsig[pair==id & sig_count>0],
                     aes(x = sig_count, xend = sig_count,
                         y = as.numeric(Signature), yend = as.numeric(Signature) + 0.9),
                     color = "red", size = 1, alpha = 0.5) +
        geom_label(data = allsig[pair == id],
                   aes(x = sig_count, y = as.numeric(Signature) + 0.8,
                       label = paste0("qt ", format(perc * 100, digits = 2), "%")),
                   nudge_x = 0.2,
                   hjust = "left",
                   color = "black",
                   label.size = 0,
                   alpha = 0.5) +
        scale_x_continuous(trans = "log1p",
                           breaks = c(0, 1, 10, 100, 1000, 10000, 100000),
                           labels = c(0, 1, 10,
                                      expression(10^2),
                                      expression(10^3),
                                      expression(10^4),
                                      expression(10^5)),
                           limits = c(0, 100000)) +
        labs(title = paste0("Signatures vs. ", cohort.type, " background"), x = "Burden") +
        theme_minimal() +
        theme(legend.position = "none",
              title = element_text(size = 9, family = "sans"),
              axis.title = element_text(size = 10, family = "sans"),
              axis.text.x = element_text(size = 10, family = "sans"),
              axis.text.y = element_text(size = 7, family = "sans"))
    
    return(sigbar)
}
    

#' @name rna.waterfall.plot
#' @title rna.waterfall.plot
#'
#' @description
#'
#' create waterfall plot
#' 
#' @param tpm.quantiles.fn file path data.table with column gene and tpm
#' @param rna.change.fn (character) genes with chagne in RNA expression
#' @param pair character pair of interest. must be a column in tpm.cohort if tpm.pair is NULL
#'
#' @return creates a plot
rna.waterfall.plot = function(melted.expr.fn,
                              rna.change.fn,
                              pair) {

    rna.change = fread(rna.change.fn, header = TRUE)
    if (!("gene" %in% colnames(rna.change))) {
        warning("Genes to label must contain column $gene")
        genes.to.label = c()
    } else {
        genes.to.label = rna.change[, gene]
    }

    melted.expr = fread(melted.expr.fn, header = TRUE)

    ## compute zscores
    melted.expr[, value := as.double(value)]
    melted.expr[, logvalue := log1p(value)]
    melted.expr[, m := mean(logvalue, na.rm = TRUE), by = "gene"]
    melted.expr[, v := var(logvalue, na.rm = TRUE), by = "gene"]
    melted.expr[, zs := (logvalue - m)/sqrt(v), by = "gene"]

    sel = melted.expr$pair == pair
    ds = melted.expr[sel, .(gene, zs)]

    od = order(ds$zs, decreasing = TRUE)
    flevels = ds$gene[od]
    ds[, gene := factor(gene, levels = flevels)]

    ## text for genes to label
    if (!is.null(genes.to.label)) {
        ds[, gene.label := ifelse(as.character(gene) %in% genes.to.label,
                                  as.character(gene), NA)]
        ds[gene %in% genes.to.label, role := "onc"]
        ds[gene %in% genes.to.label, label.y := ifelse(zs > 0, zs + 0.1, zs - 0.1)]
        
    } else {
        ds[, gene.label := NA]
        ds[, role := NA]
        ds[, label.y := NA]
    }
    
    pt = ggplot(ds, aes(x = gene, y = zs)) +
        geom_bar(stat = "identity", width = 1) +
        geom_label_repel(mapping = aes(label = gene.label,
                                       x = as.numeric(gene),
                                       y = label.y),
                         data = ds[!is.na(gene.label)],
                         alpha = 0.8,
                         max.overlaps = 100)

    pad = nrow(ds) * 0.05
    pt = pt +
        expand_limits(x = c(-pad, nrow(ds) + pad)) +
        labs(x = "gene", y = "z-score") +
        theme_bw() +
        theme(legend.position="none",
              axis.title.x = element_text(size = 25, family = "sans"),
              axis.title.y = element_text(size = 25, family = "sans"),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = 20, family = "sans"),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.margin=unit(c(1,1,1,1),"cm"))

    return(pt)
}



####################################################################
#' ppgrid
#'
#' least squares grid search for purity and ploidy modes
#'
#' @param segstats GRanges object of intervals with meta data fields "mean" and "sd" (i.e. output of segstats function)
#' @param allelic logical flag, if TRUE will also look for mean_high, sd_high, mean_low, sd_low variables and choose among top solutions from top copy number according to the best allelic fit
#' @param purity.min min purity value allowed
#' @param purity.max max purity value allowed
#' @param ploidy.min min ploidy value allowed
#' @param ploidy.max max ploidy value allowed
#' @param ploidy.step grid length of ploidy values
#' @param purity.step grid length of purity values
#' @param plot whether to plot the results to file
#' @param verbose print intermediate outputs
#' @param mc.cores integer number of cores to use (default 1)
#' @param return.all return the log likelihood matrix
#' @return data.frame with top purity and ploidy solutions and associated gamma and beta values, for use in downstream jbaMI
############################################
ppgrid = function(segstats,
                  allelic = FALSE,
                  purity.min = 0.01,
                  purity.max = 1.0,
                  ploidy.step = 0.01,
                  purity.step = 0.01,
                  ploidy.min = 1.2, # ploidy bounds (can be generous)
                  ploidy.max = 6,
                  plot = F,
                  verbose = F,
                  mc.cores = 1,
                  return.nll = FALSE){
    if (verbose)
        jmessage('setting up ppgrid matrices .. \n')

    if (is.na(ploidy.min)) ploidy.min = 1.2
    if (is.na(ploidy.max)) ploidy.max = 6
    if (is.na(purity.min)) purity.min = 0.01
    if (is.na(purity.max)) purity.max = 1

    ##  purity.guesses = seq(0, 1, purity.step)
    purity.guesses = seq(pmax(0, purity.min), pmin(1.00, purity.max), purity.step)
    ## ploidy.guesses = seq(pmin(0.5, ploidy.min), pmax(10, ploidy.max), ploidy.step)
    ploidy.guesses = seq(pmax(0.5, ploidy.min), pmax(0.5, ploidy.max), ploidy.step)

    if (allelic)
        if (!all(c('mean_high', 'mean_low', 'sd_high', 'sd_low') %in% names(values(segstats))))
        {
            jwarning('If allelic = TRUE then must have meta data fields mean_high, mean_low, sd_high, sd_low in input segstats')
            allelic = FALSE
        }

    if (is.null(segstats$mean))
        jerror('segstats must have field $mean')

    segstats = segstats[!is.na(segstats$mean) & !is.na(segstats$sd)]

    if (!is.null(segstats$ncn))
        segstats = segstats[segstats$ncn==2, ]

    ## if (is.null(segstats$ncn))
    ##     ncn = rep(2, length(mu))
    ## else
    ##     ncn = segstats$ncn

    if (any(tmpix <-is.infinite(segstats$mean) | is.infinite(segstats$sd)))
    {
        segstats$sd[tmpix] = segstats$mean[tmpix] = NA
    }
    segstats = segstats[!is.na(segstats$mean) & !is.na(segstats$sd), ]
    if (length(segstats)==0)
        jerror('No non NA segments provided')

    mu = segstats$mean
    w = as.numeric(width(segstats))
    Sw = sum(as.numeric(width(segstats)))
    sd = segstats$sd
    m0 = sum(as.numeric(mu*w))/Sw

    if (verbose)
        cat(paste(c(rep('.', length(purity.guesses)), '\n'), collapse = ''))

    NLL = matrix(unlist(parallel::mclapply(seq_along(purity.guesses), function(i)
    {
        if (verbose)
            cat('.')
        nll = rep(NA, length(ploidy.guesses))
        for (j in seq_along(ploidy.guesses))
        {
            alpha = purity.guesses[i]
            tau = ploidy.guesses[j]
            gamma = 2/alpha - 2
            beta = (tau + gamma)/m0 
            v = pmax(0, round(beta*mu-gamma))
            nll[j] = sum((v-beta*mu+gamma)^2/((sd)^2))
        }
        return(nll)
    }, mc.cores = mc.cores)), nrow = length(purity.guesses), byrow = T)

    dimnames(NLL) = list(as.character(purity.guesses), as.character(ploidy.guesses))

    if (verbose)
        cat('\n')

    ## rix = as.numeric(rownames(NLL))>=purity.min & as.numeric(rownames(NLL))<=purity.max
    ## cix = as.numeric(colnames(NLL))>=ploidy.min & as.numeric(colnames(NLL))<=ploidy.max
    ## NLL = NLL[rix, cix, drop = FALSE]

    a = rep(NA, nrow(NLL));
    b = rep(NA, ncol(NLL)+2)
    b.inf = rep(Inf, ncol(NLL)+2)
                                        #  a = rep(Inf, nrow(NLL));
                                        #  b = rep(Inf, ncol(NLL)+2)
    NLLc = rbind(b, cbind(a, NLL, a), b) ## padded NLL and all of its shifts
    NLLul = rbind(cbind(NLL, a, a), b.inf, b)
    NLLuc = rbind(cbind(a, NLL, a), b.inf, b)
    NLLur = rbind(cbind(a, a, NLL), b.inf, b)
    NLLcl = rbind(b, cbind(NLL, a, a), b)
    NLLcr = rbind(b, cbind(a, a, NLL), b)
    NLLll = rbind(b, b, cbind(NLL, a, a))
    NLLlc = rbind(b, b, cbind(a, NLL, a))
    NLLlr = rbind(b, b, cbind(a, a, NLL))

    if (min(c(ncol(NLL), nrow(NLL)))>1) ## up up down down left right left right ba ba start
        M = (NLLc < NLLul &
             NLLc < NLLuc &
             NLLc < NLLur &
             NLLc < NLLcl &
             NLLc < NLLcr &
             NLLc < NLLll &
             NLLc < NLLlc &
             NLLc < NLLlr)[-c(1, nrow(NLLc)),
                           -c(1, ncol(NLLc)),
                           drop = FALSE]
    else if (ncol(NLL)==1) ## one column, only go up and down
        M = (NLLc < NLLuc & NLLc < NLLlc)[-c(1, nrow(NLLc)), -c(1, ncol(NLLc)), drop = FALSE]
    else ## only row, only go left right
        M = (NLLc < NLLcl & NLLc < NLLcr)[-c(1, nrow(NLLc)), -c(1, ncol(NLLc)), drop = FALSE]

    if (length(M)>1)
    {
        ix = Matrix::which(M, arr.ind=T);
        if (nrow(ix)>1)
        {
            C = hclust(d = dist(ix), method = 'single')
            cl = cutree(C, h = min(c(nrow(NLL), ncol(NLL), 2)))
            minima = ix[vaggregate(1:nrow(ix), by = list(cl), function(x) x[which.min(NLL[ix[x, drop = FALSE]])]), , drop = FALSE]
        }
        else
            minima = ix[1,, drop = FALSE]
    }
    else
        minima = cbind(1,1)

    out = data.frame(purity = as.numeric(rownames(NLL)[minima[,1]]), ploidy = as.numeric(colnames(NLL)[minima[,2]]), NLL = NLL[minima],
                     i = minima[,1], j = minima[,2])

    out = out[order(out$NLL), , drop = FALSE]
    rownames(out) = 1:nrow(out)
    ## Saturday, Sep 02, 2017 10:33:26 PM
    ## Noted floating point error, use the epsilon trick to replace '>='
    ## out = out[out$purity>=purity.min & out$purity<=purity.max & out$ploidy>=ploidy.min & out$ploidy<=ploidy.max, ]
    eps = 1e9
    out = out[out$purity - purity.min >= -eps &
              out$purity - purity.max <= eps &
              out$ploidy - ploidy.min >= -eps &
              out$ploidy - ploidy.max <= eps, ]
    out$gamma = 2/out$purity -2
    out$beta = (out$ploidy + out$gamma)/m0
    out$mincn = mapply(function(gamma, beta) min(round(beta*mu-gamma)), out$gamma, out$beta)
    out$maxcn = mapply(function(gamma, beta) max(round(beta*mu-gamma)), out$gamma, out$beta)

    ## group solutions with (nearly the same) slope (i.e. 1/beta), these should have almost identical
    ## NLL (also take into account in-list distance just be safe)
    if (nrow(out)>1)
        out$group = cutree(hclust(d = dist(cbind(100/out$beta, 1:nrow(out)), method = 'manhattan'), method = 'single'), h = 2)
    else
        out$group = 1
    out = out[out$group<=3, ,drop = FALSE] ## only pick top 3 groups

    if (allelic) ## if allelic then use allelic distance to rank best solution in group
    {
        ## remove all NA allelic samples
        segstats = segstats[!is.na(segstats$mean_high) & !is.na(segstats$sd_high) & !is.na(segstats$mean_low) & !is.na(segstats$sd_low)]
        out$NLL.allelic = NA
        mu = cbind(segstats$mean_high, segstats$mean_low)
        w = matrix(rep(as.numeric(width(segstats)), 2), ncol = 2, byrow = TRUE)
        Sw = sum(as.numeric(width(segstats)))*2
        sd = cbind(segstats$sd_high, segstats$sd_low)
        m0 = sum(as.numeric(mu*w))/Sw

        if (verbose)
            cat(paste(c(rep('.', length(purity.guesses)), '\n'), collapse = ''))

        for (i in 1:nrow(out))
        {
            if (verbose)
            {
                jmessage(sprintf('Evaluating alleles for solution %s of %s\n', i, nrow(out)))
            }
            alpha = out$purity[i]
            tau = out$ploidy[i]
                                        #                  gamma = 2/alpha - 2
            gamma = 1/alpha - 1 ## 1 since we are looking at hets
            beta = (tau + gamma)/m0 ## replaced with below 9/10/14
                                        #          beta = ( tau + tau_normal * gamma /2 ) / m0
                                        #          v = pmax(0, round(beta*mu-ncn*gamma/2))
            v = pmax(0, round(beta*mu-gamma))

            vtot = round(out$beta[i]*segstats$mean-out$gamma[i])
            vlow.mle = rep(NA, length(vtot))

            for (j in seq_along(vlow.mle))
            {
                if (vtot[j]==0)
                    vlow.mle[j] = 0
                else
                {
                    vlow = 0:floor(vtot[j]/2)
                    vhigh = vtot[j]-vlow
                    tmp.nll = cbind((vlow-beta*mu[j,2]+gamma)^2/(sd[j,2])^2, (vhigh-beta*mu[j, 1]+gamma)^2/((sd[j,1])^2))
                    vlow.mle[j] = vlow[which.min(rowSums(tmp.nll))]
                }
            }

            vlow.mle = apply(cbind(mu, sd, vtot), 1, function(x) {
                tot = x[5]
                if (tot == 0)
                    return(0)
                else
                {
                    vlow = 0:floor(tot/2)
                    vhigh = tot-vlow
                    muh = x[1]
                    mul = x[2]
                    sdh = x[3]
                    sdl = x[4]
                    tmp.nll = cbind((vlow-beta*mul+gamma)^2/(sdl)^2, (vhigh-beta*muh+gamma)^2/((sdh)^2))
                    return(vlow[which.min(rowSums(tmp.nll))])
                }
            })

            out$NLL.allelic[i] = sum((cbind(vtot-vlow.mle, vlow.mle)-beta*mu+gamma)^2/sd^2)
        }

        out$NLL.tot = out$NLL
        out$NLL = out$NLL.tot + out$NLL.allelic
        out.all = out
        ix = vaggregate(1:nrow(out), by = list(out$group), FUN = function(x) x[order(abs(out$NLL[x]))][1])
    } else ## otherwise choose the one that gives the lowest magnitude copy number
    {
        out.all = out
        ix = vaggregate(1:nrow(out), by = list(out$group), FUN = function(x) x[order(abs(out$mincn[x]), out$mincn[x]<0)][1])
    }

    out = out[ix, , drop = FALSE]
    out$NLL = vaggregate(out$NLL, by = list(out$group), FUN = min)

    out.all$keep = 1:nrow(out.all) %in% ix ## keep track of other ploidy group peaks for drawing purposes
    out.all = out.all[out.all$group %in% out$group, ] ## only draw the groups in the top solution
    out = out.all;
    out = out[order(out$group, !out$keep, out$NLL), ]
    out$rank = NA
    out$rank[out$keep] = 1:sum(out$keep)
    out$keep = out$i = out$j = NULL
    rownames(out) = NULL
    if (return.nll){
        out = list(out = out, NLL = NLL)
    }
    return(out)
}


## diagnostic function used by karyograph_stub
#' @name ggplot_ppfit
#' @rdname internal
ggplot_ppfit = function(seg, purity, ploidy, beta = NULL, gamma = NULL, ploidy_normal = 2, xlim = c(-Inf, Inf)){
    tmp = seg  ## only plot seg that we haven't fixed SD for and that have normal cn 1, to minimize confusion
    dupval = sort(table(tmp$mean), decreasing = TRUE)[1]
    if (!is.na(dupval))
        if (dupval>5)
            tmp = tmp[-which(as.character(tmp$mean) == names(dupval))]

    if (length(tmp)==0)
        return()

    ## sampling random loci to plot not segments
    ## segsamp = pmin(sample(tmp$mean, 1e6, replace = T, prob = width(tmp)), xlim[2])
    tmp = tmp %Q% (!is.na(mean))
    segsamp = sample(tmp[, "mean"], 1e6, replace = T, prob = width(tmp)) %>% gr2dt
    segsamp[, mean := pmax(xlim[1], pmin(mean, xlim[2]))]
    segsamp = segsamp[!is.na(mean)]

    mu = seg$mean
    w = width(seg)
    Sw = sum(as.numeric(width(seg)))
    sd = seg$sd
    m0 = sum(as.numeric(mu*w))/Sw
    alpha = purity
    tau = ploidy
    ncn = rep(2, length(mu))
    
    beta = (tau + ploidy_normal * gamma / 2 ) / m0
    gamma = 2/alpha - 2
    grids = 1/beta*(0:1000) + gamma/beta
    grids = grids[which(grids < segsamp[, max(mean)])]

    
    h = ggplot(segsamp, aes(x = mean)) +
        geom_histogram(bins = 1000) +
        geom_vline(xintercept = grids, lty = 2, color = alpha("red", 0.8), lwd = 0.5) +
        ## scale_x_continuous(limits = c(0, segsamp[, quantile(mean, 0.99)])) +
        labs(x = "Segment intensity",
             title = sprintf('Purity: %.2f Ploidy: %.2f Beta: %.2f Gamma: %.2f', purity, ploidy, beta, gamma)) +
        theme_minimal(12)
    
    return(h)
    ## ppdf(print(h))
    
    ## hist(
    ##     pmax(xlim[1], pmin(xlim[2], segsamp)),
    ##     1000, xlab = 'Segment intensity',
    ##     main = sprintf('Purity: %s Ploidy: %s Beta: %s Gamma: %s', kag$purity, kag$ploidy, round(kag$beta,2), round(kag$gamma,2)),
    ##     xlim = c(pmax(0, xlim[1]), pmin(xlim[2], max(segsamp, na.rm = T))))
    ## abline(v = 1/kag$beta*(0:1000) + kag$gamma/kag$beta, col = alpha('red', 0.3), lty = c(4, rep(2, 1000)))
}

#' @name rna_reformat
#' @title rna_reformat
#'
#' @description
#'
#' preprocess kallisto/salmon raw tsv files for downstream analysis
#' 
#' @param kallisto.fname (character) path to kallisto/salmon output
#' @param pair (character) id
#' @param gngt.fname (character) gencode gTrack file name
#'
#' @return data.table with columns gene, pair
rna_reformat = function(kallisto.fname,
                        pair = NULL,
                        gngt.fname = NULL) {

    ## check files all valid
    if (is.null(kallisto.fname) || !file.exists(kallisto.fname)) {
        stop("kallisto.fname not valid")
    }

    if (is.null(pair)) {
        stop("pair cannot be NULL")
    }

    ## grab correct gene for target id
    if (is.character(gngt.fname) && file.good(gngt.fname)){
        ge.data = stack(readRDS(gngt.fname)@data[[1]])
    } else if (inherits(gngt.fname, "GRanges")){
        ge.data = gngt.fname
    } else {
        stop("Invalid genocde gTrack prov. Please check your --gencode argument.")
    }

    kallisto.dt = fread(kallisto.fname, header = TRUE)

    ## check if empty
    if (nrow(kallisto.dt) == 0) {
        out.dt = data.table(gene = c(), tmp = c())
        return(setnames(out.dt, "tmp", pair))
    }

    ## check if file is already ok
    outcols = c("gene", pair)
    if (all(outcols %in% colnames(kallisto.dt))) {
        return (kallisto.dt[, ..outcols])
    }

    ## check if target_id and tpm are in kallisto
    salmon.reqcols = c("Name", "TPM") ## required columns in Salmon output
    kallisto.reqcols = c("target_id", "tpm") ## required columns in kallisto output
    ## browser()
    if (all(salmon.reqcols %in% colnames(kallisto.dt))) {

        ## preprocess salmon input
        message("RNA input type: Salmon")

        outcols = c("target_id", outcols)

        gene.symbols = unlist(lapply(strsplit(kallisto.dt[, Name], "\\|"), function(x) {x[[5]]}))
        kallisto.dt[, gene := gene.symbols]
        kallisto.dt = kallisto.dt[TPM > 0 & !is.na(TPM) & !is.na(gene),
                                  .(TPM = max(TPM)),
                                  by = gene]
        setnames(kallisto.dt, "TPM", pair)
        kallisto.dt[, target_id := gene]
    } else if (all(kallisto.reqcols %in% colnames(kallisto.dt))) {

        ## preprocess Kallisto input
        message("RNA input type: Kallisto")

        outcols = c("target_id", outcols)

        kallisto.dt[, gene := ge.data$gene_name[match(target_id, ge.data$transcript_id)]]

        ## subset for high expression and max transcript for that gene
        ## not perfect, but potentially helps to prevent double-counting co-expressed exons
        kallisto.dt = kallisto.dt[tpm > 0 & !is.na(tpm) & !is.na(gene),][, .(tpm = max(tpm), target_id = target_id[which.max(tpm)]), by = gene]

        ## reset tpm column
        setnames(kallisto.dt, "tpm", pair)
        
    } else {
        stop("required columns target_id and tpm missing from kallisto output")
    }

    return(kallisto.dt[!is.na(gene), ..outcols])
}

#' @name grab.gene.ranges
#' @title grab.gene.ranges
#'
#' @param gngt.fname (character) gencode gTrack file name
#' @param genes (character) vector of gene symbols
#'
#' @return GRanges with metadata column gene_name
grab.gene.ranges = function(gngt.fname, genes = as.character()) {

    if (!file.exists(gngt.fname)) {
        stop("Gencode gTrack does not exist")
    }
    if (is.null(genes) || length(genes) == 0) {
        stop("genes empty")
    }

    gngt = readRDS(gngt.fname)

    ge.stacked = gngt@data[[1]] %>% range %>% stack
    drivers.ge.dt = gr2dt(ge.stacked %Q% (name %in% genes))
    drivers.ge.dt[, gene_name := name]
    drivers.gr = dt2gr(drivers.ge.dt)

    return(drivers.gr)
}


#' @name filter.snpeff
#' @title filter.snpeff
#'
#' @description
#'
#' Overlaps annotated VCF with genes provided in delimited file and returns results as data.table
#' If VCF doesn't exist or is not provided returns an empty data.table
#' 
#' @param vcf (character) vcf file name
#' @param gngt.fname (character) gencode gTrack file name
#' @param onc (character) path to .rds file with character vector of oncogenes
#' @param tsg (character) path to .rds file with character vector of TSGs
#' @param ref.name (character) one of hg19 or hg38
#' @param verbose (logical) default FALSE
#'
#' @return data.table with columns gene, seqnames, pos, REF, ALT, variant.p, vartype, annotation
filter.snpeff = function(vcf,
                         gngt.fname,
                         onc, tsg, drivers.fname, ref.name = "hg19", verbose = FALSE,
                         type = NULL) {

    if (NROW(type) == 0) {
        grab_snv = FALSE
        grab_indel = FALSE
    } else {
        type = type[type %in% c("snv", "indel")]
        if (length(type) > 1) {
            message("multiple types provided: ", paste(type, collapse = ", "))
            message("using: ", type[1])
            type = type[1]
        }
        grab_snv = "snv" == type
        grab_indel = "indel" == type
    }


    ## parsing drivers
    txt.ptrn = "(.tsv|.txt|.csv|.tab)(.xz|.bz2|.gz){0,}$"
    no_ext <- function (x, compression = FALSE) 
    {
        if (compression) 
            x <- sub("[.](gz|bz2|xz)$", "", x)
        sub("([^.]+)\\.[[:alnum:]]+$", "\\1", x)
    }
    if (check_file(drivers.fname)) {
        if (grepl(".rds$", drivers.fname, ignore.case = TRUE)) {
            drivers = readRDS(drivers.fname)
            if (!(is.character(drivers) ||
                  inherits(drivers, c("matrix", "data.frame"))))
                drivers = NULL
            
        } else if (grepl(txt.ptrn, drivers.fname)) {
            drivers = fread(drivers.fname)
        }

        if (!is.null(drivers)) {

            if (NCOL(drivers) == 1) {
                ## annotating drivers by file name
                drivers = as.data.frame(drivers)
            }
        }
        drivers = drivers[[1]]
    } else {
        drivers = NULL
    }
    
    if (verbose) message("Querying for ", type)

    dummy.out = data.table(
        gene = character(),
        seqnames = character(),
        pos = numeric(),
        REF = character(),
        ALT = character(),
        variant.p = character(),
        vartype = character(),
        impact = character(),
        annotation = character()
    )
    
    if (is.null(vcf) || is.na(vcf) || !file.exists(vcf)) {
        warning("vcf file missing")
        return(dummy.out)
    }
    if (!file.exists(gngt.fname)) {
        stop("gencode gtrack file missing")
    }
    if (!ref.name %in% c("hg19", "hg38")) {
        stop("invalid ref.name")
    }

    if (verbose) {
        message("Reading VCF/BCF input")
    }

    if (grepl("bcf", vcf)) {
        vcf.gr = grok_bcf(vcf, long = TRUE, indel = grab_indel, snv = grab_snv)
    } else if (grepl("vcf", vcf)) {
        vcf.gr = grok_vcf(vcf, long = TRUE)
        if (grab_snv) {
            vcf.gr = vcf.gr[vcf.gr$vartype %in% "SNV"]
        } else if (grab_indel) {
            vcf.gr = vcf.gr[vcf.gr$vartype %in% c("INS", "DEL")]
        }
    } else {
        stop("Invalid file type")
    }

    if (!("impact" %in% names(values(vcf.gr)))) {
        warning("Missing impact column! Returning all variants")
    } else {
        vcf.gr = vcf.gr %Q% (!impact %in% c("LOW", "MODIFIER"))
    }

    if (length(vcf.gr) == 0) {
        if (verbose) {
            message("No high impact variants")
        }
        return (dummy.out)
    }

    vcf.gr = vcf.gr %Q%
        (!duplicated(paste(CHROM, POS)))

    ## if (length(vcf.gr) == 0) {
    ##     if (verbose) {
    ##         message("No high impact variants")
    ##     }
    ##     return (dummy.out)
    ## }

    if (verbose) {
        message("Found ", length(vcf.gr), " variants, overlapping with genes")
    }

    genes = c(readRDS(onc), readRDS(tsg), drivers)

    if (length(genes) == 0) {
        if (verbose) {
            message("No high impact variants")
        }
        return (dummy.out)
    }

    genes.gr = grab.gene.ranges(gngt.fname, genes)
    vcf.gr = vcf.gr %*% genes.gr

    if (length(vcf.gr) == 0) {

        if (verbose) {
            message("No SNPs in cancer genes")
        }
        return (dummy.out)
    }

    vcf.dt = as.data.table(vcf.gr)[, .(gene = gene_name,
                                       seqnames, pos = start, impact,
                                       REF, ALT, variant.p, vartype, annotation)]
    return(unique(vcf.dt))
}

#' @name create.summary
#' @title create.summary
#'
#' @description
#'
#' calculate stuff for summary table
#'
#' @param jabba_rds (character) jabba rds
#' @param snv_vcf (character) path to snv vcf/bcf
#' @param indel_vcf (character) path to indel vcf/bcf
#' @param verbose (logical) default FALSE
#' @param chrs (character vector) seqnames default 1:22, x, y
#' @param amp.thresh (numeric) default 4
#' @param del.thresh (numeric) default 0.5
#'
#' @return list with names
#' - purity
#' - ploidy
#' - junction_burden
#' - del_mbp
#' - amp_mbp
#' - cna_mbp
#' - cna_frac
#' - mut_count
#' - mut_per_mbp
#' @export
create.summary = function(jabba_rds,
                          snv_vcf,
                          indel_vcf,
                          verbose = FALSE,
                          chrs = c(as.character(1:22), "X", "Y"),
                          amp.thresh = 4,
                          del.thresh = 0.5) {

    jab = readRDS(jabba_rds)
    gg = gG(jabba = jab)

    ## get standard chromosomes
    chrs = grep("(^(chr)*[0-9XY]+$)", unique(gg$nodes$dt$seqnames), value = TRUE)

    ## segments as data table
    segs.dt = gg$nodes$dt[as.character(seqnames) %in% chrs, ]

    out = list(purity = jab$purity,
               ploidy = jab$ploidy,
               junction_burden = length(gg$junctions[type == "ALT" & cn > 0]))

    if (verbose) {
        message("Computing total width of deleted and amplified segments...")
    }
    
    #out$del_mbp = sum(segs.dt[cn <= (del.thresh * jab$ploidy), width], na.rm = TRUE) / 1e6
    out$del_mbp = sum(segs.dt[cn <= (0.5 * jab$ploidy), width], na.rm = TRUE) / 1e6
    #out$amp_mbp = sum(segs.dt[cn >= (amp.thresh * jab$ploidy), width], na.rm = TRUE) / 1e6
    out$amp_mbp = sum(segs.dt[cn >= (1.5 * jab$ploidy), width], na.rm = TRUE) / 1e6
    out$cna_mbp = out$del_mbp + out$amp_mbp

    if (verbose) {
        message("Mbp affected by CNA: ", out$cna_mbp)
    }

    total_width = sum(segs.dt[, width], na.rm = TRUE) / 1e6
    out$cna_frac = out$cna_mbp / total_width

    if (verbose) {
        message("Fraction affected by CNA: ", out$cna_frac)
    }

    if (file.good(snv_vcf)) {
        snv.vcf = readVcf(file = snv_vcf)
        if ("FILTER" %in% names(values(snv.vcf))) {
            snv.count = sum(values(snv.vcf)$FILTER == "PASS", na.rm = TRUE)
        } else {
            snv.count = length(snv.vcf)
        }
    } else {
        snv.count = 0
    }

    if (file.good(indel_vcf)) {
        indel.vcf = readVcf(file = indel_vcf)
        if ("FILTER" %in% names(values(indel.vcf))) {
            indel.count = sum(values(indel.vcf)$FILTER == "PASS", na.rm = TRUE)
        } else {
            indel.count = length(indel.vcf)
        }
    } else {
        indel.count = 0
    }

    out$mut_count = snv.count + indel.count
    out$mut_per_mbp = out$mut_count / total_width

    if (verbose) {
        message("Number of indels/SNVs: ", out$mut_count)
        message("Number of indels/SNVs per mbp: ", out$mut_per_mbp)
    }

    return(out)
}

#' @name grab.hets
#' @title grab.hets
#'
#' @description
#'
#' returns allele gtrack given sites.txt from het pileup
#'
#' @param agt.fname (character) path to sites.txt
#' @param min.frac (numeric) between 0 and 1, min frequency in normal to count as het site
#' @param max.frac (numeric) between 0 and 1, max frequency in normal to count as het site
#' 
#' @return allele gTrack
grab.hets = function(agt.fname = NULL,
                     min.frac = 0.2,
                     max.frac = 0.8)
{
    if (is.null(agt.fname) || !file.exists(agt.fname)) {
        stop("agt.fname does not exist")
    }

    ## prepare and filter
    agt.dt = fread(agt.fname)[alt.frac.n > min.frac & alt.frac.n < max.frac,]
    ## add major and minor
    agt.dt[, which.major := ifelse(alt.count.t > ref.count.t, "alt", "ref")]
    agt.dt[, major.count := ifelse(which.major == "alt", alt.count.t, ref.count.t)]
    agt.dt[, minor.count := ifelse(which.major == "alt", ref.count.t, alt.count.t)]

    ## melt the data frame
    agt.melted = rbind(agt.dt[, .(seqnames, start, end, count = major.count, allele = "major")],
                       agt.dt[, .(seqnames, start, end, count = minor.count, allele = "minor")]
                       )

    ## make GRanges
    agt.gr = dt2gr(agt.melted[, .(seqnames, start, end, count, allele)])

    return (agt.gr)
}


#' @name pp_plot
#' @title pp_plot
#'
#' @details
#'
#' create histogram of estimated copy number given purity and ploidy
#'
#' @param jabba_rds (character) JaBbA output
#' @param cov.fname (character) path to coverage GRanges (supply if allele = FALSE)
#' @param hets.fname (character) path to sites.txt (supply if allele = TRUE)
#' @param allele (logical) allelic CN? default FALSE
#' @param field (character) default ratio if allele is FALSE and count if allele is TRUE
#' @param plot.min (numeric) minimum CN default -2
#' @param plot.max (numeric) max CN (factor times ploidy) default 2
#' @param bins (numeric) number of histogram bins default 500
#' @param scatter (logical) default FALSE
#' @param height (numeric) plot height
#' @param width (numeric) plot width
#' @param output.fname (character) path of output directory
#' @param purity (numeric) if not provided then the value is read from the input JaBbA
#' @param ploidy (numeric) if not provided then the value is read from the input JaBbA
#' @param verbose (logical)
#'
#' @return output.fname
#' @export
pp_plot = function(jabba_rds = NULL,
                   cov.fname = NULL,
                   hets.fname = NULL,
                   allele = FALSE,
                   field = NULL,
                   plot.min = -2,
                   plot.max = 2,
                   scatter = FALSE,
                   bins = 500,
                   height = 800,
                   width = 800,
                   output.fname = "./plot.png",
                   purity = NA,
                   ploidy = NA,
                   verbose = FALSE) {

    if (is.null(jabba_rds) || !file.exists(jabba_rds)) {
        stop("jabba_rds does not exist")
    }
    jab = readRDS(jabba_rds)
    purity = ifelse(!is.na(purity) && is.numeric(purity), purity, jab$purity)
    ploidy = ifelse(!is.na(ploidy) && is.numeric(ploidy), ploidy, jab$ploidy)
    if (!allele) {
        if (is.null(cov.fname) || !file.exists(cov.fname)) {
            stop("cov.fname not supplied and allele = TRUE")
        }
        if (!grepl(pattern = "rds", x = cov.fname)) {
            stop("cov.fname must be .rds file containing GRanges object")
        }
        cov = readRDS(cov.fname)
        if (!inherits(cov, "GRanges")) {
            stop("cov is not GRanges")
        }
        if (length(cov) == 0) {
            stop("empty GRanges")
        }
        if (is.null(field)) {
            field = "ratio"
        }
        if (!(field %in% names(values(cov)))) {
            stop("cov missing required field")
        }
        if (verbose) {
            message("Grabbing coverage and converting rel2abs")
        }
        cov$cn = rel2abs(cov, field = field, purity = purity, ploidy = ploidy, allele = FALSE)
        ## get mean CN over JaBbA segments
        if (verbose) {
            message("computing mean over jabba segments")
        }
        segs = gr.stripstrand(jab$segstats %Q% (strand(jab$segstats)=="+"))[, c()]
        segs = gr.val(segs, cov[, "cn"], val = "cn", mean = TRUE, na.rm = TRUE)
        if (verbose) {
            message("tiling")
        }
        tiles = gr.tile(gr = segs, width = 1e4)
        tiles = gr.val(tiles, segs[, "cn"], val = "cn", mean = TRUE, na.rm = TRUE)
        if (verbose) {
            message("Grabbing transformation slope and intercept")
        }
        eqn = rel2abs(cov, field = field, purity = purity, ploidy = ploidy, allele = FALSE, return.params = TRUE)
        dt = as.data.table(tiles)
    } else {
        if (is.null(hets.fname) || !file.exists(hets.fname)) {
            stop("hets.fname not supplied")
        }
        hets = grab.hets(hets.fname)
        if (is.null(field)) {
            field = "count"
        }
        if (!field %in% names(values(hets))) {
            stop("hets missing required field")
        }
        if (verbose) {
            message("Grabbing hets and converting rel2abs")
        }
        hets$cn = rel2abs(hets, field = field, purity = purity, ploidy = ploidy, allele = TRUE)
        eqn = rel2abs(hets, field = field, purity = purity, ploidy = ploidy, allele = TRUE, return.params = TRUE)
        if (verbose) {
            message("computing mean over jabba segments")
        }
        segs = gr.stripstrand(jab$segstats %Q% (strand(jab$segstats)=="+"))[, c()]
        major.segs = gr.val(segs, hets %Q% (allele == "major"), val = "cn", mean = TRUE, na.rm = TRUE)
        minor.segs = gr.val(segs, hets %Q% (allele == "minor"), val = "cn", mean = TRUE, na.rm = TRUE)
        if (verbose) {
            message("Tiling")
        }
        tiles = gr.tile(gr = segs, width = 1e4)
        major.tiles = gr.val(tiles, major.segs, val = "cn", mean = TRUE, na.rm = TRUE)
        minor.tiles = gr.val(tiles, minor.segs, val = "cn", mean = TRUE, na.rm = TRUE)
        dt = rbind(as.data.table(major.tiles)[, .(seqnames, start, end, allele = "major", cn)],
                   as.data.table(minor.tiles)[, .(seqnames, start, end, allele = "minor", cn)])
    }

    maxval = plot.max * ploidy # max dosage
    minval = plot.min ## min dosage

    ## remove things with weird ploidy
    dt = dt[cn < maxval & cn > minval & grepl("[0-9]", seqnames)==TRUE]

    if (verbose) {
        message("Making plot for ", nrow(dt), " points")
    }
    
    if (!allele) {

        pt = ggplot(dt, aes(x = cn)) +
            geom_histogram(fill = "gray", bins = bins, alpha = 0.8) +
            scale_x_continuous(breaks = 0:floor(maxval),
                               labels = 0:floor(maxval) %>% as.character,
                               sec.axis = sec_axis(trans = ~(. - eqn["intercept"])/eqn["slope"],
                                                   name = field)) +
            geom_vline(xintercept = 0:floor(maxval), color = "red", linetype = "longdash") +
            labs(x = "Estimated CN", y = "count") +
            theme_bw() +
            theme(legend.position = "none",
                  axis.title = element_text(size = 20, family = "sans"),
                  axis.text.x = element_text(size = 20, family = "sans"),
                  axis.text.y = element_text(size = 14, family = "sans"))

    } else {

        if (scatter) {

            dt = cbind(as.data.table(major.tiles)[, .(seqnames, start, end, major.cn = cn)],
                       as.data.table(minor.tiles)[, .(minor.cn = cn)])
            dt = dt[major.cn < maxval & minor.cn < maxval &
                    major.cn > minval & minor.cn > minval &
                    grepl("[0-9]", seqnames)==TRUE,]

            pt = ggplot(dt, aes(x = major.cn, y = minor.cn)) +
                geom_point(size = 2, alpha = 0.1, color = "gray") +
                scale_x_continuous(breaks = 0:floor(maxval),
                                   labels = 0:floor(maxval) %>% as.character,
                                   sec.axis = sec_axis(trans = ~(. - eqn["intercept"])/eqn["slope"],
                                                       name = paste("Major", field))) +
                scale_y_continuous(breaks = 0:floor(maxval),
                                   labels = 0:floor(maxval) %>% as.character,
                                   sec.axis = sec_axis(trans = ~(. - eqn["intercept"])/eqn["slope"],
                                                       name = paste("Minor", field))) +
                labs(x = "Major CN", y = "Minor CN") +
                theme_bw() +
                theme(legend.position = "none",
                      axis.title = element_text(size = 20, family = "sans"),
                      axis.text.x = element_text(size = 20, family = "sans"),
                      axis.text.y = element_text(size = 14, family = "sans"))

            pt = ggMarginal(pt, type = "histogram",
                            xparams = list(bins = bins),
                            yparams = list(bins = bins))
            
        } else {

            pt = ggplot(dt, aes(x = cn)) +
                geom_histogram(fill = "gray", bins = bins, alpha = 0.8) +
                scale_x_continuous(breaks = 0:floor(maxval),
                                   labels = 0:floor(maxval) %>% as.character,
                                   sec.axis = sec_axis(trans = ~(. - eqn["intercept"])/eqn["slope"],
                                                       name = field)) +
                geom_vline(xintercept = 0:floor(maxval), color = "red", linetype = "longdash") +
                labs(x = "Estimated CN", y = "count") +
                facet_grid(row = vars(allele)) +
                theme_bw() +
                theme(legend.position = "none",
                      axis.title = element_text(size = 20, family = "sans"),
                      axis.text.x = element_text(size = 20, family = "sans"),
                      axis.text.y = element_text(size = 14, family = "sans"),
                      strip.text.y = element_text(size = 20, family = "sans"))
        }

    }

    return(pt) ##ppng(print(pt), filename = normalizePath(output.fname), height = height, width = width)
}

#' @name oncotable
#' @title oncotable
#' @description
#'
#' Very slightly modified version of oncotable copied from skitools
#' 
#' Takes as input (keyed) "tumors" (aka pairs) table which a metadata table with specific
#' columns pointing to paths corresponding to one or more of the following pipeline outputs:
#'
#' $annotated_bcf  Path to annotated.bcf file that is the primary output of SnpEff module from which TMB and basic mutation
#' descriptions are extracted along with their basic categories (these will comprising the core of the oncoplot are computed)
#' 
#' $fusions  Path to fusion.rds file that is the primary output of the Fusions modjle, from which protein coding fusions will
#' be computed for
#' 
#' $jabba_rds  Path to jabba.simple.rds output representing primary output of JaBbA module from which SCNA and overall
#' junction burden are computed
#' 
#' $complex    Path to complex.rds gGnome cached object that is the primary output of Events module, from which simple
#' and complex event burdens are computed
#' 
#' $signature_counts Path to signature_counts.txt that is the primary output of Signatures module from which SNV signature
#' counts are computed
#' 
#' $hrd_results Path to the output of HRDetect, from which the basic parameters of the HRDetect model are taken from.
#'
#' The function then outputs a melted data.table of "interesting" features that can be saved and/or immediately output
#' into oncoprint.  This data.table will at the very least have fields $id $type (event type), $track, and  $source
#' populated in addition to a few other data type specific columns.
#'
#' The $source column is the name of the column of tumors from which that data was extracted, and track is a grouping
#' variable that allows separation of the various data types. 
#'
#' All the paths above can be NA or non existent, in which case a dummy row is inserted into the table so that downstream
#' applications know that data is missing for that sample. 
#'
#' @param tumors keyed data.table i.e. keyed by unique tumor id with specific columns corresponding to  paths to pipeline outputs(see description)
#' @param gencode path to gencode .gtf or .rds with GRanges object, or a GRanges object i.e. resulting from importing the (appropriate) GENCODE .gtf via rtracklayer, note: this input is only used in CNA to gene mapping
#' @param amp.thresh SCNA amplification threshold to call an amp as a function of ploidy (4)
#' @param del.thresh SCNA deletion threshold for (het) del as a function of ploidy (by default cn = 1 will be called del, but this allows additoinal regions in high ploidy tumors to be considered het dels)
#' @param max.reldist (numeric) maximum relativedistance to be considered enh hijacking, default 0.25
#' @param min.refdist (numeric) minimum reference distance to be considered enh hijacking, default 5e6
#' @param mc.cores number of cores for multithreading
#' @param verbose logical flag 
#' @author Marcin Imielinski
#' @export
oncotable = function(tumors, gencode = NULL, verbose = TRUE,
                     amp.thresh = 4,
                     filter = 'PASS',
                     del.thresh = 0.5,
                     max.reldist = 0.25,
                     min.refdist = 5e6,
                     mc.cores = 1) {
    require(skitools)
    .oncotable = function(dat, x = dat[[key(dat)]][1], pge = NULL, verbose = TRUE, amp.thresh = 4, del.thresh = 0.5, filter = 'PASS', max.reldist = 0.25, min.refdist = 5e6)
    {
        out = data.table()

        ## collect RNA data if it exists
        if (file.good(dat[x, rna])) {
            if (verbose) {
                message("pulling RNA for ", x)
            }
            rna.dt = fread(dat[x, rna])##readRDS(dat[x, rna])
            if (nrow(rna.dt)) {
                rna = rna.dt[, .(gene, value, role, quantile = qt, id = x,
                                 type = ifelse(grepl("over", direction, ignore.case = TRUE), "overexpression",
                                        ifelse(grepl("under", direction, ignore.case = TRUE),
                                               "underexpression",
                                               NA_character_)),
                                 track = "expression", source = "rna")]
                out = rbind(out, rna, fill = TRUE, use.names = TRUE)
            } else {
                out = rbind(out, data.table(id = x, type = NA, source = "rna"), fill = TRUE, use.names = TRUE)
            }
        } else {
            out = rbind(out, data.table(id = x, type = NA, source = "rna"), fill = TRUE, use.names = TRUE)
        }

        ## collect gene fusions
        if (file.good(dat$fusions))
        {
            if (verbose)
                message('pulling $fusions for ', x)
            fus = readRDS(dat[x, fusions])$meta
            if (nrow(fus))
            {
                fus = fus[silent == FALSE, ][!duplicated(genes), ]
                fus[, vartype := ifelse(in.frame == TRUE, 'fusion', 'outframe_fusion')] # annotate out of frame fusions
                fus = fus[, .(gene = strsplit(genes, ',') %>% unlist, vartype = rep(vartype, sapply(strsplit(genes, ','), length)))][, id := x][, track := 'variants'][, type := vartype][, source := 'fusions']
                out = rbind(out, fus, fill = TRUE, use.names = TRUE)

            }
        } 
        else ## signal missing result
            out = rbind(out, data.table(id = x, type = NA, source = 'fusions'), fill = TRUE, use.names = TRUE)

        ## collect complex events
        if (file.good(dat$complex))
        {
            if (verbose)
                message('pulling $complex events for ', x)
            sv = readRDS(dat[x, complex])$meta$events
            if (nrow(sv))
            {
                ## modified to keep individual complex SV's separate
                sv = sv[, .(value = .N), by = type][, id := x][, track := ifelse(type %in% c('del', 'dup', 'invdup', 'tra', 'inv'), 'simple sv', 'complex sv')][, source := 'complex']
                ## sel.cols = intersect(c("type", "ev.id", "footprint"), colnames(sv))
                ## sv = sv[, ..sel.cols][, ev.type := type][, id := x][, track := ifelse(type %in% c('del', 'dup', 'invdup', 'tra', 'inv'), 'simple sv', 'complex sv')][, source := 'complex']
                out = rbind(out, sv, fill = TRUE, use.names = TRUE)
            }
        }
        else
            out = rbind(out, data.table(id = x, type = NA, source = 'complex', track = ""), fill = TRUE, use.names = TRUE)

        ## collect copy number / jabba
        if (!is.null(dat$jabba_rds) && file.exists(dat[x, jabba_rds]))
        {
            if (verbose)
                message('pulling $jabba_rds to get SCNA and purity / ploidy for ', x)
            jab = readRDS(dat[x, jabba_rds])
            out = rbind(out,
                        data.table(id = x, value = c(jab$purity, jab$ploidy), type = c('purity', 'ploidy'), track = 'pp'),
                        fill = TRUE, use.names = TRUE)

            ## grabs SCNA data from pre-computed table
            if (!is.null(dat$scna) && file.exists(dat[x, scna])) {
                scna.dt = readRDS(dat[x, scna])

                ## subset for previously annotated variants if data table is nonempty
                if (nrow(scna.dt) && "cnv" %in% colnames(scna.dt)) {
                    if ('loh' %in% colnames(scna.dt)){
                        scna.dt = scna.dt[cnv %in% c("amp", "del", "homdel", "hetdel") | loh == TRUE, ]
                    } else {
                        scna.dt = scna.dt[cnv %in% c("amp", "del", "homdel", "hetdel"),]
                        scna.dt[, loh := NA]
                    }
                }

                ## if there are any CN variants, rbind them to existing output
                if (nrow(scna.dt)) {
                    sel.cols = intersect(c("gene_name", "gene", "cnv", "loh",
                                           "min_cn", "cn.low", "cn.high", "min_normalized_cn", "max_cn", "max_normalized_cn",
                                           "seqnames", "start", "end", "ncn", "gene_id"),
                                         colnames(scna.dt))
                    sel.cols = intersect(sel.cols, names(scna.dt))
                    scna = scna.dt[, ..sel.cols]

                    ## make sure that gene is named as gene instead of gene_name
                    if ("gene_name" %in% colnames(scna)) {
                        setnames(scna, "gene_name", "gene")
                    }

                    ## make sure ncn copy number is present
                    if (!("ncn" %in% colnames(scna))) {
                        scna[, ncn := 2]
                    }

                    ## set 'type' to cnv annotation if present
                    if ("cnv" %in% colnames(scna)) {
                        # TODO: currently the type will not include loh info
                        #       so if we want the oncoprint to show this information in the future we will need to make some adjustments
                        setnames(scna, "cnv", "type")
                    } else {
                        scna[min_normalized_cn >= amp.thresh, type := 'amp']
                        scna[min_cn > 1 & min_normalized_cn <= del.thresh, type := 'del']
                        scna[min_cn == 1 & min_cn < ncn, type := 'hetdel']
                        scna[min_cn == 0, type := 'homdel']
                    }
                    
                    out = rbind(out, scna[, ":="(id = x,
                                                 track = "variants",
                                                 vartype = "scna",
                                                 source = "jabba_rds")],
                                fill = TRUE, use.names = TRUE)

                    # add an loh vartype
                    out = rbind(out, scna[, ":="(id = x,
                                                 track = "variants",
                                                 type = loh,
                                                 vartype = "loh",
                                                 source = "jabba_rds")],
                                fill = TRUE, use.names = TRUE)
                } else {
                    if (verbose) {
                        message("No SCNAs found")
                    }
                }
            }
        } else {
            out = rbind(out, data.table(id = x, type = NA, source = 'jabba_rds', track = "variants"), fill = TRUE, use.names = TRUE)
        }

        if ('proximity' %in% names(dat) && file.good(dat[x, proximity]) && nrow(readRDS(dat[x, proximity])$dt)) {
            if (verbose)
                message('Processing proximity results.')
            proximity.dt = readRDS(dat[x, proximity])$dt[reldist < max.reldist & refdist > min.refdist,]
            if (proximity.dt[,.N] > 0){
                out = rbind(out, proximity.dt[, .(gene = gene_name, value = reldist,
                                                  reldist, altdist, refdist, walk.id,
                                                  type = "proximity", source = "proximity", track = "proximity")],
                            use.names = TRUE,
                            fill = TRUE)
                if (verbose)
                    message('Adding ', proximity.dt[,.N], ' proximity entries.')
            } else {
                if (verbose)
                    message('No proximity hits found')
            }
        } else {
            out = rbind(out, data.table(id = x, type = NA, source = "proximity", track = "proximity"), fill = TRUE, use.names = TRUE)
        }


        if (file.good(dat[x, deconstruct_variants])) {
            if (verbose)
                message('pulling $deconstruct_variants for ', x)
            sig = fread(dat[x, deconstruct_variants])
            sig.dt = sig[, .(value = .N), by = max.post][max.post %like% "Signature"]
            setnames(sig.dt, "max.post", "type")
            sig.dt[, ":="(id = x,
                          etiology = NA,
                          frac = value / sum(value, na.rm = T),
                          track = "signature",
                          source = "deconstruct_variants")]
            out = rbind(out, sig.dt, fill = TRUE, use.names = TRUE)
        } else {
            out = rbind(out, data.table(id = x, type = NA, source = 'deconstruct_sigs', track = "signature"), fill = TRUE, use.names = TRUE)
        }

        ## collect gene mutations for SNVs and optionally indels
        if (file.good(dat$annotated_snv_bcf)) {
            if (verbose)
                message('pulling $annotated_bcf for ', x, ' using FILTER=', filter)
            vars = annotated_bcf_to_oncotable(dat$annotated_snv_bcf,
                                              filter = filter,
                                              id = x,
                                              verbose = verbose)
            out = rbind(out, vars, fill = TRUE, use.names = TRUE)
        } else {
            out = rbind(out,
                        data.table(id = x, type = NA, source = 'annotated_bcf', track = "variants"),
                        fill = TRUE, use.names = TRUE)
        }
        
        if (file.good(dat$annotated_indel_bcf)) {
            if (verbose)
                message('pulling $annotated_indel_bcf for ', x, ' using FILTER=', filter)
            vars = annotated_bcf_to_oncotable(dat$annotated_indel_bcf,
                                              filter = filter,
                                              id = x,
                                              verbose = verbose)
            out = rbind(out, vars, fill = TRUE, use.names = TRUE)
        } else {
            if (verbose) {
                message("$annotated_indel_bcf not supplied")
            }
            out = rbind(out,
                        data.table(id = x, type = NA, source = 'annotated_bcf', track = "variants"),
                        fill = TRUE, use.names = TRUE)
        }

	    ## adding hrd results
    	if ( !is.null(dat$hrd_results) && file.good(dat[x, hrd_results])){
	    if(verbose){
	    	message('adding HRDetect results for ',x)
		}
	hrd.res = readRDS(dat[x, hrd_results])
        hrd.out = as.data.table(hrd.res$hrdetect_output) %>% data.table::melt()
   	colnames(hrd.out)[1]="vartype"
	hrd.out$source=rep("HRDetect",nrow(hrd.out))
	hrd.out$vartype=as.character(hrd.out$vartype)
	hrd.out$vartype[6]="hrd_index"
	hrd.out$vartype[8]="hrd_probability"
	hrd.out$type=c('HRD_Intercept','microdel_prop','SNV3_sig','SV3_sig','SV5_sig','hrd_indx','SNV8_sig','hrdprob_BRCA')
	hrd.out$track=rep("hrd_res",nrow(hrd.out))
	out=rbind(out, hrd.out, fill = TRUE, use.names = TRUE)
    }
	else{
		out=rbind(out, data.table(id = x, type = NA, source = 'HRDetect'), fill = TRUE, use.names = TRUE)
	}

        if (verbose)
            message('done ', x)

        return(out)
    }

    if (is.null(key(tumors)))
    {
        if (is.null(tumors$id))
            stop('Input tumors table must be keyed or have column $id')
        else
            setkey(tumors, id)
    }

    out = mclapply(tumors[[key(tumors)]], .oncotable,
                   dat = tumors, amp.thresh = amp.thresh, filter = filter, del.thresh = del.thresh, verbose = verbose, mc.cores = mc.cores, max.reldist = max.reldist, min.refdist = min.refdist)
    out = rbindlist(out, fill = TRUE, use.names = TRUE)

    setnames(out, 'id', key(tumors))

    ## make sure no NA tracks or source
    out[is.na(source), source := ""]
    out[is.na(track), track := ""]
    
    return(out)
}

#' @name oncoprint
#' @title oncoprint
#' @description
#' 
#' Slighly modified version of oncoprint function from skitools.
#'
#' @param tumors  keyed table of tumors (aka pairs table) with field $oncotable which points to a cached .rds file of an oncotable e.g. produced by oncotable function or Oncotable module / task
#' @param oncotab output from oncotable function with field $id
#' @param genes character vector of genes
#' @param columns additional columns of tumors matrix to plot as horizontal tracks below the main track
#' @param split character of name of column in tumors table to split on (NULL)
#' @param sort  logical flag whether to automatically sort rows i.e. genes and columns i.e. tumors in a "stair step" pattern or default to the provided (TRUE)
#' @param noncoding logical flag whether to show non protein coding mutations
#' @param sort.genes logical flag whether to sort rows i.e. genes with respect to their frequency (TRUE)
#' @param sort.tumors logical flag whether to sort columns i.e. patients in a stairstep pattern with respect to the provided gene order (TRUE)
#' @param sv.stack  logical flag whether to stack bar plot simple and complex SV event counts (FALSE)
#' @param signatures logical flag whether to show signatures (if data is provided / available) (TRUE)
#' @param svevents logical flag whether to show events (if data is provided / available) (TRUE)
#' @param tmb logical flag whether to show TMB bar plot (TRUE)
#' @param tmb.log  logical flag whether to log TMB + 1 (TRUE)
#' @param pp logical flag whether to show purity / ploidy (if data is provided / available) (TRUE)
#' @param ppdf whether to print to pdf via ppdf
#' @param track.height height of tracks in 'cm'
#' @param split.gap  gap between splits
#' @param signature.main integer indices of main COSMIC signatures to keep
#' @param signature.thresh lower threshold for non main signature fraction in at least one sample to plot
#' @param outframe.fusions show fusions that are out-of-frame (FALSE)
#' @param hrd_res logical flag whether to show HRD results (FALSE)
#' @param cex length 1 or 2 vector canvas expansion factor to apply to the oncoprint itself (relative to 10 x 10 cm) (c(1,3))
#' @param return.mat whether to return.mat
#' @param wes logical flag whether to use wesanderson coolors
#' @param mc.cores multicore threads to use for $oncotable loading from tumors table (not relevant if oncotab provided)
#' @param ... other arguments to ppdf
#' @return ComplexHeatmap object (if ppdf = FALSE), and genotype matrix (if return)
#' @author Marcin Imielinski
#' @export 
oncoprint = function(tumors = NULL,
                     oncotab = NULL,
                     genes = c('KRAS', 'EGFR', 'BRAF', 'TP53', 'TERT', 'CCND1', 'MYC', 'PIK3CA', 'PTEN', 'CDKN2A', 'ARID1A', 'SMARCA4'),
                     split = NULL, 
                     sort = TRUE, sort.genes = sort, sort.tumors = sort,
                     columns = NULL,
                     noncoding = FALSE,
                     cna = TRUE, tmb = TRUE, pp = TRUE, signature = TRUE, hrd_res=TRUE, svevents = TRUE, basic = FALSE, 
                     ppdf = TRUE,
                     return.oncotab = FALSE,
                     return.mat = FALSE,                     
                     wes = TRUE,
                     drop = TRUE,
                     drop.genes = FALSE, 
                     track.height = 1,
                     signature.thresh = 0.2,
                     signature.main = c(1:5,7,9,13),
                     outframe.fusions = FALSE,
                     track.gap = track.height/2,
                     split.gap = 1,
                     colnames.fontsize = 10,
                     rownames.fontsize = 10,
                     track.fontsize = 10,
                     mc.cores = 1,
                     verbose = FALSE,
                     height = 20,
                     width = 20,
                     ...)
{
  require(skitools)

  if (basic)
    tmb = svevents = signature = FALSE

  if (!length(genes))
    stop('genes must be provided either as a vector or named list of gene identifiers')

  if (is.list(genes))
    genes = dunlist(genes)[, .(genes = V1, group = listid)]
  else
    genes = data.table(genes = genes, group = NA)

  genes = genes[!duplicated(genes), ]

  if (!is.null(tumors))
  {
    if (!is.null(key(tumors)))
      tumors$id = tumors[[key(tumors)]]

    if (is.null(tumors$id))
      stop('tumors be either keyed or have $id field, if you are resorting e.g. manually sorting your input table the key may get lost so then you should set an $id field explicitly')
    
    if (any(duplicated(tumors$id)))
      stop('check key field in tumors table: duplicated ids present. The key should be unique per row, and matched to the $id field of oncotab')
  }

  missing = c()
  if (is.null(oncotab))
  {
    errmsg = 'Either oncotab or tumors argument must be provided, where tumors is a keyed data table (where each row is a tumor) with column $oncotable of file paths pointing to the cached rds Oncotable results for each tumors'
    if (is.null(tumors) || is.null(tumors$oncotable))
      stop(errmsg)

    fe = file.exists(tumors$oncotable)
    missing = union(missing, tumors$id[!fe])

    if (any(!fe))
      warning(paste(sum(!fe), 'of', length(fe), 'tumors with missing oncotab, will remove if drop = TRUE, otherwise mark'))

    if (!nrow(tumors))
      stop('No tumors with $oncotable field pointing to existing path')

    if (verbose)
      message('Scraping $oncotable paths for oncotable .rds files.  To speed up, consider multi-threading with mc.cores and if you will be creating multiple plots.  Also consider running this with return.oncotab = TRUE and use that for subsequent calls via oncotab = argument.')

    oncotab = mclapply(which(fe), function(x) {y = readRDS(tumors$oncotable[x]); if (nrow(y)) y[, id := tumors$id[x]]; return(y)}, mc.cores = mc.cores) %>% rbindlist(fill = TRUE)
    oncotab$id = factor(oncotab$id, tumors$id)    
  }

  if (!is.null(tumors))
    {
      oncotab$id = factor(oncotab$id, tumors$id)
      missing = union(missing, setdiff(tumors$id, oncotab$id))
    }
  else
    oncotab$id = factor(oncotab$id)
  
  oncotab = oncotab[!is.na(id), ]

  if (!nrow(oncotab))
  {
    if (!is.null(tumors))
      stop('empty oncotable provided, check tumors table, there may be an id mismatch or no non empty files')
    else
      stop('empty oncotable provided, please check inputs')
  }

  vars = oncotab[track == 'variants', ][gene %in% genes$genes, ][type != 'synonymous', ]

  ## keep track of missing samples ie those that had either SNV, jabba, fusions
  ## will get a gray column in the plot
  missing = union(missing, vars[track == 'variants' & is.na(type), id])

  if (!noncoding)
    vars = vars[!(type %in% c('promoter', 'noncoding', 'regulatory')), ]

  if (!cna)
    vars = vars[!(type %in% c('amp', 'del', 'hetdel', 'homdel')), ]  

  vars[, gene := factor(gene, genes$genes)]
  vars = vars[!is.na(gene), ]

  ## convert to matrix format for complex heatmap
  if (nrow(vars))
    {
      varc = dcast.data.table(data = vars, gene ~ id, value.var = "type", fill = 'WT', drop = FALSE, fun.aggregate = function(x) paste(x, collapse = ','))
      varm = as.matrix(varc[, -1])
      rownames(varm) = varc$gene
    }
  else
  {
    varm = matrix('WT', nrow = length(levels(vars$gene)), ncol = length(levels(vars$id)),
           dimnames = list(levels(vars$gene), levels(vars$id)))
  }

  ## prune / label missing genotypes (ie either due to missing or incomplete oncotable entries)
  if (length(missing))
    {
      if (!drop)
        varm[, intersect(colnames(varm), missing)] = 'missing'
      else
        varm = varm[, setdiff(colnames(varm), missing)]
    }

  ## then gene binary order
  if (sort.genes)
    {
      ##ix = skitools::border(varm!='') %>% rev
      ix = rev(order(rowSums(varm!='WT' & varm != 'missing', na.rm = TRUE)))
      varm = varm[ix, , drop = FALSE]
    }
  
  ## then sample binary mutation order
  if (sort.tumors)
    {
      jx = rev(skitools::border(t(varm)!='WT' & t(varm) != 'missing'))
      varm = varm[, jx, drop = FALSE]
    }
    
  ## customize appeagrid appearance with mix of rectangles and circles
  ord = c("amp", 'del', "hetdel", "homdel", 'trunc', 'splice', 'inframe_indel', 'fusion', 'missense', 'promoter', 'regulatory')
  if (outframe.fusions == TRUE){
      ord = c("amp", 'del', "hetdel", "homdel", 'trunc', 'splice', 'inframe_indel', 'outframe_fusion', 'fusion', 'missense', 'promoter', 'regulatory')
  }
  alter_fun = function(x, y, w, h, v) {
    CSIZE = 0.25
    w = convertWidth(w, "cm")*0.7
    h = convertHeight(h, "cm")*0.7
    l = min(unit.c(w, h))
    grid.rect(x, y, w, h, gp = gpar(fill = alpha("grey90", 0.4), col = NA))
    v = v[ord]
    for (i in which(v)) {
      if (names(v)[i] %in% c('amp', 'del', "hetdel", "homdel", 'fusion', 'outframe_fusion'))
        grid.rect(x,y,w,h, gp = gpar(fill = varcol[names(v)[i]], col = NA))
      else if (grepl("missing", names(v)[i]))
        grid.rect(x, y, w, h, gp = gpar(fill = varcol[names(v)[i]], col = NA))
      else if (grepl("trunc", names(v)[i]))
        {
          grid.segments(x - w*0.5, y - h*0.5, x + w*0.5, y + h*0.5,
                        gp = gpar(lwd = 2, col = varcol[names(v)[i]]))
          grid.segments(x - w*0.5, y + h*0.5, x + w*0.5, y - h*0.5,
                        gp = gpar(lwd = 2, col = varcol[names(v)[i]]))
        }
      else if (grepl("(missense)|(promoter)|(regulatory)", names(v)[i]))
      {
        grid.circle(x,y,l*CSIZE, gp = gpar(fill = varcol[names(v)[i]], col = NA))
      }
      else {
        if (grepl("indel", names(v)[i]))
          grid.rect(x,y,w*0.9,h*0.4, gp = gpar(fill = varcol[names(v)[i]], col = NA))
      }
    }
  }

  varcol = c(
    WT = alpha('gray', 0),
    fusion = alpha('green', 0.5),
    outframe_fusion = alpha('greenyellow', 0.5),
    hetdel = 'lightblue',
    missing = 'gray',            
    amp = "red",
    drop = FALSE,
    homdel = "darkblue",
    del = 'cyan',
    missense = 'gray40',
    inframe_indel = 'darkgreen',
    promoter  = alpha('red', 0.5),
    regulatory  = alpha('red', 0.2),
    trunc = alpha("blue", 0.8),
    mir = alpha('purple', 0.4),
    splice = "purple"
  )
  
  ids = colnames(varm)
  out.mat = varm ## in case we want to return.mat

  ## generate additional plots if requested / available
  bottom_data = top_data = list()
  if (tmb & any(oncotab$track == 'tmb'))
  {
    tmbd = oncotab[track == 'tmb' & type == 'density', structure(value, names = as.character(id))][ids]
    
    if (!all(tmbd == 0)){
      top_data$TMB = tmbd
      out.mat = rbind(TMB = tmbd, out.mat)
    }
  }

  if (pp & any(oncotab$track == 'pp'))
  {
    top_data$Purity = oncotab[track == 'pp' & type == 'purity', structure(value, names = as.character(id))][ids]
    top_data$Ploidy = oncotab[track == 'pp' & type == 'ploidy', structure(value, names = as.character(id))][ids]

    out.mat = rbind(Purity = top_data$Purity, Ploidy = top_data$Ploidy, out.mat)
  }

  if (hrd_res & any(oncotab$track == 'hrd_res'))
  {
    top_data$Microhomology_Deletion_Proportion = oncotab[track == 'hrd_res' & type == 'microdel_prop', structure(value, names = as.character(id))][ids]
    top_data$SNV3_Score = oncotab[track == 'hrd_res' & type == 'SNV3_sig', structure(value, names = as.character(id))][ids]
    top_data$SV3_Score = oncotab[track == 'hrd_res' & type == 'SV3_sig', structure(value, names = as.character(id))][ids]
    top_data$SV5_Score = oncotab[track == 'hrd_res' & type == 'SV5_sig', structure(value, names = as.character(id))][ids]
    top_data$SNV8_Score = oncotab[track == 'hrd_res' & type == 'SNV8_sig', structure(value, names = as.character(id))][ids]
    top_data$HRD_Index = oncotab[track == 'hrd_res' & type == 'hrd_indx', structure(value, names = as.character(id))][ids]
    top_data$HRD_Probability_BRCA = oncotab[track == 'hrd_res' & type == 'hrdprob_BRCA', structure(value, names = as.character(id))][ids]

    out.mat = rbind(HRD_Index = top_data$HRD_Index, HRD_Probability_BRCA = top_data$HRD_Probability_BRCA,  SNV3_Score = top_data$SNV3_Score, SV3_Score = top_data$SV3_Score, SV5_Score = top_data$SV5_Score, SNV9_Score= top_data$SNV8_Score, out.mat)
  }
	
  sink("~/test.txt")
  print(out.mat)
  sink()
  ## put together top track from all topdata
  ab = anno_oncoprint_barplot(border = FALSE, height = unit(track.height, "cm"))                
  toptracks = HeatmapAnnotation(column_barplot = ab)
  if (length(top_data))
  {
    topcols = brewer.master(names(top_data), wes = wes)
    tmp = lapply(names(top_data),
                 function(x) anno_barplot(top_data[[x]],
                                          border = FALSE,
                                          axis_param = list(gp = gpar(fontsize = track.fontsize)),
                                          height = unit(track.height, 'cm'),
                                          gp = gpar(fill = topcols[x], col = topcols[x])))
    names(tmp) = names(top_data)
    tmp$gap = unit(track.gap, 'cm')
    toptracks = do.call(HeatmapAnnotation, c(tmp, list(column_barplot = ab)))
  }

  packed_legends = list()
  bottomtracks = list()
  if (isTRUE(signature) && any(oncotab$track == 'signature') && 'frac' %in% names(oncotab))
  {
    sigd = oncotab[track == 'signature', ][type != 'Residual', ]

    ## keep any signature outside of keep that has at least signature.thresh in at least
    ## one tumor
    signature.keep = paste('Signature', signature.main, sep = '_') %>%
      union(sigd[frac>signature.thresh, type])
    sigd[, type := ifelse(type %in% signature.keep, as.character(gsub('Signature_', '', type)), 'other')]
    sigdc = dcast.data.table(sigd, id ~ type, value.var = 'frac', fun.aggregate = sum, drop = FALSE)
    sigdm = as.matrix(sigdc[, -1])
    rownames(sigdm) = sigdc$id
    sigdm = sigdm[ids,, drop = FALSE]
    sigdm = sigdm[, suppressWarnings(order(as.numeric(colnames(sigdm)))), drop = FALSE]
    out.mat = rbind(out.mat, t(sigdm))
    if (wes)
      sigcols = brewer.master(colnames(sigdm), 'BottleRocket1', wes = TRUE)
    else
      sigcols = brewer.master(colnames(sigdm), 'Dark2')

    sigcols['other'] = 'gray'
    bottomtracks$COSMIC = anno_barplot(
      sigdm,
      legend = TRUE,
      axis_param = list(gp = gpar(fontsize = track.fontsize)),
      height = unit(3*track.height, 'cm'),
      border = FALSE,
      gp = gpar(fill = sigcols, col = sigcols)
    )
    packed_legends = c(packed_legends,
      list(Legend(labels = names(sigcols), ncol = 2, legend_gp = gpar(fill = sigcols), title = 'COSMIC')))
  }

  if (svevents & any(oncotab$track %in% c('complex sv', 'simple sv')))
  {
    cx = dcast.data.table(oncotab[track == 'complex sv', ][, type := as.character(type)][, id := factor(id, ids)], id ~ type, fill = 0, drop = FALSE, value.var = 'value')
    simple = dcast.data.table(oncotab[track == 'simple sv', ][, type := as.character(type)][, id := factor(id, ids)], id ~ type, fill = 0, drop = FALSE, value.var = 'value')
    out.mat = rbind(out.mat, t(as.matrix(cx[,-1])), t(as.matrix(simple[,-1])))

    uev = names(cx)[-1]
    if (wes)
    {
      cxcols = brewer.master(names(cx)[-1], 'IsleOfDogs1', wes = TRUE)
      simplecols = brewer.master(names(simple)[-1], 'Zissou1', wes = TRUE)
    }
    else
    {
      cxcols = brewer.master(names(cx)[-1], 'Accent', wes = FALSE)
      simplecols = brewer.master(names(simple)[-1], 'Pastel1', wes = FALSE)
    }

    cxtracks = lapply(names(cx)[-1], function(x)
      anno_barplot(
        cx[[x]],
        legend = TRUE,
        axis_param = list(gp = gpar(fontsize = track.fontsize)),
        height = unit(track.height, 'cm'),
        border = FALSE,
        gp = gpar(fill = cxcols[x], col = NA)
      ))
    names(cxtracks) = names(cx)[-1]

    simpletracks = lapply(names(simple)[-1], function(x)
      anno_barplot(
        simple[[x]],
        legend = TRUE,
        axis_param = list(gp = gpar(fontsize = track.fontsize)),
        height = unit(track.height, 'cm'),
        border = FALSE,
        gp = gpar(fill = simplecols[x], col = NA)
        ))
    names(simpletracks) = names(simple)[-1]

    bottomtracks = c(bottomtracks, simpletracks, cxtracks)
  }

  ## process custom columns if any 
  if (!is.null(tumors) && length(intersect(columns, names(tumors))))
  {
    columns = intersect(columns, names(tumors))
    custom = tumors[match(ids, id), columns, with = FALSE]
    out.mat = rbind(out.mat, t(as.matrix(custom)))
    customcols = brewer.master(columns, wes = wes)
    customtracks = lapply(columns, function(x)
    {
      ## discrete data simple plot ie heatmap
      if (is.character(custom[[x]]) | is.factor(custom[[x]]) | is.logical(custom[[x]]))
      {
        if (is.logical(custom[[x]]))
          cols = c("FALSE" = 'gray', "TRUE" = 'red')
        else
          cols = brewer.master(unique(custom[[x]]), wes = wes)
        list(
          anno = anno_simple(
            as.character(custom[[x]]),
            height = unit(track.height/2, 'cm'),
            col = cols),
          legend = Legend(labels = names(cols),
                          ncol = 2, legend_gp = gpar(fill = cols, col = NA),
                          title = x)
        )
      }
      else ## numeric data barplot
        list(anno = 
               anno_barplot(
                 custom[[x]],
                 legend = TRUE,
                 axis_param = list(gp = gpar(fontsize = track.fontsize)),
                 height = unit(track.height, 'cm'),
                 border = FALSE,
                 gp = gpar(fill = customcols[x], col = NA)
               ))
    })

    customanno = lapply(customtracks, function(x) x$anno)
    names(customanno) = columns
    bottomtracks = c(bottomtracks, customanno)

    ix = lengths(customtracks)==2
    if (any(ix))
      packed_legends = c(packed_legends,
                         lapply(customtracks[ix], function(x) x$legend))
  }
  
  if (length(bottomtracks))
  {
    bottomtracks$gap = unit(track.gap, 'cm')
    bottomtracks = do.call(HeatmapAnnotation, bottomtracks)
  }

  if (length(packed_legends))
    packed_legends = do.call(packLegend, packed_legends)


  if (!is.null(split))
  {
    if (is.null(tumors))
      warning('split variable must be provided along with keyed tumors table')

    if (split %in% names(tumors))
      split = tumors[match(ids, id), ][[split]]
    else
    {
      warning('split column not found in provided tumors table')
      split = NULL
    }
  }

  gene_split = NULL
  if (!all(is.na(genes$group)))
    gene_split = genes[match(rownames(varm), genes), group]

  ## to overcome empty plot issue and also plot pct correctly
  show_pct = TRUE
  if (any(varm!='WT'))
    varm[varm == 'WT'] = ''
  else
    show_pct = FALSE ## if plot has no alterations we keep the WT so oncoPrint doesn't freak 

  if (!length(toptracks))
    toptracks = NULL

  if (!length(bottomtracks))
    bottomtracks = NULL

  op = ComplexHeatmap::oncoPrint(varm,
                      get_type = function(x) unlist(strsplit(x, ",")), ##get type = separating each cell in matrix into vector
                      alter_fun = alter_fun,
                      top_annotation = toptracks,
                      bottom_annotation = bottomtracks,
                      row_split = gene_split,
                      show_pct = show_pct, 
                      row_gap = unit(split.gap, 'cm'),
                      column_split = split,
                      column_gap = unit(split.gap, 'cm'),
                      col = varcol,
                      remove_empty_columns = FALSE,
                      remove_empty_rows = drop.genes, 
                      row_order = 1:nrow(varm),
                      column_order = 1:ncol(varm),
                      pct_gp = gpar(fontsize = rownames.fontsize),
                      row_names_gp = gpar(fontsize = rownames.fontsize),
                      column_names_gp = gpar(fontsize = colnames.fontsize),
                      show_column_names = TRUE,
                      show_heatmap_legend = TRUE
                      )


  if (ppdf)
    if (length(packed_legends))
      skitools::ppdf(draw(op, annotation_legend_list = packed_legends), height = height, width = width, ...)
    else
      skitools::ppdf(draw(op), height = height, width = width, ...)

  if (return.oncotab)
    oncotab
  else if (return.mat)
    out.mat
  else
    op
} 


#' @name annotated_bcf_to_oncotable
#' @title annotated_bcf_to_oncotable
#'
#' @description
#'
#' Collect results from annotated BCF for oncotable
#'
#' @param annotated_bcf (character) path to file
#' @param filter (character) default PASS
#' @param keepeff (character) variants to keep
#' @param id (character) ID of sample
#' @param verbose (logical) default FALSE
annotated_bcf_to_oncotable = function(annotated_bcf,
                                      filter = 'PASS',
                                      keepeff = c('trunc', 'cnadel', 'cnadup', 'complexsv', 'splice', 'inframe_indel', 'fusion', 'missense', 'promoter', 'regulatory','mir'),
                                      id = "",
                                      verbose = FALSE) {
    bcf = grok_bcf(annotated_bcf, label = id, long = TRUE, filter = filter)
    if (verbose) {
        message(length(bcf), ' variants pass filter')
    }
    
    bcf = bcf[bcf$short %in% keepeff]
    if (verbose) {
        message(length(bcf), ' variants pass keepeff')
    }
    if (length(bcf)) {
        bcf$variant.g = paste0(seqnames(bcf), ':', start(bcf), '-', end(bcf), ' ', bcf$ALT, '>', bcf$REF)
        vars = gr2dt(bcf)[, .(id = id,
                              gene, vartype,
                              variant.g,
                              variant.p, distance,
                              annotation,
                              type = short,
                              track = 'variants',
                              source = 'annotated_bcf')] %>% unique
    } else {
        vars = data.table(id = id, type = NA, source = 'annotated_bcf', track = "variants")
    }
    return(vars)
}
    

#' @name compute_rna_quantiles
#' @title compute_rna_quantiles
#'
#' @description
#'
#' Creates a data.table containing the expression level of each gene relative to a cohort.
#' If supplied, also adds role of gene (e.g. ONC/TSG)
#' 
#' @param tpm (character) path to file with RNA TPM input (can be kallisto, salmon, or prespecified format)
#' @param tpm_cohort (character) path to file with cohort expression
#' @param onc (character) character vector of oncogenes
#' @param tsg (character) character vector of TSGs
#' @param surface (character) character vector of surface genes
#' @param id (character) id of current sample
#' @param gencode.gtrack (character) path to gencode gTrack
#' @param quantile.thresh (numeric) default 0.05
#' @param verbose (logical) 
#'
#' @return data.table with columns gene, pair, value, qt, role, and direction of over/under expression
compute_rna_quantiles = function(tpm = NULL,
                                 tpm.cohort = NULL,
                                 onc = as.character(),
                                 tsg = as.character(),
                                 surface = as.character(),
                                 id = NULL,
                                 gencode.gtrack = NULL,
                                 quantile.thresh = 0.05,                                 
                                 verbose = FALSE) {

    if (is.null(id) | is.na(id)) {
        stop("id must be supplied")
    }

    if (!file.good(tpm)) {
        stop("tpm not valid")
    }

    if (!file.good(tpm.cohort)) {
        stop("tpm_cohort not valid")
    }

    if (verbose) {
        message("Reading TPM and cohort inputs")
    }

    tpm.dt = rna_reformat(tpm, pair = id, gngt.fname = gencode.gtrack)
    tpm.cohort.dt = data.table::fread(tpm.cohort, header = TRUE)

    if (!is.element("gene", colnames(tpm.cohort.dt))) {
        stop("There must be a 'gene' column in the cohort RNA expression levels containing gene symbol")
    }

    tpm.cohort.dt = tpm.cohort.dt[!duplicated(gene)]

    if (verbose) {
        message("Computing RNA expression quantiles")
    }

    melted.expr = rna_quantile(tpm.cohort.dt, id, tpm.dt)
    melted.expr = melted.expr[!is.na(value),]

    if (verbose) {
        message("Annotating with supplied TSG and oncogenes")
    }

    melted.expr[, role := ""]
    melted.expr[(gene %in% onc) & (gene %in% tsg), role := "ONC|TSG"]
    melted.expr[(gene %in% onc) & (gene %in% surface), role := "ONC|SURF"]
    melted.expr[(gene %in% onc) & !(gene %in% tsg), role := "ONC"]
    melted.expr[!(gene %in% onc) & (gene %in% surface), role := "SURF"]
    melted.expr[!(gene %in% onc) & (gene %in% tsg), role := "TSG"]

    if (verbose) {
        message("Identifying over/under expressed genes")
    }

    melted.expr[role %like% "ONC" & qt >= (1 - quantile.thresh), direction := "over"]
    melted.expr[role %like% "SURF" & qt >= (1 - quantile.thresh), direction := "over"]
    melted.expr[role %like% "TSG" & qt < quantile.thresh, direction := "under"]
	
    return(melted.expr)
}

#' @name plot_expression_histograms
#' @title plot_expression_histograms
#'
#' @description
#'
#' Produce histograms for RNA expression TPM
#' 
#' @param rna.change.fn (character) path to over/underexpressed oncogenes file
#' @param melted.expr.fn (character) path to melted expression file
#' @param pair (character)
#' @param outdir (character) where to save plots
#' @param ... additional arguments to ppng
#'
#' @return data.table with columns $gene and $expr.hist.fname
plot_expression_histograms = function(rna.change.fn = NULL,
                                      melted.expr.fn = NULL,
                                      pair = NULL,
                                      outdir = ".",
                                      ...) {

    rna.change = fread(rna.change.fn, header = TRUE)
    melted.expr = fread(melted.expr.fn, header = TRUE)

    if (!nrow(rna.change)) {
        return(data.table(gene = as.character(), expr.hist.fname = as.character()))
    }
    
    plot.dt = lapply(rna.change[, gene],
                     function (g) {
                         direction = rna.change[gene == g, direction]
                         gene.dt = melted.expr[gene == g & !is.na(value)]
                         fn = paste0(outdir, "/", g, ".", direction, ".expr.png")
                         p = ggplot(gene.dt, aes(x = value, y = gene)) +
                             geom_density_ridges2(                            
                                 bandwidth = 0.1,
                                 alpha = 0.5,
                                 scale = 0.9,
                                 rel_min_height = 0.01,
                                 color = NA,
                                 jittered_points = TRUE,
                                 position = position_points_jitter(width = 0.01, height = 0),
                                 point_shape = '|',
                                 point_size = 3,
                                 point_alpha = 0.3,
                                 point_colour = "black") +
                             geom_vline(
                                 xintercept = rna.change[gene == g, value],
                                 color = ifelse(direction=="over", "red", "blue"),
                                 lty = "dashed", lwd = 3) +
                             scale_x_continuous(trans = "log1p",
                                                breaks = c(1, 10, 100, 1000, 10000)) +
                             theme_minimal(32) +
                             theme(axis.text.x = element_text(angle = 45, vjust = 1)) +
                             labs(x = "TPM")
                         ppng(print(p), filename = fn, ...)
                         return(data.table(gene = g, expr.hist.fname = fn))
                     })
    plot.dt = rbindlist(plot.dt, fill = TRUE)
    return(plot.dt)
}

#' @name wgs_gtrack
#' @title wgs_gtrack
#'
#' @param jabba_rds (character) path to jabba
#' @param cvgt.fname (character) coverage gTrack
#' @param agt.fname (character) allele gtrack
#'
#' @return gTrack object with nice formatting
#' @export
wgs_gtrack = function(jabba_rds, cvgt.fname, agt.fname = NULL) {

    gg = gG(jabba = jabba_rds)
    cvgt = readRDS(cvgt.fname)

    y0 = 0
    gg.max.cn = max(gg$nodes$dt[!is.na(cn) & !is.infinite(cn), cn], na.rm = TRUE)
    y1 = gg.max.cn + 1
    
    ## gGraph gTrack formatting
    gg.gt = gg$gt
    gg.gt$yaxis.pretty = 4
    gg.gt$ylab = "CN"
    gg.gt$y0 = 0
    gg.gt$y1 = y1
    gg.gt$gap = 1e7 ## gap between chromosomes

    ## coverage gTrack formatting
    cvgt$yaxis.pretty = 4
    cvgt$y0 = 0
    cvgt$y1 = y1
    
    if (!file.good(agt.fname)) {
        gt = c(cvgt, gg.gt)
    } else {

        agt = readRDS(agt.fname)

        ## agtrack formatting
        agt$y0 = 0
        agt$y1 = y1
        agt$yaxis.pretty = 4

        ## concatenate with agt
        gt = c(agt, cvgt, gg.gt)
    }
    return(gt)
}


#' @name makeSummaryTables
#' @title makeSummaryTables
#' @description
#' 
#' Function for generating summary tables of interesting driver and surface genes for casereport.
#'
#' @param cnv_table file path to casereport driver copy number variants table
#' @param surface_cnv file path to casereport surface copy number variants table
#' @param fusions_table file path to casereport fusions table
#' @param expression_table file path to casereport driver over/under expression table
#' @param surface_expression file path to casereport surface over/under expression table
#' @param mutations_table file path to casereport driver mutations table
#' @param onco_table file path to casereport oncotable
#' @param onc file path to list of oncogenes to consider in the table
#' @param tsg file path to list of tsgs to consider in the table
#' @param surface file path to list of surface genes to consider in the table
#' @return summary table of driver genes.
makeSummaryTables = function(cnv_table,surface_cnv,fusions_table,expression_table,surface_expression,mutations_table,onco_table,onc,tsg,surface){
	genelist=vector()
	if(file.good(cnv_table)){
		genelist=c(genelist,fread(cnv_table)$gene_name)
	}
    if(file.good(surface_cnv)){
        genelist=c(genelist,fread(surface_cnv)$gene_name)
    }
	if(file.good(fusions_table)){
		genelist=c(genelist,fread(fusions_table)$driver.name)
	}
	if(file.good(expression_table)){
		genelist=c(genelist,fread(expression_table)$gene)
	}
    if(file.good(surface_expression)){
        genelist=c(genelist,fread(surface_expression)$gene)
    }
	if(file.good(mutations_table)){
		genelist=c(genelist,fread(mutations_table)$gene)
	}
	
    genelist=unique(genelist)
    oncotable=readRDS(onco_table)
    onc=readRDS(onc)
    tsg=readRDS(tsg)
    surface=readRDS(surface)

    summaryTable=NA
    pmkbTier=get_pmkb_tier_table(NA)
    if (length(genelist) == 0){
        return(data.table(gene = character(), role = character(),
             type = character(), tier = character(),
             source = character()))
    }
    for(i in 1:length(genelist)){
        thisGene=oncotable[oncotable$gene==genelist[i] & !is.na(oncotable$gene),]
        if(genelist[i] %in% pmkbTier$gene){
            #thisTier=min(pmkbTier[pmkbTier$gene==thisGene$gene[1],]$Tier)  
            thisTier=pmkbTier[pmkbTier$gene==thisGene$gene[1],]$tier
        }else{
            thisTier=NA
        }

        if(is.na(summaryTable)){
            summaryTable=data.table(gene=thisGene$gene[1],role=toString(unique(thisGene$role)),type=toString(unique(thisGene$type)),tier=thisTier,source=toString(unique(thisGene$source)))
        }else{
            summaryTable=rbind(summaryTable,data.table(gene=thisGene$gene[1],role=toString(unique(thisGene$role)),type=toString(unique(thisGene$type)),tier=thisTier,source=toString(unique(thisGene$source))))     
		}
    }
    
    summaryTable=summaryTable[!is.na(summaryTable$type) & summaryTable$type!=" ",]
    summaryTable=summaryTable[summaryTable$type!="NA",]

    summaryTable$type=str_replace_all(summaryTable$type,"NA, ","")
    summaryTable$role=str_replace_all(summaryTable$role,"NA, ","")
    summaryTable$type=str_replace_all(summaryTable$type,", NA","")
    summaryTable$role=str_replace_all(summaryTable$role,", NA","")
    summaryTable$type=str_replace_all(summaryTable$type,", $","")
    summaryTable$role=str_replace_all(summaryTable$role,", $","")
       summaryTable$role=str_replace_all(summaryTable$role,", ,",",")

    summaryTable$withHetdel=ifelse(grepl("hetdel",summaryTable$type),"True","False")    

    summaryTable=summaryTable[order(summaryTable$tier,summaryTable$withHetdel),]
    summaryTable$withHetdel=NULL

    summaryTable=summaryTable[!(summaryTable$type=="del" & summaryTable$role=="ONC"),]
    summaryTable=summaryTable[!(summaryTable$type=="amp" & summaryTable$role=="TSG"),]
    summaryTable=summaryTable[!(grepl("underexpression",summaryTable$type,fixed=TRUE) & summaryTable$role=="ONC"),]
    summaryTable=summaryTable[!(grepl("overexpression",summaryTable$type,fixed=TRUE) & summaryTable$role=="TSG"),]

    summaryTable$type=str_replace_all(summaryTable$type,"del"," loss")
	summaryTable$type=str_replace_all(summaryTable$type,"het loss","hetdel")
    driverTable=summaryTable[(gene %in% onc | gene %in% tsg) & !(gene %in% surface)]
    surfaceTable=summaryTable[!(gene %in% onc | gene %in% tsg) & (gene %in% surface)]
    driverTable$gene=paste0('<a href=https://www.oncokb.org/gene', driverTable$gene, ' target=_blank rel=noopener noreferrer >', driverTable$gene, '</a>')
    surfaceTable$gene=paste0('<a href=https://www.oncokb.org/gene', surfaceTable$gene, ' target=_blank rel=noopener noreferrer >', surfaceTable$gene, '</a>')
    
	return(list(driverTable,surfaceTable))
}


#' @name summarize_cases
#' @title summarize_cases
#' @description
#' 
#' Very ad-hoc function for mskilab use only!
#'
#' takes a case report flow module and produces a data.table with links to the reports
#' this is assuming that the reports are somewhere under /gpfs/commons/projects/imielinski_web
#'
#' @param jb Flow Job object or a data.frame/data.table containing the contents of outputs(jb), where jb is the case report Flow Job object.
#' @param output_file output TXT in which the tabular data will be saved. If not file is provided then the data is saved to a temporary file.
#' @param libdir path to the casereport repository clone
#' @param html_dir path to the directory in which to put the html version of the table
#' @param metadata data.frame or data.table with additional metadata for the samples in jb. The metadata table must contain a header. The first column of the metadata should contain sample names that match the sample names in the Job object and should only contain unique values, otherwise this table is ignored and no metadata is added.
#' @return data.table
#' @export
summarize_cases = function(jb, output_file = NULL, libdir = '~/git/casereport', html_dir = NULL, metadata = NULL){
    if (is.character(jb) && grepl('rds$', jb)){
        jb = readRDS(jb)
    }
    if (!inherits(jb, 'Job') && !inherits(jb, 'data.frame')){
        stop('jb must be a data.frame or an object of class Flow::Job, but you provided: ', class(jb))
    }
    if (!is.null(output_file) && !is.character(output_file)){
        stop('Invalid value for output_file. output_file must be of class character, but you provided: ', class(output_file))
    }
    if (inherits(jb, 'data.frame')){
        dt = as.data.table(jb)
        if (is.null(key(dt))){
            setkeyv(dt, names(dt)[1])
        }
    } else {
        dt = outputs(jb)
    }
    if (!('wgs_casereport' %in% names(dt))){
        stop('Invalid input. The input job information must contain column "wgs_casereport"')
    }
    summary.fn = paste0(dirname(dt[, wgs_casereport]),'/summary.rds')
    dt[, link := paste0('<a href=', wgs_casereport,
                        '>report</a>')]
    dt[, link:=gsub('/gpfs/commons/projects/imielinski_web/','//mskiweb.nygenome.org/', link)]
    k = dt %>% key
    iids = dt[,get(k)]
    names(summary.fn) = iids
    summ = lapply(iids, function(ix){
        fn = summary.fn[ix]
        if (file.good(fn)){
            summary.dt = readRDS(fn) %>% as.data.table
            summary.dt[, id := ix]
            return(summary.dt)
        }
        sdt = data.table(id = ix)
        return(sdt)
    })
    summ = rbindlist(summ, fill = T)
    dt = merge.data.table(dt[,.(id = get(k), link)], summ, by = 'id')
    if (!is.null(output_file) && is.character(output_file) && dir.exists(dirname(output_file))){
        message('Writing table to: ', output_file)
        fwrite(dt, output_file)
    }
    if (!is.null(metadata) && is.data.frame(metadata) && all(!duplicated(as.data.table(metadata[, 1])))){
        metadata = as.data.table(metadata)
        setnames(metadata, names(metadata)[1], 'id')
        dt = merge.data.table(metadata, dt, all.y = T, by = 'id')
    }

    if (!is.null(html_dir) && dir.exists(html_dir)){
        html_path = normalizePath(paste0(html_dir, "/case.reports.html"))
        message("Generating html output to: ", html_path)
        if (!file.good(output_file)){
            output_file = tempfile()
            message('Writing table to: ', output_file)
            fwrite(dt, output_file)
        }
        rmarkdown::render(
            input = normalizePath(system.file('extdata', 'case_report_module/wgs.report.table.rmd', package = 'casereport')),
            output_format = "html_document",
            output_file = html_path,
            knit_root_dir = normalizePath(html_dir),
            ## params = report.config,
            params = list(summary_table = normalizePath(output_file)),
            quiet = FALSE)
    }
    return(dt)
}


#' @name set_param
#' @title set_param
#' @description
#' 
#' check if a field exists in a list and if not then set it
#'
#' @param l input list
#' @param name the name of the field
#' @param value the value to assign if the field is NULL (by default NA_character_)
#' @return list
set_param = function(l = list(), name = '', value = NA_character_){
    if (!(name %in% names(l))){
        l[name] = value
    }
    return(l)
}

#' @name try2
#' @title wrapper around tryCatch - robust to parallel:: functions
#'
#' A slightly more robust version of try that works within the parallel:: set of functions
#' that pre-deploy a cluster.
#'
#' @param expr an R expression to try.
#' @param finally expression to be evaluated before returning or exiting.
try2 = function(expr, ..., finally) {
    tryCatch(expr,
             error = function(e) {
                 msg = structure(paste(conditionMessage(e), conditionCall(e), sep = "\n"), class = "err")
                 cat("Error: ", msg, "\n\n")
                 return(msg)
             },
             finally = finally,
             ... = ...)
}

#' @name rand.string
#' @title make a random string
#'
#' @return random string
#' @author Someone from Stackoverflow
#' @export rand.string
rand.string <- function(n=1, length=12)
{
    randomString <- c(1:n)                  # initialize vector
    for (i in 1:n)
    {
        randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                        length, replace=TRUE),
                                 collapse="")
    }
    return(randomString)
}
