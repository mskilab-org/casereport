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
