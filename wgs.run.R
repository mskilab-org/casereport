library(optparse)

if (!exists("opt")){
    option_list = list(
        make_option(c("--libdir"), type = "character", help = "dir that contains this file and other source codes"),
        make_option(c("--outdir"), type = "character", default = './', help = "Directory to dump output into"),
        make_option(c("--pair"), type = "character", help = "ID of this sample"),
        make_option(c("--jabba_rds"), type = "character", help = "jabba output, 'jabba.simple.rds'"),
        make_option(c("--cbs_cov_rds"), type = "character", default = NA_character_, help = "cbs_cov output"),
        make_option(c("--het_pileups_wgs"), type = "character", default = NA_character_, help = "cbs_cov output"),
        make_option(c("--complex"), type = "character", help = "complex event caller, RDS"),
        make_option(c("--fusions"), type = "character", help = "fusions module output, RDS"),
        make_option(c("--proximity"), type = "character", help = "proximity module output, RDS"),
        make_option(c("--deconstruct_sigs"), type = "character", help = "deconstruct_sigs module output, RDS"),    
        make_option(c("--gencode"), type = "character", default = "~/DB/GENCODE/hg19/gencode.v19.annotation.gtf", help = "GENCODE gene models in GTF/GFF3 formats"),
        make_option(c("--chrom_sizes"), type = "character", default = "~/DB/UCSC/hg19.broad.chrom.sizes", help = "chrom.sizes file of the reference genome"),
        make_option(c("--knit_only"), type = "logical", default = FALSE, action = "store_true", help = "if true, skip module and just knit"),
        make_option(c("--amp_thresh"), type = "numeric", default = 4,
                    help = "Threshold over ploidy to call amplification"),
        make_option(c("--del_thresh"), type = "numeric", default = 0.5,
                    help = "Threshold over ploidy to call deletion"),
        make_option(c("--server"), type = "character", default = "https://mskilab.com/gGraph/", help = "URL of the gGnome.js browser"),
        make_option(c("--tumor_type"), type = "character", default = "", help = "tumor type"),
        make_option(c("--overwrite"), type = "logical", default = FALSE, action = "store_true", help = "overwrite existing data in the output dir")
    )
    parseobj = OptionParser(option_list = option_list)
    opt = parse_args(parseobj)
    opt$outdir = normalizePath(opt$outdir)
    opt$libdir = normalizePath(opt$libdir)
    saveRDS(opt, paste0(opt$outdir, '/cmd.args.rds'))
}

message("Loading Libraries -- Please wait...")
suppressMessages(expr = {
    suppressPackageStartupMessages(expr = {
        library(stringr)
        library(gGnome)
        library(gTrack)
        library(gUtils)
        library(tidyr)
        library(dplyr)
        library(ggforce)
        library(ggridges)
        library(httr)
        library(jsonlite)
        library(knitr)
        library(rmarkdown)
        library(wesanderson)
        message("Loading critical dependencies from KevUtils")
        source(paste0(opt$libdir, "/utils.R"))
        source(paste0(opt$libdir, "/config.R"))
        source(file.path(opt$libdir, "sv.gallery.R"))
        source(file.path(opt$libdir, "fusion.gallery.R"))
    })
})

if (!opt$knit_only) {
    message("Preparing data and plots")

    message("Returning Purity, Ploidy, and run 'events' if not already provided")
    jabba = readRDS(opt$jabba_rds)
    gg_fn = paste0(opt$outdir, "/complex.rds")
    if (check_file(gg_fn, overwrite = opt$overwrite)){
        gg = readRDS(gg_fn)
    } else {
        if (file.exists(opt$complex) & file.size(opt$complex)>0){
            file.copy(opt$complex, paste0(opt$outdir, "/complex.rds"))
            gg = readRDS(opt$complex)
        } else {
            gg = events(gG(jabba = jabba))
            saveRDS(gg, paste0(opt$outdir, "/complex.rds"))
        }
    }

    message('Calling CNVs for oncogenes and tumor suppressor genes')
    # get the ncn data from jabba
    genes_cn.fn = paste0(opt$outdir, '/genes_cn.rds')
    oncogenes.fn = file.path(opt$libdir, "data", "onc.rds")
    tsg.fn = file.path(opt$libdir, "data", "tsg.rds")
    if (check_file(genes_cn.fn, overwrite = opt$overwrite)){
        genes_cn = readRDS(genes_cn.fn)
    } else {
        kag_rds = gsub("jabba.simple.rds", "karyograph.rds", opt$jabba_rds)
        nseg = NULL
        if (file.exists(kag_rds) & file.size(kag_rds) > 0){
            message('Loading JaBbA karyograph from ', kag_rds)
            kag = readRDS(kag_rds)
            if ('ncn' %in% names(mcols(kag$segstats))){
                nseg = kag$segstats[,c('ncn')]
            } else {
                message('JaBbA karyograph does not contain the "ncn" field so CN = 2 will be assumed for the normal copy number of all seqnames.')
            }
        } else {
            # TODO: perhaps we should allow to directly provide a nseg input
            message('JaBbA karyograph was not found at the expected location (', kag_rds, ') so we will use CN = 2 for the normal copy number of all chromosomes.')
        }

        genes_cn = get_gene_copy_numbers(gg, pge = pge, nseg = nseg)
        genes_cn_annotated = get_gene_ampdel_annotations(genes_cn, amp.thresh = opt$amp_thresh,
                                       del.thresh = opt$del_thresh)

        saveRDS(genes_cn_annotated, genes_cn.fn)
    }

    driver.genes.cnv.fn = paste0(opt$outdir, '/driver.genes.cnv.txt')
    if (!check_file(driver.genes.cnv.fn, overwrite = opt$overwrite)){
        if (genes_cn[, .N] > 0){
            onc = readRDS(oncogenes.fn)
            tsg = readRDS(tsg.fn)
            driver.genes_cn = genes_cn_annotated[gene_name %in% c(onc, tsg)]
            fields = c("gene_name", "cnv", "min_cn", "max_cn", "min_normalized_cn", "max_normalized_cn", "number_of_cn_segments", "ncn", "seqnames", "start", "end", "width", "gene_id", "gene_type", "source",  "level", "hgnc_id", "havana_gene")
            fields = intersect(fields, names(driver.genes_cn))
            fwrite(driver.genes_cn[, ..fields], driver.genes.cnv.fn)
        }
    }

    message("Prepare coverage data")
    cvgt_fn = paste0(opt$outdir, "/coverage.gtrack.rds")
    if (check_file(cvgt_fn, overwrite = opt$overwrite)){
        cvgt = readRDS(cvgt_fn)
    } else {
        ## pull coverage file from jabba_rds
        cov.file = readRDS(file.path(dirname(opt$jabba_rds), "cmd.args.rds"))$coverage
        cvgt = covcbs(cov.file, purity = jabba$purity, ploidy = jabba$ploidy, rebin = 5e3)
        saveRDS(cvgt, cvgt_fn)
    }

    wgs.gtrack.fname = file.path(opt$outdir, "wgs.gtrack.png")
    wgs.circos.fname = file.path(opt$outdir, "wgs.circos.png")
    if (opt$overwrite | !file.exists(wgs.gtrack.fname)) {
        message("Generating whole-genome gTrack plots")
        ppng(plot(c(cvgt, gg$gt), c(as.character(1:22), "X", "Y")),
             filename  = wgs.gtrack.fname,
             height = 1000,
             width = 3000)
    } else {
        message("Whole genome gTracks already exist")
    }

    if (opt$overwrite | !file.exists(wgs.circos.fname)) {
        message("Generating whole-genome circos plot")
        ppng(circos(junctions = gg$junctions[type == "ALT"],
                    cov = cvgt@data[[1]],
                    field = "cn",
                    link.h.ratio = 0.1,
                    cex.points = 0.1),
             filename = wgs.circos.fname,
             height = 1000,
             width = 1000)
    } else {
        message("Whole genome circos plot already exists")
    }
    
    if (opt$overwrite | !file.exists(file.path(opt$outdir, "fusions.driver.txt"))) {
        message("Preparing fusion genes report")
        fusions.slickr.dt = fusion.wrapper(fusions.fname = opt$fusions,
                                           complex.fname = opt$complex,
                                           cvgt.fname = file.path(opt$outdir, "coverage.gtrack.rds"),
                                           cgc.fname = file.path(opt$libdir, "data", "cgc.tsv"),
                                           file.path(opt$libdir, "data", "gt.ge.hg19.rds"),
                                           pad = 1e5,
                                           height = 1000,
                                           width = 1000,
                                           outdir = opt$outdir)

        ## save data table for drivers and non-drivers separately
        fwrite(fusions.slickr.dt[driver == TRUE,], file.path(opt$outdir, "fusions.driver.txt"))
        fwrite(fusions.slickr.dt[driver == FALSE,], file.path(opt$outdir, "fusions.other.txt"))
    } else {
        message("Fusion files already exist")
    }
                                        

    if (opt$overwrite | !file.exists(file.path(opt$outdir, "sv.gallery.txt"))) {
        message("Preparing SV gallery")
        sv.slickr.dt = gallery.wrapper(complex.fname = opt$complex,
                                       background.fname = file.path(opt$libdir, "data", "sv.burden.txt"),
                                       cvgt.fname = file.path(opt$outdir, "coverage.gtrack.rds"),
                                       server = opt$server,
                                       pair = opt$pair,
                                       pad = 5e5,
                                       height = 1000, ## png image height
                                       width = 1000, ## png image width
                                       outdir = opt$outdir)

        fwrite(sv.slickr.dt, file.path(opt$outdir, "sv.gallery.txt"))
    } else {
        message("SV gallery files already exist")
    }
}



message("Start knitting")
rmarkdown::render(
    input = paste0(opt$libdir, "/wgs.report.rmd"),
    output_format = "html_document",
    output_file = paste0(opt$outdir, "/", opt$pair,".wgs.report.html"),
    knit_root_dir = opt$outdir,
    params = list(set_title = paste0(opt$pair),
                  pair = opt$pair,
                  jabba_rds = normalizePath(opt$jabba_rds),
                  outdir = normalizePath(opt$outdir),
                  tumor_type = opt$tumor_type,
                  server = opt$server),
    quiet = FALSE)
