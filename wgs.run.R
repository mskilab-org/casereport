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
        make_option(c("--enhancer"), type = "character", help = "annotation of putative active enhancers in the tissue type, used for proximity analysis"),
        make_option(c("--deconstruct_sigs"), type = "character", default = NA_character_, help = "deconstruct_sigs module output, RDS"),
        make_option(c("--deconstruct_variants"), type = "character", default = NA_character_, help = "deconstruct_sigs module variant output, TXT"),
        make_option(c("--sigs_cohort"), type = "character", default = NA_character_, help = "variant count for each signature in a cohort"),
        make_option(c("--tpm"), type = "character", default = NA_character_, help = "Textual file containing the TPM values of genes in this sample (raw kallisto output acceptable)"),
        make_option(c("--tpm_cohort"), type = "character", default = NA_character_, help = "Textual file containing the TPM values of genes in a reference cohort"),
        make_option(c("--hrd_results"), type = "character", default = NA_character_, help = "The comprehensive HRDetect module results"),
        make_option(c("--snpeff_snv"), type = "character", default = NA_character_, help = "snpeff snv results"),
        make_option(c("--snpeff_indel"), type = "character", default = NA_character_, help = "snpeff indel results"),
        make_option(c("--gencode"), type = "character", default = "~/DB/GENCODE/hg19/gencode.v19.annotation.gtf", help = "GENCODE gene models in GTF/GFF3 formats"),
        make_option(c("--genes"), type = "character", default = 'http://mskilab.com/fishHook/hg19/gencode.v19.genes.gtf', help = "GENCODE gene models collapsed so that each gene is represented by a single range. This is simply a collapsed version of --gencode."),
        make_option(c("--drivers"), type = "character", default = NA_character_, help = "path to file with gene symbols (see /data/cgc.tsv for example)"),
        make_option(c("--chrom_sizes"), type = "character", default = "~/DB/UCSC/hg19.broad.chrom.sizes", help = "chrom.sizes file of the reference genome"),
        make_option(c("--knit_only"), type = "logical", default = FALSE, action = "store_true", help = "if true, skip module and just knit"),
        make_option(c("--amp_thresh"), type = "numeric", default = 4,
                    help = "Threshold over ploidy to call amplification"),
        make_option(c("--del_thresh"), type = "numeric", default = 0.5,
                    help = "Threshold over ploidy to call deletion"),
        make_option(c("--server"), type = "character", default = "https://mskilab.com/gGraph/", help = "URL of the gGnome.js browser"),
        make_option(c("--tumor_type"), type = "character", default = "", help = "tumor type"),
        make_option(c("--ref"), type = "character", default = "hg19", help = "one of 'hg19', 'hg38'"),
        make_option(c("--snpeff_config"), type = "character", default = "~/modules/SnpEff/snpEff.config", help = "snpeff.config file path"),
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
        library(forcats)
        library(stringr)
        library(gGnome)
        library(gTrack)
        library(gUtils)
        library(tidyr)
        library(dplyr)
        library(ggforce)
        library(ggridges)
        library(ggrepel)
        library(httr)
        library(jsonlite)
        library(knitr)
        library(rmarkdown)
        library(ComplexHeatmap)
        library(deconstructSigs)
        library(DT)
        message("Loading critical dependencies from KevUtils")
        source(paste0(opt$libdir, "/utils.R"))
        source(paste0(opt$libdir, "/config.R"))
        source(file.path(opt$libdir, "sv.gallery.R"))
        source(file.path(opt$libdir, "fusion.gallery.R"))
    })
})

if (!opt$knit_only){
    message("Preparing data and plots")

    message("Returning Purity, Ploidy, and run 'events' if not already provided")
    jabba = readRDS(opt$jabba_rds)

    if (file.exists(paste0(opt$outdir, "/complex.rds"))){
        gg = readRDS(paste0(opt$outdir, "/complex.rds"))
    } else if (file.exists(opt$complex) & file.size(opt$complex)>0){
        file.copy(opt$complex, paste0(opt$outdir, "/complex.rds"))
        gg = readRDS(opt$complex)
    } else {
        if (file.exists(opt$complex) & file.size(opt$complex)>0){
            file.copy(opt$complex, paste0(opt$outdir, "/complex.rds"))
            gg = readRDS(opt$complex)
        } else {
            gg = events(gG(jabba = jabba))
            saveRDS(gg, paste0(opt$outdir, "/complex.rds"))
        }
    }

    message("Checking for RNA expression input")
    if (file.good(opt$tpm_cohort) & file.good(opt$tpm)) {

        ## save reformatted TPM data
        tpm.fn = file.path(opt$outdir, "tpm.txt")

        tpm.dt = kallisto.preprocess(opt$tpm,
                                     pair = opt$pair,
                                     gngt.fname = file.path(opt$libdir, "data", "gt.ge.hg19.rds"))

        ## check for and process kallisto input
        fwrite(tpm.dt, tpm.fn)
        opt$tpm = tpm.fn

        message("Found and preprocessed TPM input")
    }

    message('Calling CNVs for oncogenes and tumor suppressor genes')
    ## # get the ncn data from jabba
    genes_cn.fn = paste0(opt$outdir, '/genes_cn.rds')
    oncogenes.fn = file.path(opt$libdir, "data", "onc.rds")
    tsg.fn = file.path(opt$libdir, "data", "tsg.rds")
    if (check_file(genes_cn.fn, overwrite = opt$overwrite)){
        genes_cn_annotated = readRDS(genes_cn.fn)
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
            ## # TODO: perhaps we should allow to directly provide a nseg input
            message('JaBbA karyograph was not found at the expected location (', kag_rds, ') so we will use CN = 2 for the normal copy number of all chromosomes.')
        }

        #' zchoo Monday, May 03, 2021 02:33:55 PM
        #' add ploidy to avoid bug
        #' simplify seqnames to avoid empty data tables
        genes_cn = get_gene_copy_numbers(gg, gene_ranges = opt$genes, nseg = nseg, ploidy = kag$ploidy, simplify_seqnames = TRUE, complex.fname = opt$complex)
        genes_cn_annotated = get_gene_ampdel_annotations(genes_cn, amp.thresh = opt$amp_thresh,
                                       del.thresh = opt$del_thresh)



        #' zchoo Wednesday, May 12, 2021 04:17:06 PM
        #' integrate RNA expression data with genes data
        if (file.good(opt$tpm_cohort) & file.good(opt$tpm)) {
            message("Computing quantiles")
            melted.expr = rna.quantile(opt$tpm_cohort,
                                       opt$pair,
                                       opt$tpm)
            message("Done computing quantiles")
            genes_cn_annotated = merge.data.table(genes_cn_annotated,
                                                  melted.expr[pair == opt$pair,
                                                              .(gene_name = gene, expr.quantile = qt, expr.value = value)],
                                                  by = "gene_name",
                                                  all.x = TRUE)
            message("Done merging")
        } else {
            genes_cn_annotated[, ":="(expr.quantile = NA, expr.value = NA)]
        }

        ## add over/under expression annotations
        qt.thres = 0.05 ## expose this parameter later
        genes_cn_annotated[expr.quantile < qt.thres, expr := "under"]
        genes_cn_annotated[expr.quantile > (1 - qt.thres), expr := "over"]

        ## save annotated genes
        saveRDS(genes_cn_annotated, genes_cn.fn)
    }

    driver.genes.cnv.fn = paste0(opt$outdir, '/driver.genes.cnv.txt')
    driver.genes.expr.fn = paste0(opt$outdir, '/driver.genes.expr.txt')
    if (!check_file(driver.genes.cnv.fn, overwrite = opt$overwrite)){
        if (genes_cn_annotated[, .N] > 0){
            onc = readRDS(oncogenes.fn)
            tsg = readRDS(tsg.fn)
            #' zchoo Tuesday, May 04, 2021 10:41:25 PM
            ## save just expression change and just CNA separately
            driver.genes_cn = genes_cn_annotated[(cnv == "amp" & gene_name %in% onc) |
                                                 (cnv == "del" & gene_name %in% tsg)]
            driver.genes_expr = genes_cn_annotated[(expr == "over" & gene_name %in% onc) |
                                                   (expr == "under" & gene_name %in% tsg)]
            #' zchoo Wednesday, May 12, 2021 05:00:24 PM
            ## subset these to make less overwhelming...
            ## fields = c("gene_name", "cnv", "min_cn", "max_cn", "min_normalized_cn", "max_normalized_cn", "number_of_cn_segments", "ncn", "seqnames", "start", "end", "width", "gene_id", "gene_type", "source",  "level", "hgnc_id", "havana_gene", "ev.id", "ev.type")
            fields = c("gene_name", "cnv", "expr", "min_cn", "max_cn", "min_normalized_cn", "max_normalized_cn", "expr.value", "expr.quantile", "seqnames", "start", "end", "width", "ev.id", "ev.type")
            fields = intersect(fields, names(driver.genes_cn))
            fwrite(driver.genes_cn[, ..fields], driver.genes.cnv.fn)
            fwrite(driver.genes_expr[, ..fields], driver.genes.expr.fn)
        }
    }

    ## prepare driver gallery
    opt$complex = paste0(opt$outdir, "/complex.rds")

    message("Prepare coverage data")
    cvgt_fn = paste0(opt$outdir, "/coverage.gtrack.rds")
    if (check_file(cvgt_fn, overwrite = opt$overwrite)){
        cvgt = readRDS(cvgt_fn)
    } else {
        ## pull coverage file from jabba_rds
        cov.file = readRDS(file.path(dirname(opt$jabba_rds), "cmd.args.rds"))$coverage
        cvgt = covcbs(cov.file, purity = jabba$purity, ploidy = jabba$ploidy, rebin = 5e3,
                      ylab = "CN", y.cap = FALSE, xaxis.chronly = TRUE)
        saveRDS(cvgt, cvgt_fn)
    }

    agt_fn = paste0(opt$outdir, "/agtrack.rds")
    if (check_file(agt_fn, overwrite = opt$overwrite)) {
        agt = readRDS(agt_fn)
    } else {
        message("Checking for hets")
        if (is.null(opt$het_pileups_wgs) || !file.exists(opt$het_pileups_wgs)) {
            message("no het pileups provided, skipping")
            agt = NULL
        } else {
            agt = grab.agtrack(opt$het_pileups_wgs,
                               purity = jabba$purity,
                               ploidy = jabba$ploidy)
            saveRDS(agt, agt_fn)
        }
    }

    ## xtYao legacy code
    ## if (!file.exists(paste0(opt$outdir, "/hets.gtrack.rds"))){
    ##     hgt = covcbs(opt$het_pileups_wgs,
    ##                  purity = jabba$purity, ploidy = jabba$ploidy, rebin = 5e3)
    ##     saveRDS(cvgt, paste0(opt$outdir, "/coverage.gtrack.rds"))
    ## } else {
    ##     cvgt = readRDS(paste0(opt$outdir, "/coverage.gtrack.rds"))
    ## }
    ## message("Generate full genome gTrack plot")
    ## gts = c(cvgt, gg$gtrack(name = opt$pair, height = 30))
    ## png(filename = paste0(opt$outdir, "/genome.wide.gtrack.png"), width = 2700, height = 900)
    ## plot(gts, hg, gap = 1e7, y0 = 0, cex.label = 2, yaxis.cex = 0.5)
    ## dev.off()

    ## TODO: make this generic
    ## cgc = fread(paste0(opt$libdir, "/data/cgc.tsv"))
    onc = readRDS(paste0(opt$libdir, "/data/onc.rds"))
    tsg = readRDS(paste0(opt$libdir, "/data/tsg.rds"))

    wgs.gtrack.fname = file.path(opt$outdir, "wgs.gtrack.png")
    wgs.circos.fname = file.path(opt$outdir, "wgs.circos.png")
    if (opt$overwrite | !file.exists(wgs.gtrack.fname)) {
        message("Generating whole-genome gTrack plots")
        ## synchronize y0 and y1 of cvgt and agt
        y0 = 0
        gg.max.cn = max(gg$nodes$dt[!is.na(cn) & !is.infinite(cn), cn], na.rm = TRUE)
        y1 = round(gg.max.cn/10) * 10
        
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
        
        if (is.null(agt)) {
            gt = c(cvgt, gg.gt)
        } else {

            ## agtrack formatting
            agt$y0 = 0
            agt$y1 = y1
            agt$yaxis.pretty = 4

            ## concatenate with agt
            gt = c(agt, cvgt, gg.gt)
        }
        ppng(plot(gt, c(as.character(1:22), "X", "Y")),
             filename  = wgs.gtrack.fname,
             height = 1500,
             width = 5000)
    } else {
        message("Whole genome gTracks already exist")
    }

    if (opt$overwrite | !file.exists(wgs.circos.fname)) {
        message("Generating whole-genome circos plot")
        ppng(wgs.circos(junctions = gg$junctions[type == "ALT"],
                        cov = cvgt@data[[1]],
                        field = "cn",
                        link.h.ratio = 0.1,
                        cex.points = 0.1,
                        cytoband.path = file.path(opt$libdir, "data", "hg19.cytoband.txt")),
             filename = wgs.circos.fname,
             height = 1000,
             width = 1000)
    } else {
        message("Whole genome circos plot already exists")
    }

    ## ##################
    ## Driver gene mutations
    ## ##################
    driver.mutations.fname = file.path(opt$outdir, "driver.mutations.txt")
    snv.bcf.fname = file.path(opt$outdir, "snpeff", "snv", "annotated.bcf")
    indel.bcf.fname = file.path(opt$outdir, "snpeff", "indel", "annotated.bcf")
    if (opt$overwrite | !file.exists(driver.mutations.fname)) {

        ## ################
        ## read cgc gene list
        ## ################
        
        cgc.fname = ifelse(is.null(opt$drivers) || is.na(opt$drivers),
                           file.path(opt$libdir, "data", "cancerGeneList.tsv"),
                           opt$drivers)

        ## ################
        ## run SnpEff
        ## ################
        
        if ( (!is.null(opt$snpeff_snv)) && file.exists(opt$snpeff_snv)) {

            message("Running SNV SnpEff")
            snpeff.libdir = normalizePath(file.path(opt$libdir, "SnpEff_module"))
            snpeff.ref = opt$ref
            snpeff.vcf = normalizePath(opt$snpeff_snv)
            snpeff.outdir = normalizePath(file.path(opt$outdir, "snpeff", "snv"))
            snpeff.config = normalizePath(opt$snpeff_config)

            snpeff.cm = paste("sh", file.path(snpeff.libdir, "run.sh"),
                              snpeff.libdir,
                              snpeff.ref,
                              snpeff.vcf,
                              snpeff.outdir,
                              snpeff.config)

            system(snpeff.cm)

            ## grab annotations
            snv.dt = filter.snpeff(vcf = snv.bcf.fname, ##opt$snpeff_snv,
                                   gngt.fname = file.path(opt$libdir, "data", "gt.ge.hg19.rds"),
                                   cgc.fname = cgc.fname,
                                   ref.name = opt$ref,
                                   verbose = TRUE)            
        } else {
            message("SNV vcf does not exist")
            snv.dt = NULL
        }
        if ( (!is.null(opt$snpeff_indel)) && file.exists(opt$snpeff_indel)) {

            message("Running indel SnpEff")
            snpeff.libdir = file.path(opt$libdir, "SnpEff_module")
            snpeff.ref = opt$ref
            snpeff.vcf = opt$snpeff_indel
            snpeff.outdir = file.path(opt$outdir, "snpeff", "indel")
            snpeff.config = normalizePath(opt$snpeff_config)

            snpeff.cm = paste("sh", file.path(snpeff.libdir, "run.sh"),
                              snpeff.libdir,
                              snpeff.ref,
                              snpeff.vcf,
                              snpeff.outdir,
                              snpeff.config)

            system(snpeff.cm)

            indel.dt = filter.snpeff(vcf = indel.bcf.fname, ## opt$snpeff_indel,
                                     gngt.fname = file.path(opt$libdir, "data", "gt.ge.hg19.rds"),
                                     cgc.fname = cgc.fname,
                                     ref.name = opt$ref,
                                     verbose = TRUE)

        } else {
            message("indel snv does not exist")
            indel.dt = NULL
        }

        if (!is.null(snv.dt) & !is.null(indel.dt)) {
            driver.mutations.dt = rbind(snv.dt, indel.dt, use.names =  TRUE, fill = TRUE)
        } else if (!is.null(snv.dt)) {
            driver.mutations.dt = snv.dt
        } else {
            driver.mutations.dt = data.table()
        }
        fwrite(driver.mutations.dt, driver.mutations.fname)
    } else {
        message("Driver mutations file already exists")
    }

    ## ##################
    ## Fusions
    ## ##################
    fusions.driver.fname = file.path(opt$outdir, "fusions.driver.txt")
    fusions.other.fname = file.path(opt$outdir, "fusions.other.txt")
    if (opt$overwrite | !file.exists(fusions.driver.fname) | !file.exists(fusions.other.fname)) {
        ## if opt$fusions not available, run it
        if (!file.exists(opt$fusions) | file.size(opt$fusions)==0){
            if (file.exists(paste0(opt$outdir, "/fusions.rds"))){
                opt$fusions = paste0(opt$outdir, "/fusions.rds")
            } else {
                if (!exists("gff")){
                    gff = skidb::read_gencode(fn = opt$gencode)
                }
                fu = fusions(gg, gff)
                saveRDS(fu, paste0(opt$outdir, "/fusions.rds"))
                opt$fusions = paste0(opt$outdir, "/fusions.rds")
                saveRDS(opt, "cmd.args.rds")
            }
        }
        message("Preparing fusion genes report")
        ## grab name of driver genes file
        cgc.fname = ifelse(is.null(opt$drivers) || is.na(opt$drivers), file.path(opt$libdir, "data", "cgc.tsv"), opt$drivers)
        fusions.slickr.dt = fusion.wrapper(fusions.fname = opt$fusions,
                                           complex.fname = opt$complex,
                                           cvgt.fname = cvgt_fn,
                                           agt.fname = agt_fn,
                                           cgc.fname = cgc.fname,
                                           gngt.fname = file.path(opt$libdir, "data", "gt.ge.hg19.rds"),
                                           pad = 0.5,
                                           height = 2000,
                                           width = 1000,
                                           outdir = opt$outdir)

        ## make sure data table is not empty
        if (nrow(fusions.slickr.dt) == 0) {

            ## write empty data table, rmd file will deal with this.
            fwrite(fusions.slickr.dt, fusions.driver.fname)
            fwrite(fusions.slickr.dt, fusions.other.fname)
            
        } else {

            ## save data table for drivers and non-drivers separately
            fwrite(fusions.slickr.dt[driver == TRUE,], fusions.driver.fname)
            fwrite(fusions.slickr.dt[driver == FALSE,], fusions.other.fname)

        }
        
    } else {
        message("Fusion files already exist")
    }
    
    ## ##################
    ## SV gallery code
    ## ##################
    ## generate gTrack with just cgc genes
    cgc.fname = ifelse(is.null(opt$drivers) || is.na(opt$drivers),
                       file.path(opt$libdir, "data", "cgc.tsv"),
                       opt$drivers)
    cgc.gtrack.fname = cgc.gtrack(cgc.fname = cgc.fname,
                                  gencode.fname = opt$gencode,
                                  outdir = opt$outdir)

    if (opt$overwrite | !file.exists(file.path(opt$outdir, "sv.gallery.txt"))) {
        message("Preparing SV gallery")

        
        sv.slickr.dt = gallery.wrapper(complex.fname = opt$complex,
                                       background.fname = file.path(opt$libdir, "data", "sv.burden.txt"),
                                       cvgt.fname = cvgt_fn,
                                       gngt.fname = file.path(opt$libdir, "data", "gt.ge.hg19.rds"),
                                       cgcgt.fname = cgc.gtrack.fname,
                                       agt.fname = agt_fn,
                                       server = opt$server,
                                       pair = opt$pair,
                                       pad = 0.5,
                                       height = 1200, ## png image height
                                       width = 1000, ## png image width
                                       outdir = opt$outdir)

        fwrite(sv.slickr.dt, file.path(opt$outdir, "sv.gallery.txt"))
    } else {
        message("SV gallery files already exist")
    }

    ## #################
    ## CN gallery
    ## #################
    cn.gallery.fn = file.path(opt$outdir, "cn.gallery.txt")
    if (!check_file(cn.gallery.fn, overwrite = opt$overwrite)) {
        message("preparing CN gallery")
        ## grab ploidy
        pl = readRDS(opt$jabba_rds)$ploidy
        cn.slickr.dt = cn.plot(drivers.fname = driver.genes.cnv.fn,
                               opt$complex,
                               cvgt.fname = cvgt_fn,
                               gngt.fname = file.path(opt$libdir, "data", "gt.ge.hg19.rds"),
                               cgcgt.fname = cgc.gtrack.fname,
                               agt.fname = agt_fn,
                               server = opt$server,
                               pair = opt$pair,
                               amp.thresh = opt$amp_thresh,
                               ploidy = pl,
                               pad = 0.5,
                               height = 1600,
                               width = 1000,
                               outdir = opt$outdir)
        fwrite(cn.slickr.dt, cn.gallery.fn)
    } else {
        message("CN gallery files already exist")
    }

    ## ##################
    ## RNA expression level over a cohort
    ## ##################
    if (file.good(opt$tpm_cohort)){
        tpm_cohort = fread(opt$tpm_cohort, header = TRUE)
        if (is.element(opt$pair, colnames(tpm_cohort))){
            message("Found this sample in the cohort expression matrix")
            if (file.good(opt$tpm)){
                tpm = fread(opt$tpm, header = TRUE)
                
                message("Found this sample's input expression matrix, overwriting...")
                tpm_cohort[[opt$pair]] = NULL
                tpm_cohort = data.table::merge.data.table(
                    tpm_cohort, tpm, by = "gene", all.x = TRUE)
            }
        }
        ## limit to annotated ONC/TSG
        tpm_cohort = tpm_cohort[gene %in% c(onc, tsg)]
        tpm_cohort[, role := case_when(gene %in% onc ~ "ONC",
                                       gene %in% tsg ~ "TSG",
                                       TRUE ~ NA_character_)]
        ## melt to analyze
        mexp = data.table::melt(tpm_cohort, id.vars = c("gene", "role"),
                                variable.name = "pair")
        mexp[, ":="(qt = rank(as.double(.SD$value))/.N), by = gene]
        ## TODO: make the quantile threshold adjustable
        cool.exp = mexp[pair==opt$pair][(role=="TSG" & qt<0.05) | (role=="ONC" & qt>0.95)]
        ## cool.exp[order(qt)]

        cool.exp.fn = file.path(opt$outdir, "cool.expr.rds")
        if (nrow(cool.exp)>0){
            cool.exp[, direction := ifelse(qt>0.95, "over", "under")]
            cool.exp[, gf := paste0(opt$outdir, "/", gene, ".", direction, ".expr.png")]
            setkey(cool.exp, "gene")
            for (g in cool.exp$gene){
                if (!file.good(cool.exp[g, gf])){
                    png(filename = cool.exp[g, gf], width = 1600, height = 900)
                    d = mexp[gene==g][!is.na(value)]
                    message(g)
                    ## ppng({
                    p = ggplot(d, aes(x = value, y = gene)) +
                        geom_density_ridges2() +
                        geom_vline(
                            xintercept = cool.exp[g, value],
                            color = ifelse(cool.exp[g, direction]=="over", "red", "blue"),
                            lty = "dashed", lwd = 3) +
                        scale_x_continuous(trans = "log1p",
                                           breaks = c(1, 10, 100, 1000, 10000)) +
                        ## scale_y_continuous(expand = c(0, 0)) +
                        theme_minimal(32) +
                        theme(axis.text.x = element_text(angle = 45, vjust = 1)) +
                        labs(x = "TPM")
                    print(p)
                    ## }, width = 1600, height = 900)
                    dev.off()
                }
            }
            saveRDS(cool.exp, cool.exp.fn)
        }

        ## expression change gallery
        expr.gallery.fn = file.path(opt$outdir, "expr.gallery.txt")
        if (!check_file(expr.gallery.fn, overwrite = opt$overwrite) & file.exists(cool.exp.fn)) {
            message("preparing expression gallery")
            #grab ploidy
            pl = readRDS(opt$jabba_rds)$ploidy
            expr.slickr.dt = cn.plot(drivers.fname = cool.exp.fn,
                                     opt$complex,
                                     cvgt.fname = cvgt_fn,
                                     gngt.fname = file.path(opt$libdir, "data", "gt.ge.hg19.rds"),
                                     cgcgt.fname = cgc.gtrack.fname,
                                     agt.fname = agt_fn,
                                     server = opt$server,
                                     pair = opt$pair,
                                     amp.thresh = opt$amp_thresh,
                                     ploidy = pl,
                                     pad = 0.5,
                                     height = 1600,
                                     width = 1000,
                                     outdir = opt$outdir)
            fwrite(expr.slickr.dt, expr.gallery.fn)
        } else {
            message("expression gallery files already exist")
        }
        waterfall.fn = file.path(opt$outdir, "waterfall.png")
        if (!check_file(waterfall.fn, overwrite = opt$overwrite) & file.exists(cool.exp.fn)) {
            message("generating waterfall plot")
            gns = readRDS(cool.exp.fn)$gene ## genes with changes in expression
            rna.waterfall.plot(tpm.cohort = opt$tpm_cohort,
                               pair = opt$pair,
                               tpm.pair = opt$tpm,
                               out.fn = waterfall.fn,
                               genes.to.label = gns)
        } else {
            message("Expression gallery files already exist")
        }
    }

    ## ##################
    ## SNV signatures
    ## ##################
    if (file.good(opt$deconstruct_sigs)){
        ## signatures
        sig = readRDS(opt$deconstruct_sigs)
        sigd = as.data.table(sig$weights) %>% data.table::melt(variable.name = "Signature", value.name = "Proportion") %>% data.table
        sigd[, Signature := factor(Signature, levels = rev(levels(Signature)))]
        sigd[, pair := opt$pair]

        ## counts of variants
        var = fread(opt$deconstruct_variants)[nchar(max.post)>0]
        sct = var[, .(sig_count = .N), by = .(Signature = max.post)]
        sct[, pair := opt$pair]

        if (!file.exists(paste0(opt$outdir, "/deconstruct_sigs.png")) | opt$overwrite){
            ## ppng({
            png(filename = paste0(opt$outdir, "/deconstruct_sigs.png"), width = 1000, height = 1000)
            deconstructSigs::plotSignatures(sig)
            dev.off()
            ## makePie(sig)
            ## }, "deconstruct_sigs.png", width = 1600, height = 1200)
        }

        if (!file.exists(paste0(opt$outdir, "/sig.composition.png")) | opt$overwrite){

            sig.fn = file.path(opt$libdir, "data", "all.signatures.txt")
            background.type = "Cell cohort"
            
            if (file.good(opt$sigs_cohort)){
                sig.fn = opt$sigs_cohort
                background.type = "supplied" ## what is the background dist, used for plot title
            } else {
                if (!is.null(opt$tumor_type) & !is.na(opt$tumor_type)) {

                    ## tumor type specific signature
                    tumor.sig.fn = file.path(opt$libdir, "data", paste0(opt$tumor_type, ".signatures.txt"))

                    ## check if file exists for specific tumor type
                    if (file.exists(tumor.sig.fn)) {
                        message("Using supplied background signature burden for tumor type: ", opt$tumor_type)
                        sig.fn = tumor.sig.fn
                        background.type = opt$tumor_type
                    } else {
                        message("Using background signature burden for whole Cell cohort")
                    }
                } else {
                    message("Using background signature burden for whole Cell cohort")
                }
            }

            
            allsig = data.table::fread(sig.fn)
            allsig = data.table::melt(allsig, id = "pair", variable.name = "Signature", value.name = "sig_count")

            ## recast as character to avoid releveling for now
            allsig[, Signature := as.character(Signature)]
            sct[, Signature := as.character(Signature)]

            if (opt$pair %in% allsig[, pair]){
                ## remove existing record from the cohort first
                allsig = allsig[pair!=opt$pair]
            }

            ## sct[, Signature := factor(Signature, levels = levels(allsig$Signature))][!is.na(Signature)]
            allsig = rbind(allsig, sct)
            
            ## calculate percentile
            allsig[, perc := rank(sig_count)/.N, by = Signature]

            ## highlight the tracks where the signature is non-zero in this sample
            ## allsig[, highlight := as.character(Signature) %in% as.character(sct[, Signature])]
            allsig = allsig[!is.na(Signature) & (Signature %in% sct$Signature)]

            ## order levels...
            keep.slevels = allsig[pair == opt$pair, Signature] %>% as.character
            ## allsig[, Signature := factor(Signature, levels = new.slevels)]

            allsig = allsig[!is.na(Signature) & Signature %in% keep.slevels,][, Signature := gsub("Signature.", "", Signature)]
            ## just level one time
            new.slevels = allsig[pair == opt$pair,][order(sig_count), Signature]
            allsig[, Signature := factor(Signature, levels = new.slevels)]
            
            sigbar = ggplot(allsig,aes(y = Signature, x = sig_count, fill = Signature)) +
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
                geom_segment(data = allsig[pair==opt$pair & sig_count>0],
                             aes(x = sig_count, xend = sig_count,
                                 y = as.numeric(Signature), yend = as.numeric(Signature) + 0.9),
                             color = "red", size = 1, alpha = 0.5) +
                geom_label(data = allsig[pair==opt$pair],
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
                labs(title = paste0("Signatures vs. ", background.type, " background"), x = "Burden") +
                theme_minimal() +
                theme(legend.position = "none",
                      title = element_text(size = 20, family = "sans"),
                      axis.title = element_text(size = 20, family = "sans"),
                      axis.text.x = element_text(size = 15, family = "sans"),
                      axis.text.y = element_text(size = 20, family = "sans"))

            ## theme(
            ##     text = element_text(size = 32),
            ##     axis.text.x
            ##     legend.position = "none"
            ##     ## axis.text.x = element_text(angle = 45, hjust = 1)
            ## )

            ppng(print(sigbar), filename = paste0(opt$outdir, "/sig.composition.png"), height = 800, width = 800)
            
            
        } else {
            png(filename = "sig.composition.png", width = 800, height= 1200)
            ## individual sample compositions
            sct[, Signature := forcats::fct_reorder(Signature, sig_count)]
            sigbar = ggplot(
                sct[order(sig_count)],
                aes(x = Signature, y = sig_count)) +
                geom_bar(stat = "identity") +
                theme_minimal() +
                theme(text = element_text(size = 32),
                      axis.text.x = element_text(angle = 45, hjust = 1)) +
                coord_flip()
            print(sigbar)
            dev.off()
        }
    }

    ## ##################
    ## HRDetect results
    ## ##################
    if (!file.good(paste0(opt$outdir, "/hrdetect.rds")) | opt$overwrite){
        if (file.good(opt$hrd_results)){
            hrd.res = readRDS(opt$hrd_results)
            ## melt the input data matrix
            hrd.dat = as.data.table(hrd.res$data_matrix) %>% data.table::melt()
            hrd.dat[, pair := opt$pair]
            ## melt the output results
            hrd.out = as.data.table(hrd.res$hrdetect_output) %>% data.table::melt()
            hrd.out[, pair := opt$pair]
            saveRDS(hrd.out, paste0(opt$outdir, "/hrdetect.rds"))

            hrd = merge(hrd.out, hrd.dat[, .(variable, data = value)], by = "variable")

            ## original training data
            hrd_cohort = fread(paste0(opt$libdir, "/data/hrdetect.og.txt"))
            hrd_cohort = data.table::melt(hrd_cohort, id.vars = c("pair", "is.hrd"))
            ## if (opt$pair %in% hrd_cohort$pair){
            hrd_cohort = rbindlist(list(hrd_cohort[pair!=opt$pair], hrd.dat), fill = TRUE)
            ## }

            hrd_cohort[, range(value), by = .(variable, is.hrd)]

            hrd.yes = hrd.out[variable=="Probability", value>0.7]

            ## these are the dimensions that need to be logged
            ldat = hrd_cohort[variable != "del.mh.prop"]
            ldat[, variable := factor(variable, levels = setdiff(levels(variable), "del.mh.prop"))]
            ## plot the distribution of raw count
            hrd.plot = ggplot(ldat, aes(x = value, y = variable)) +
                geom_density_ridges(aes(point_colour = is.hrd,
                                        fill = is.hrd),
                                    bandwidth = 0.1,
                                    alpha = 0.5,
                                    scale = 0.9,
                                    rel_min_height = 0.01,
                                    color = NA,
                                    jittered_points = TRUE,
                                    position = position_points_jitter(width = 0.01, height = 0),
                                    point_shape = '|',
                                    point_size = 3,
                                    point_alpha = 0.3) +
                scale_discrete_manual("point_colour", values = binary.cols) +
                scale_fill_manual(values = binary.cols) +
                geom_segment(data = ldat[pair==opt$pair],
                             aes(x = value, xend = value,
                                 y = as.numeric(variable), yend = as.numeric(variable) + 0.9),
                             color = ifelse(hrd.yes, "#D7261E", "#1f78b4"),
                             alpha = 1, lwd = 3) +
                scale_x_continuous(trans = "log1p",
                                   breaks = c(0, 1, 10, 100, 1000, 10000, 100000),
                                   labels = c(0, 1, 10,
                                              expression(10^2),
                                              expression(10^3),
                                              expression(10^4),
                                              expression(10^5)),
                                   limits = c(0, 100000),
                                   expand = c(0, 0)) +
                labs(title = paste0("HRDetect features compared to Davies et al. 2017"), x = "Burden") +
                theme_minimal() +
                theme(legend.position = "none",
                      title = element_text(size = 20, family = "sans"),
                      axis.title = element_text(size = 20, family = "sans"),
                      axis.text.x = element_text(size = 15, family = "sans"),
                      axis.text.y = element_text(size = 20, family = "sans"))

            ## these are the dimensions that need to be logged
            rdat = hrd_cohort[variable == "del.mh.prop"]
            rdat[, variable := factor(variable, levels = "del.mh.prop")]
            ## plot the distribution of raw count
            hrd.plot.2 = ggplot(rdat, aes(x = value, y = variable)) +
                geom_density_ridges(aes(point_colour = is.hrd,
                                        fill = is.hrd),
                                    bandwidth = 0.1,
                                    alpha = 0.5,
                                    scale = 0.9,
                                    rel_min_height = 0.01,
                                    color = NA,
                                    jittered_points = TRUE,
                                    position = position_points_jitter(width = 0.01, height = 0),
                                    point_shape = '|',
                                    point_size = 3,
                                    point_alpha = 0.3) +
                scale_discrete_manual("point_colour", values = binary.cols, labels = NULL) +
                scale_fill_manual(values = binary.cols,
                                  labels = c("Non HR deficient", "HR deficient")) +
                ## scale_y_discrete(labels = hrdetect.dims) +
                labs(x = "Proportion") +
                guides(point_colour = FALSE) +
                geom_segment(data = rdat[pair==opt$pair],
                             aes(x = value, xend = value,
                                 y = as.numeric(variable), yend = as.numeric(variable) + 3),
                             color = ifelse(hrd.yes, "#D7261E", "#1f78b4"),
                             alpha = 1, lwd = 3) +
                ## scale_x_continuous(trans = "log1p",
                ##                    breaks = c(0, 1, 10, 100, 1000, 10000, 100000),
                ##                    labels = c(0, 1, 10,
                ##                               expression(10^2),
                ##                               expression(10^3),
                ##                               expression(10^4),
                ##                               expression(10^5)),
                ##                    limits = c(0, 100000)) +
                ## labs(title = paste0("HRDetect features compared to Davies et al. 2017"), x = "Burden") +
                theme_minimal() +
                theme(legend.position = "bottom",
                      legend.text = element_text(size=20),
                      title = element_text(size = 20, family = "sans"),
                      axis.title = element_text(size = 20, family = "sans"),
                      axis.text.x = element_text(size = 15, family = "sans"),
                      axis.text.y = element_text(size = 20, family = "sans"))

            ## finally draw the output probability in a waterfall plot
            ## hrd.plot.3 = ggplot(hrd_cohort, aes(x = rank)) %>%
            

            ## draw the plots
            png(paste0(opt$outdir, "/hrdetect.log.dat.png"), width = 800, height = 800)
            print(hrd.plot)
            dev.off()

            png(paste0(opt$outdir, "/hrdetect.prop.dat.png"), width = 800, height = 200)
            print(hrd.plot.2)
            dev.off()
        }
    }

    ## ##################
    ## Enhancer/gene proximity
    ## ##################
    if (file.good(opt$proximity) && file.good(cool.exp.fn <- paste0(opt$outdir, "/cool.expr.rds"))){
        prox = readRDS(opt$proximity)
        pdt = prox$dt
        pdt = pdt[reldist<0.25 & refdist>5e6]
        cool.exp = readRDS(cool.exp.fn)

        if (any(cool.exp[direction=="over", gene] %in% pdt$gene_name)){
            ## filter by overexpressed genes
            pdt = merge.data.table(pdt, cool.exp[direction=="over"], by.x = "gene_name", by.y = "gene", all.x = TRUE)
            pgs = pdt[(direction=="over"), unique(gene_name)]
            gt.cgc = readRDS(cgc.gtrack.fname)
            gt.cgc$name = "CGC"
            enh = readRDS(opt$enhancer)

            ## ## create a dir for the plots
            ## prox.plot.dir = paste0(opt$outdir, "/proximity.")
            
            for (g in pgs){
                this.png = paste0(opt$outdir, "/", g, ".proximity.png")
                if (!file.exists(this.png) | opt$overwrite){
                    w = prox[pdt[gene_name==g, walk.id]]
                    this.enh = copy(enh)
                    this.enh$col = binary.cols[as.character(seq_along(enh) %in% w$dt$qid)]
                    gt.enh = gTrack(this.enh, name = "enhancer", height = 5)
                    pgt = c(cvgt, gg$gtrack(height = 30), w$gtrack(name = "shortest walk"), gt.cgc, gt.enh)
                    png(this.png, height = 1200, width = 1600)
                    plot(pgt, w$footprint + 1e6, legend.params = list(plot = FALSE), y0 = 0)
                    dev.off()
                }
            }
            
        }
    }

    ## ##################
    ## Purity ploidy sliders
    ## ##################
    ## get the purity ploidy
    if (file.good(opt$jabba_rds)){
        kagf = gsub("jabba.simple", "karyograph", opt$jabba_rds)
        if (file.good(kagf)){

            ## start gathering the LL matrix
            kag = readRDS(kagf)
            pu0 = gg$meta$purity
            pl0 = gg$meta$ploidy
            seg = kag$segstats %Q% (strand=="+" & !bad) %>% gr.stripstrand
            ppg = ppgrid(segstats = seg, return.nll = TRUE)
            nll = ppg$NLL

            ## first plot, histogram with grids
            tmp = seg  ## only plot seg that we haven't fixed SD for and that have normal cn 1, to minimize confusion
            dupval = sort(table(tmp$mean), decreasing = TRUE)[1]
            if (!is.na(dupval))
                if (dupval>5)
                    tmp = tmp[-which(as.character(tmp$mean) == names(dupval))]

            ## sampling random loci to plot not segments
            ## segsamp = pmin(sample(tmp$mean, 1e6, replace = T, prob = width(tmp)), xlim[2])
            tmp = tmp %Q% (!is.na(mean))
            segsamp = sample(tmp[, "mean"], 1e6, replace = T, prob = width(tmp)) %>% gr2dt
            segsamp[, mean := pmax(xlim[1], pmin(mean, xlim[2]))]
            segsamp = segsamp[!is.na(mean)]            

            mu = tmp$mean
            w = width(tmp)
            Sw = sum(as.numeric(width(tmp)))
            sd = tmp$sd
            m0 = sum(as.numeric(mu*w))/Sw
            ncn = rep(2, length(mu))
            ploidy_normal = 2

            alpha = pu0
            tau = pl0
            beta = (tau + ploidy_normal * gamma / 2 ) / m0
            gamma = 2/alpha - 2
            grids = 1/beta*(0:100) + gamma/beta
            grids = grids[which(grids < segsamp[, max(mean)])]
            ## make static image
            h = ggplot(segsamp, aes(x = mean)) +
                geom_histogram(bins = 1000) +
                geom_vline(xintercept = grids, color = alpha("red", 0.5), lty = 2) +
                scale_x_continuous(limits = c(0, segsamp[, quantile(mean, 0.999)])) +
                labs(title = sprintf('Purity: %.2f Ploidy: %.2f Beta: %.2f Gamma: %.2f', alpha, tau, beta, gamma),
                     x = "Intensity") +
                theme_minimal(12)

            png(paste0(opt$outdir, "/ppfit.hist.png"), width = 1200, height = 800)
            print(h)
            dev.off()
            ## hp = ggplotly(h)



            ## line <- list(
            ##     type = "line",
            ##     line = list(color = alpha("#e6194b", 0.5), lty = 2),
            ##     xref = "x",
            ##     yref = "y"
            ## )
            ## lines <- list()
            ## for (i in seq_along(grids)) {
            ##     line[["x0"]] <- grids[i]
            ##     line[["x1"]] <- grids[i]
            ##     line[c("y0", "y1")] <- c(0, 1e5)
            ##     lines <- c(lines, list(line))
            ## }
            ## ##
            ## hp3 = layout(hp2,
            ##              title = sprintf('Purity: %.2f Ploidy: %.2f Beta: %.2f Gamma: %.2f', alpha, tau, beta, gamma),
            ##              shapes = lines)
            ## wij(ggplotly(hp3))

            ## PP grid
            ## purity.min = 0.01
            ## purity.max = 1.0
            ## ploidy.step = 0.4
            ## purity.step = 0.2
            ## ploidy.min = 1.8 # ploidy bounds (can be generous)
            ## ploidy.max = 4.2
            ##  purity.guesses = seq(0, 1, purity.step)
            ## pu.guesses = seq(pmax(0, purity.min), pmin(1.00, purity.max), purity.step)
            ## ploidy.guesses = seq(pmin(0.5, ploidy.min), pmax(10, ploidy.max), ploidy.step)
            ## pl.guesses = seq(pmax(0.5, ploidy.min), pmax(0.5, ploidy.max), ploidy.step)
            ## ppgr = as.data.table(expand.grid(list(purity = pu.guesses, ploidy = pl.guesses)))
            ## ppgr[, ":="(visible = FALSE, name = sprintf("Purity: %.2f Ploidy: %.2f", purity, ploidy))]
            ## ppgr[(abs(purity-pu0)<0.001 & abs(ploidy-pl0)<0.001), visible := TRUE]
            

            ## hp2 = plot_ly(x = segsamp$mean,
            ##               type = "histogram")
            ## ## create sliders for purity and ploidy values
            ## pu.scroll = list()
            ## for (i in seq_along(pu.guesses)){
            ##     pu.scroll[[i]] <- list(
            ##         method = "relayout",
            ##         ## args = list('visible', ppgr[, ifelse(abs(purity-pu.guesses[[i]])<0.001, TRUE, FALSE)])
            ##         args = list('visible', ifelse(seq_along(pu.guesses)==i, TRUE, FALSE))
            ##     )
            ## }
            ## pl.scroll = list()
            ## for (j in seq_along(pl.guesses)){
            ##     pl.scroll[[j]] <- list(
            ##         method = "relayout",
            ##         ## args = list('visible', ppgr[, ifelse(abs(ploidy-pl.guesses[[j]])<0.001, TRUE, FALSE)])
            ##         args = list('visible', ifelse(seq_along(pl.guesses)==j, TRUE, FALSE))
            ##     )
            ## }
            
            ## for (i in seq_along(pu.guesses)){
            ##     for (j in seq_along(pl.guesses)){
            ##         alpha = pu.guesses[[i]]
            ##         tau = pl.guesses[[j]]
            ##         beta = (tau + ploidy_normal * gamma / 2 ) / m0
            ##         gamma = 2/alpha - 2
            ##         grids = 1/beta*(0:100) + gamma/beta
            ##         grids = grids[which(grids < segsamp[, max(mean)])]
            ##         hp2 <- add_segments(
            ##             hp2, x = grids, y = 0, xend = grids, yend = 1e5,
            ##             ## visible = pu.scroll[[i]]$args[[2]] & pl.scroll[[i]]$args[[2]],
            ##             visible = FALSE,
            ##             name = sprintf('Purity: %.2f Ploidy: %.2f Beta: %.2f Gamma: %.2f', alpha, tau, beta, gamma),
            ##             hoverinfo = "name",
            ##             showlegend = FALSE)
            ##     }
            ## }
            
            ## ## for (i in seq_len(nrow(ppgr))){
            ## ## for (i in 1:101){
            ##     ## alpha = ppgr[i, purity]
            ##     ## tau = ppgr[i, ploidy]
            ##     ## add value to the scroll
            ##     ## pu.sc = list(args = list('visible', rep(FALSE, nrow(ppgr))),
            ##     ##              method = "restyle")
            ##     ## pu.sc$args[[2]][i] = TRUE
            ##     ## pp.scroll[[i]] = pu.sc
            ##     ## pl.sc = list(args = list('visible', rep(FALSE, nrow(ppgr))),
            ##     ##              method = "restyle")
            ##     ## pl.sc$args[[2]][i] = TRUE
            ##     ## pl.scroll[[i]] = pl.sc
            ## ##     message(i)
            ## ## }

            

            ## ## layout final plot with scroll bars
            ## hp3 <- hp2 %>%
            ##     layout(
            ##         sliders = list(
            ##             list(
            ##                 active = which(abs(purity.guesses-pu0)<0.001),
            ##                 currentvalue = list(prefix = "Purity = "),
            ##                 pad = list(t = 20),
            ##                 steps = pu.scroll),
            ##             list(
            ##                 active = which(abs(ploidy.guesses-pl0)<0.001),
            ##                 currentvalue = list(prefix = "Ploidy = "),
            ##                 pad = list(t = 100),
            ##                 steps = pl.scroll)
            ##         ))
            ## wij(hp3)
            
        } else {
            message("Did not find karyograph in the same directory as JaBbA output.")
        }
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
