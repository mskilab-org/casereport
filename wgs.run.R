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
        make_option(c("--deconstruct_sigs"), type = "character", default = NA_character_, help = "deconstruct_sigs module output, RDS"),
        make_option(c("--deconstruct_variants"), type = "character", default = NA_character_, help = "deconstruct_sigs module variant output, TXT"),
        make_option(c("--sigs_cohort"), type = "character", default = NA_character_, help = "variant count for each signature in a cohort"),
        make_option(c("--tpm"), type = "character", default = NA_character_, help = "Textual file containing the TPM values of genes in this sample"),
        make_option(c("--tpm_cohort"), type = "character", default = NA_character_, help = "Textual file containing the TPM values of genes in a reference cohort"),
        make_option(c("--hrd_output"), type = "character", default = NA_character_, help = "The one line CSV of HRDetect final call"),
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

        saveRDS(genes_cn_annotated, genes_cn.fn)
    }

    driver.genes.cnv.fn = paste0(opt$outdir, '/driver.genes.cnv.txt')
    if (!check_file(driver.genes.cnv.fn, overwrite = opt$overwrite)){
        if (genes_cn_annotated[, .N] > 0){
            onc = readRDS(oncogenes.fn)
            tsg = readRDS(tsg.fn)
            #' zchoo Tuesday, May 04, 2021 10:41:25 PM
            ## subsetted so that there are just dels in tsgs and amps in oncogenes
            ## driver.genes_cn = genes_cn_annotated[gene_name %in% c(onc, tsg)]
            driver.genes_cn = genes_cn_annotated[(cnv == "amp" & gene_name %in% onc) |
                                                 (cnv == "del" & gene_name %in% tsg)]
            fields = c("gene_name", "cnv", "min_cn", "max_cn", "min_normalized_cn", "max_normalized_cn", "number_of_cn_segments", "ncn", "seqnames", "start", "end", "width", "gene_id", "gene_type", "source",  "level", "hgnc_id", "havana_gene", "ev.id", "ev.type")
            fields = intersect(fields, names(driver.genes_cn))
            fwrite(driver.genes_cn[, ..fields], driver.genes.cnv.fn)
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
                                           cvgt = cvgt_fn,
                                           agt.fname = agt_fn,
                                           cgc.fname = cgc.fname,
                                           gngt = file.path(opt$libdir, "data", "gt.ge.hg19.rds"),
                                           pad = 0.5,
                                           height = 2000,
                                           width = 1000,
                                           outdir = opt$outdir)

        ## save data table for drivers and non-drivers separately
        fwrite(fusions.slickr.dt[driver == TRUE,], fusions.driver.fname)
        fwrite(fusions.slickr.dt[driver == FALSE,], fusions.other.fname)
        
    } else {
        message("Fusion files already exist")
    }
    
    ## ##################
    ## SV gallery code
    ## ##################
    if (opt$overwrite | !file.exists(file.path(opt$outdir, "sv.gallery.txt"))) {
        message("Preparing SV gallery")

        ## generate gTrack with just cgc genes
        cgc.fname = ifelse(is.null(opt$drivers) || is.na(opt$drivers),
                           file.path(opt$libdir, "data", "cgc.tsv"),
                           opt$drivers)
        cgc.gtrack.fname = cgc.gtrack(cgc.fname = cgc.fname,
                                      gencode.fname = opt$gencode,
                                      outdir = opt$outdir)
        
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
        cn.slickr.dt = cn.plot(drivers.fname = driver.genes.cnv.fn,
                               opt$complex,
                               cvgt.fname = cvgt_fn,
                               gngt.fname = file.path(opt$libdir, "data", "gt.ge.hg19.rds"),
                               cgcgt.fname = cgc.gtrack.fname,
                               agt.fname = agt_fn,
                               server = opt$server,
                               pair = opt$pair,
                               pad = 0.5,
                               height = 1200,
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
            saveRDS(cool.exp, paste0(opt$outdir, "/cool.expr.rds"))
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
    if (file.good(opt$hrd_output)){
        hrd = fread(opt$hrd_output)
        hrd[, pair := opt$pair]
        hrd = data.table::melt(hrd, id.var = "pair")

        ## ## If there's no cohort to compare to
        ## hrd.res = ggplot(hrd[!variable %in% c("intercept", "Probability")],
        ##                  aes(x = variable, y = value)) +
        ##     geom_bar(stat = "identity")
        saveRDS(hrd, paste0(opt$outdir, "/hrdetect.rds"))
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
