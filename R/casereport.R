run.wgs.report = function(opt){
    message("Loading Libraries -- Please wait...")
    suppressMessages(expr = {
        suppressPackageStartupMessages(expr = {
            library(forcats)
            library(stringr)
            library(gGnome)
            library(gTrack)
            library(gUtils)
            library(tidyr)
            library(tidyverse)
            library(ggforce)
            library(ggridges)
            library(ggrepel)
            library(ggExtra)
            library(httr)
            library(jsonlite)
            library(knitr)
            library(rmarkdown)
            library(ComplexHeatmap)
            library(deconstructSigs)
            library(DT)
            library(VariantAnnotation)
            library(immunedeconv)
            library(readr)
            library(data.table)
            message("Loading critical dependencies from KevUtils")
            source(paste0(opt$libdir, "/utils.R"))
            source(paste0(opt$libdir, "/config.R"))
            source(paste0(opt$libdir, "/pmkb-utils.R"))
            source(file.path(opt$libdir, "sv.gallery.R"))
            source(file.path(opt$libdir, "fusion.gallery.R"))
            data.table::setDTthreads(1)
        })
    })

    ##################
    ## Initialize a report config
    ##################

    if (file.good(paste0(opt$outdir, "/", "report.config.rds"))) {
        report.config = readRDS(paste0(opt$outdir, "/", "report.config.rds"))
    } else {

        ## copy report.config from opt
        report.config = opt

        ## extra params for knitting
        report.config$set_title = opt$pair

        ## normalize paths for knitting
        report.config$jabba_rds = normalizePath(report.config$jabba_rds)
        report.config$outdir = normalizePath(report.config$outdir)

        ## add gTrack file names to report config
        report.config$coverage_gtrack = paste0(report.config$outdir, "/coverage.gtrack.rds")
        report.config$allele_gtrack = paste0(report.config$outdir, "/agtrack.rds")
        report.config$cgc_gtrack = paste0(report.config$outdir, "/", "cgc.gtrack.rds")
        report.config$gencode_gtrack = system.file("extdata", "gt.ge.hg19.rds", package = "casereport") ## this is provided

        if (check_file(opt$drivers))
            report.config$drivers = opt$drivers
        else
            report.config$drivers = NA_character_
        

        ## add CGC genes file
        report.config$cgc = cgc.fname = system.file("extdata", "cgc.tsv", package = "casereport")
        
        
        ## add oncogenes, etc. to report config
        report.config$onc = system.file("extdata", "onc.rds", package = "casereport")
        report.config$tsg = system.file("extdata", "tsg.rds", package = "casereport")
        report.config$surface = system.file("extdata", "surface.rds", package = "casereport")

        ## purity/ploidy plots
        report.config$cn_plot = paste0(report.config$outdir, "/", "cn.pp.png")
        report.config$allele_plot = paste0(report.config$outdir, "/", "allele.scatter.png")

        ## SCNA
        report.config$gene_cn = paste0(report.config$outdir, "/", "genes_cn.rds")
        report.config$driver_scna = paste0(report.config$outdir, '/driver.genes.cnv.txt')
        report.config$scna_gtracks = paste0(report.config$outdir, "/", "cn.gallery.txt")

        ## SNVS
        report.config$driver_mutations = paste0(report.config$outdir, "/", "driver.mutations.txt")

        ## SV
        report.config$sv_gtracks = paste0(report.config$outdir, "/", "sv.gallery.txt")

        ## whole genome vis
        report.config$wgs_gtrack_plot = file.path(report.config$outdir, "wgs.gtrack.png")
        report.config$wgs_circos_plot = file.path(report.config$outdir, "wgs.circos.png")

        ## fusions
        report.config$driver_fusions = file.path(report.config$outdir, "fusions.driver.txt")
        report.config$other_fusions = file.path(report.config$outdir, "fusions.other.txt")

        ## RNA expression analyses
        report.config$tpm_quantiles = paste0(report.config$outdir, "/", "tmp.quantiles.txt")
        report.config$rna_change = paste0(report.config$outdir, "/", "rna.change.txt")
        report.config$rna_change_all = paste0(report.config$outdir, "/", "rna.change.all.txt")
        report.config$expression_histograms = paste0(report.config$outdir, "/", "expr.histograms.txt")
        report.config$expression_gtracks = paste0(report.config$outdir, "/", "expr.gallery.txt")
        report.config$waterfall_plot = paste0(report.config$outdir, "/", "waterfall.png")
        report.config$rna_change_with_cn = paste0(report.config$outdir, "/", "driver.genes.expr.txt")

        ## HRD
        report.config$ot = paste0(report.config$outdir, "/ot.rds")
        report.config$ot_log = paste0(report.config$outdir, "/onenesstwoness.log.dat.png")
        report.config$ot_prob = paste0(report.config$outdir, "/onenesstwoness.prop.dat.png")
        report.config$oneness = paste0(report.config$outdir, "/Oneness.png")
        report.config$twoness = paste0(report.config$outdir, "/Twoness.png")

        ## deconstructSigs
        report.config$sig_composition = paste0(report.config$outdir, "/deconstruct_sigs.png")
        report.config$sig_histogram = paste0(report.config$outdir, "/sig.composition.png")
        
        ## Deconvolution
        report.config$deconv = paste0(report.config$outdir, "/", "deconv_results.txt")
        
        ## summary
        report.config$summary_stats = paste0(report.config$outdir, "/summary.rds")
        report.config$oncotable = paste0(report.config$outdir, "/oncotable.rds")
        report.config$summaryTable = paste0(report.config$outdir, "/summaryTable.txt")

        saveRDS(report.config, paste0(report.config$outdir, "/", "report.config.rds"))
    }

    ##################
    ## Start producing some analyses
    ##################

    if (!opt$knit_only) {
        message("Grabbing purity and ploidy")
        jabba = readRDS(opt$jabba_rds)

        if (is.null(report.config$purity)) {
            report.config$purity = jabba$purity
            saveRDS(report.config, paste0(opt$outdir, "/", "report.config.rds"))
        }

        if (is.null(report.config$ploidy)) {
            report.config$ploidy = jabba$ploidy
            saveRDS(report.config, paste0(opt$outdir, "/", "report.config.rds"))
        }
        
        
        message("Preparing gTracks")
        message("Preparing coverage gTrack")
        if (check_file(report.config$coverage_gtrack, overwrite = opt$overwrite, verbose = opt$verbose)){
            cvgt = readRDS(report.config$coverage_gtrack)
        } else {
            cov.file = opt$cbs_cov_rds
            cvgt = covcbs(cov.file, purity = jabba$purity, ploidy = jabba$ploidy, rebin = 5e3,
                          ylab = "CN", y.cap = FALSE, xaxis.chronly = TRUE)
            saveRDS(cvgt, report.config$coverage_gtrack)
        }

        if (check_file(report.config$allele_gtrack, overwrite = opt$overwrite)) {
            agt = readRDS(report.config$allele_gtrack)
        } else {
            message("Checking for hets")
            if (file.good(opt$het_pileups_wgs)) {
                agt = grab.agtrack(opt$het_pileups_wgs,
                                   purity = jabba$purity,
                                   ploidy = jabba$ploidy)
                saveRDS(agt, report.config$allele_gtrack)
            } else {
                message("no het pileups provided, skipping")
                agt = NULL            
            }
        }

        
        if (check_file(report.config$cgc_gtrack, opt$overwrite, opt$verbose)) {
            message("CGC gTrack exists")
        } else {
            cgc.gt = create_cgc_gtrack(cgc.fname = report.config$cgc,
                                       gencode.fname = opt$gencode)
            saveRDS(cgc.gt, report.config$cgc_gtrack)
        } 


        ####################
        ## Events
        ####################
        ## if (file.good(report.config$complex)) {
        ##     gg = readRDS(report.config$complex)
        ## } else if (file.good(opt$complex)) {
        ##     file.copy(opt$complex, report.config$complex)
        ##     gg = gGnome::refresh(readRDS(opt$complex))
        ## } else {
        ##     gg = gGnome::events(gG(jabba = jabba))
        ##     saveRDS(gg, report.config$complex)
        ## }
        ## gg = read

        ## gg$set(purity = jabba$purity)
        ## gg$set(ploidy = jabba$ploidy)
        

        ###################
        ## purity/ploidy QC plots
        ###################
        if (!check_file(report.config$cn_plot, opt$overwrite, verbose = opt$verbose)) {
            message("generating total CN purity/ploidy plots")
            pt = pp_plot(jabba_rds = opt$jabba_rds,
                         cov.fname = opt$cbs_cov_rds,
                         hets.fname = opt$het_pileups_wgs,
                         allele = FALSE,
                         field = "ratio",
                         plot.min = -2,
                         plot.max = 2,
                         bins = 100,
                         verbose = TRUE)
            ppng(print(pt), filename = report.config$cn_plot, height = 500, width = 500)
        } 

        if (check_file(report.config$allele_plot, opt$overwrite, verbose = opt$verbose)) {
            message("Allele scatter plot already exists")
        } else {
            if (file.good(opt$het_pileups_wgs)) {
                pt = pp_plot(jabba_rds = opt$jabba_rds,
                             cov.fname = opt$cbs_cov_rds,
                             hets.fname = opt$het_pileups_wgs,
                             allele = TRUE,
                             field = "count",
                             plot.min = -2,
                             plot.max = 2,
                             scatter = TRUE,
                             bins = 100,
                             verbose = TRUE)
                ppng(print(pt), filename = report.config$allele_plot, height = 500, width = 500)
            } else {
                message("No hets supplied so no allele CN plot is generated.")
            }
        }


        ## message("Checking for cohort metadata")
        ## if (file.good(opt$cohort_metadata)) {
        ##     message("Metadata found! Reading...")
        ##     meta = data.table::fread(opt$cohort_metadata)
        ## } else {
        ##     message("TPM cohort metadata not supplied.")
        ##     meta = data.table() ## empty data table to avoid existence errors
        ## }
        
        message("Checking for RNA expression input")
        if (check_file(report.config$tpm_quantiles, opt$overwrite, opt$verbose)) {
            message("RNA quantiles analysis already exists")
        } else {
            if (file.good(opt$tpm_cohort) & file.good(opt$tpm)) {

                if (opt$include_surface) {
                    surface = readRDS(report.config$surface)
                } else {
                    surface = c()
                }

                melted.expr = compute_rna_quantiles(tpm = opt$tpm,
                                                    tpm.cohort = opt$tpm_cohort,
                                                    onc = readRDS(report.config$onc),
                                                    tsg = readRDS(report.config$tsg),
                                                    surface = surface,
                                                    id = opt$pair,
                                                    gencode.gtrack = report.config$gencode_gtrack,
                                                    quantile.thresh = opt$quantile_thresh,
                                                    verbose = TRUE)
                if (!requireNamespace("matrixStats", quietly = TRUE)) {
                      stop("Package matrixStats needed for RNA data processing.")
                }
                cohort=fread(opt$tpm_cohort)
                W=data.table(gene=cohort$"gene",Avg=matrixStats::rowMeans(log10(cohort[,!("gene")]+1)),SD=rowSds(log10(as.matrix(cohort[,!("gene")]+1))))
                melted.expr$zscore=(log10(melted.expr$value+1)-W$Avg)/W$SD
            } else {
                message("RNA input not supplied.")
                melted.expr = data.table(gene = as.character(),
                                         pair = as.character(),
                                         value = as.numeric(),
                                         qt = as.numeric(),
                                         role = as.character(),
                                         direction = as.character(),
                         zscore = as.numeric())
                        
            }
            fwrite(melted.expr, report.config$tpm_quantiles)
        }
        if (check_file(report.config$rna_change, opt$overwrite, opt$verbose) &
            check_file(report.config$rna_change_all, opt$overwrite, opt$verbose)) {
            message("Over/underexpressed driver genes already identified")
        } else {
            melted.expr = data.table::fread(report.config$tpm_quantiles, header = TRUE)
            if (nrow(melted.expr)) {
                rna.change.dt = melted.expr[pair == opt$pair,][(role %like% "ONC" & qt > (1 - opt$quantile_thresh)) |
                                                               (role %like% "SURF" & qt > (1 - opt$quantile_thresh)) |
                                                               (role %like% "TSG" & qt < opt$quantile_thresh)]
                rna.change.all.dt = melted.expr[pair == opt$pair,]
            } else {
                rna.change.dt = melted.expr
                rna.change.all.dt = melted.expr
            }
            fwrite(rna.change.dt, report.config$rna_change)
            fwrite(rna.change.all.dt, report.config$rna_change_all)
        }

        message('Calling CNVs for oncogenes and tumor suppressor genes')
        if (check_file(report.config$gene_cn, overwrite = opt$overwrite, verbose = opt$verbose)){
            genes_cn_annotated = readRDS(report.config$gene_cn)
        } else {
            kag_rds = gsub("jabba.simple.rds", "karyograph.rds", opt$jabba_rds)
            nseg = NULL
            if (!is.na(opt$cbs_nseg_rds)){
                if (file.exists(opt$cbs_nseg_rds) & file.size(opt$cbs_nseg_rds) > 0){
                    if (verbose)
                        message('Reading normal segmentation copy number from: ', opt$cbs_nseg_rds)
                    nseg = readRDS(opt$cbs_nseg_rds)
                }
            } else {
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
                message('No cbs_nseg_rds was provided and the JaBbA karyograph was not found at the expected location (', kag_rds, ') so we will use CN = 2 for the normal copy number of all chromosomes.')
                }
            }

            gg = gG(jabba = opt$jabba_rds)
            genes_cn = get_gene_copy_numbers(gg,
                                             gene_ranges = opt$genes,
                                             nseg = nseg,
                                             ploidy = kag$ploidy,
                                             simplify_seqnames = (opt$ref == "hg19"),
                                             complex.fname = report.config$complex)
            
            genes_cn_annotated = get_gene_ampdel_annotations(genes_cn,
                                                             amp.thresh = opt$amp_thresh,
                                                             del.thresh = pmax(opt$del_thresh, 1))

            if (file.good(report.config$tpm_quantiles)) {

                message ("Merging SCNAs with RNA expression")
                
                melted.expr = data.table::fread(report.config$tpm_quantiles, header = TRUE)

                if (nrow(melted.expr)) {
                    genes_cn_annotated = merge.data.table(genes_cn_annotated,
                                                          melted.expr[pair == opt$pair,
                                                                      .(gene_name = gene,
                                                                        expr.quantile = qt,
                                                                        expr.value = value)],
                                                          by = "gene_name",
                                                          all.x = TRUE)
                    
                    genes_cn_annotated[expr.quantile < opt$quantile_thresh, expr := "under"]
                    genes_cn_annotated[expr.quantile > (1 - opt$quantile_thresh), expr := "over"]
                } else {
                    genes_cn_annotated[, ":="(expr.quantile = NA, expr.value = NA, expr = NA)]
                }
            } else {
                genes_cn_annotated[, ":="(expr.quantile = NA, expr.value = NA, expr = NA)]
            }
            
            ## save annotated genes
            saveRDS(genes_cn_annotated, report.config$gene_cn)
        }

        if (check_file(report.config$driver_scna,  overwrite = opt$overwrite)) {
            message("Driver gene SCNAs already identified")
        } else {

            ## read table of SCNAs
            if (file.good(report.config$gene_cn)) {
                genes_cn_annotated = readRDS(report.config$gene_cn)
            } else {
                stop("SCNA analysis failed.")
            }
                
            if (genes_cn_annotated[, .N] > 0){

                ## read oncogenes and TSGs from config
                onc = readRDS(report.config$onc)
                tsg = readRDS(report.config$tsg)

                if (!is.null(opt$include_surface) && opt$include_surface) {
                    surface = readRDS(report.config$surface)

                    driver.genes_cn = genes_cn_annotated[(cnv == "amp" & gene_name %in% onc) |
                                                         (cnv %in% c("hetdel", "homdel") & gene_name %in% tsg) |
                                                         (cnv == "amp" & gene_name %in% surface)]
                    driver.genes_cn[gene_name %in% surface, surface := TRUE]
                } else {
                    driver.genes_cn = genes_cn_annotated[(cnv == "amp" & gene_name %in% onc) |
                                                         (cnv %in% c("hetdel", "homdel") & gene_name %in% tsg)]
                }
                
                #' add whether gene is TSG or ONCO
                driver.genes_cn[gene_name %in% onc, annot := "ONC"]
                driver.genes_cn[gene_name %in% tsg, annot := "TSG"]

                fields = c("gene_name", "annot", "surface",
                           "cnv", "expr", "min_cn",
                           "max_cn", "min_normalized_cn", "max_normalized_cn",
                           "expr.value", "expr.quantile",
                           "seqnames", "start", "end", "width", "ev.id", "ev.type")
                if ('cn.low' %in% names(genes_cn_annotated) && 'cn.high' %in% names(genes_cn_annotated)){
                    # include allelic CN in the table
                    fields = c("gene_name", "annot", "surface",
                           "cnv", "expr", "min_cn", 'cn.high', 'cn.low',
                           "max_cn", "min_normalized_cn", "max_normalized_cn",
                           "expr.value", "expr.quantile",
                           "seqnames", "start", "end", "width", "ev.id", "ev.type")
                }
                
                cn.fields = intersect(fields, names(driver.genes_cn))
                fwrite(driver.genes_cn[, ..cn.fields], report.config$driver_scna)
            }
        }

        if (check_file(report.config$rna_change_with_cn, opt$overwrite, opt$verbose)) {
            message("Driver gene expression changes already identified")
        } else {

            if (file.good(report.config$rna_change)) {
                rna.change.dt = fread(report.config$rna_change, header = TRUE)
            } else {
                stop("RNA expression analysis failed")
            }
            
            if (file.good(report.config$gene_cn)) {
                genes_cn_annotated = readRDS(report.config$gene_cn)
            } else {
                stop("SCNA analysis failed")
            }

            fields = c("gene_name", "annot", "surface", "min_cn",
                       "max_cn", "min_normalized_cn", "max_normalized_cn",
                       "seqnames", "start", "end", "width", "ev.id", "ev.type")
            cn.fields = intersect(fields, colnames(genes_cn_annotated))
            
            driver.genes.expr.dt = merge.data.table(rna.change.dt[, .(gene, expr = direction,
                                                                      expr.value = value,
                                                                      expr.quantile = qt, role, zscore)],
                                                    genes_cn_annotated[, ..cn.fields],
                                                    by.x = "gene",
                                                    by.y = "gene_name",
                                                    all.x = TRUE)
            fwrite(driver.genes.expr.dt, report.config$rna_change_with_cn)
        }


        if (check_file(report.config$wgs_gtrack_plot, opt$overwrite, opt$verbose)) {
            message("Whole-genome gTrack plot exists")
        } else {
            message("Generating whole-genome gTrack plots")

            gt = wgs_gtrack(report.config$jabba_rds,
                            report.config$coverage_gtrack,
                            report.config$allele_gtrack)
            sl = intersect(names(seqlengths(gt@data[[1]])),
                           names(seqlengths(gt@data[[2]])))
            if (length(gt@data) > 2) {
                sl = intersect(sl, names(seqlengths(gt@data[[3]])))
            }
            plot.chrs = grep("(^(chr)*[0-9XY]+$)", sl, value = TRUE)
            
        if (length(plot.chrs) == 0){
                stop('None of the sequences in your genome graph matches the default set of sequences.')
            }

            if (opt$ref == "hg19") {

                std.chrs = c(as.character(1:22), "X", "Y")
                std.chrs = std.chrs[which(std.chrs %in% plot.chrs)]
            } else {

                std.chrs = paste0("chr", c(as.character(1:22), "X", "Y"))
                std.chrs = std.chrs[which(std.chrs %in% plot.chrs)]

            }

            ## temporary fix??
            ## why is there an error here?
            std.chrs = gsub("chr", "", std.chrs)
                

            ppng(plot(gt, std.chrs),
                 filename = report.config$wgs_gtrack_plot,
                 height = 1000,
                 width = 5000)
            
        } 

        if (!check_file(report.config$wgs_circos_plot, opt$overwrite, opt$verbose)) {
            message("Generating whole-genome circos plot")

            tryCatch( expr = {
                ppng(wgs.circos(jabba_rds = report.config$jabba_rds,
                                cov_fn = report.config$cbs_cov_rds,
                                transform = TRUE,
                                field = "ratio",
                                link.h.ratio = 0.1,
                                cex.points = 0.1,
                                cytoband.path = opt$cytoband,
                                chr.sub = (opt$ref == "hg19")),
                     filename = report.config$wgs_circos_plot,
                     height = 850,
                     width = 1000)
            }, error = function(e) {
                if (file.exists(report.config$wgs_circos_plot)) {
                    file.remove(report.config$wgs_circos_plot)
                }
                message(e)
                stop("An error was encountered during circos plot creation")
                return(NA)
            })
            
        } else {
            message("Whole genome circos plot already exists")
        }

        ## ##################
        ## Driver gene mutations
        ## ##################
        if (check_file(report.config$driver_mutations, opt$overwrite, opt$verbose)) {
            message("Driver mutation analysis already exists.")
        } else {

            ## ################
            ## read cgc gene list
            ## ################

            ## check if SNPEFF needs to be run for SNVs
            if (file.good(opt$snpeff_snv_bcf)) {
                message("snpEff snv supplied. Filtering for oncogenes/TSGs.")

                ## grab annotations
                snv.dt = filter.snpeff(vcf = opt$snpeff_snv_bcf,
                                       gngt.fname = report.config$gencode_gtrack,
                                       onc = report.config$onc,
                                       tsg = report.config$tsg,
                                       drivers.fname = report.config$drivers,
                                       ref.name = opt$ref,
                                       type = "snv",
                                       verbose = TRUE)

            } else {
            
                if (file.good(opt$snv_vcf)) {

                    message("Running SNV SnpEff")
                    snpeff.libdir = normalizePath(system.file("extdata", "SnpEff_module", package = "casereport"))
                    snpeff.ref = opt$ref
                    snpeff.vcf = normalizePath(opt$snv_vcf)
                    snpeff.outdir = normalizePath(file.path(opt$outdir, "snpeff", "snv"))
                    snpeff.config = normalizePath(opt$snpeff_config)

                    snpeff.cm = paste("sh", file.path(snpeff.libdir, "run.sh"),
                                      snpeff.libdir,
                                      snpeff.ref,
                                      snpeff.vcf,
                                      snpeff.outdir,
                                      snpeff.config)

                    system(snpeff.cm)

                    report.config$snpeff_snv_bcf = file.path(snpeff.outdir, "annotated.bcf")
                    saveRDS(report.config, "report.config.rds")

                    ## grab annotations
                    snv.dt = filter.snpeff(vcf = report.config$snpeff_snv_bcf,
                                           gngt.fname = report.config$gencode_gtrack,
                                           onc = report.config$onc,
                                           tsg = report.config$tsg,
                                           drivers.fname = report.config$drivers,
                                           ref.name = opt$ref,
                                           type = "snv",
                                           verbose = TRUE)
                    
                } else {
                    message("SNV vcf does not exist or is empty")
                    snv.dt = NULL
                }

            }

            ## check if SNPEFF needs to be run for indels
            if (file.good(report.config$snpeff_indel_bcf)) {
                message("snpEff snv supplied. Filtering for oncogenes/TSGs.")

                ## grab annotations
                indel.dt = filter.snpeff(vcf = report.config$snpeff_indel_bcf,
                                       gngt.fname = report.config$gencode_gtrack,
                                       onc = report.config$onc,
                                       tsg = report.config$tsg,
                                       drivers.fname = report.config$drivers,
                                       ref.name = opt$ref,
                                       type = "indel",
                                       verbose = TRUE)
            } else {

            
                if (file.good(opt$indel_vcf)) {

                    message("Running SNV SnpEff")
                    snpeff.libdir = normalizePath(system.file("extdata", "SnpEff_module", package = "casereport"))
                    snpeff.ref = opt$ref
                    snpeff.vcf = normalizePath(opt$indel_vcf)
                    snpeff.outdir = normalizePath(file.path(opt$outdir, "snpeff", "indel"))
                    snpeff.config = normalizePath(opt$snpeff_config)

                    snpeff.cm = paste("sh", file.path(snpeff.libdir, "run.sh"),
                                      snpeff.libdir,
                                      snpeff.ref,
                                      snpeff.vcf,
                                      snpeff.outdir,
                                      snpeff.config)

                    system(snpeff.cm)

                    report.config$snpeff_indel_bcf = file.path(snpeff.outdir, "annotated.bcf")
                    saveRDS(opt, "cmd.args.rds")

                    ## grab annotations
                    indel.dt = filter.snpeff(vcf = report.config$snpeff_indel_bcf,
                                           gngt.fname = report.config$gencode_gtrack,
                                           onc = report.config$onc,
                                           tsg = report.config$tsg,
                                           drivers.fname = report.config$drivers,
                                           ref.name = opt$ref,
                                           type = "indel",
                                           verbose = TRUE)
                    
                } else {
                    message("Indel vcf does not exist or is empty")
                    indel.dt = NULL
                }
            }

            ## browser()
            if (!is.null(snv.dt) & !is.null(indel.dt)) {
                driver.mutations.dt = rbind(snv.dt, indel.dt, use.names =  TRUE, fill = TRUE) %>% unique
            } else if (!is.null(snv.dt)) {
                driver.mutations.dt = snv.dt
            } else {
                driver.mutations.dt = data.table()
            }

            if (driver.mutations.dt[,.N] > 0) {
                mutation.tier.driver.name = 'PMKB' # we can later switch to a different annotator or cycle through multiple annotators if we wish to
                mutation.tier.annotator = mutation.tier.annotators[['PMKB']]
                driver.mutations.dt = mutation.tier.annotator(driver.mutations.dt)

                ## add ONC and TSG annotations
                ## notice that if the mutation.tier.parser added "gene.type" annotations then this step could override some of these
                cgc = fread(report.config$cgc)
                tsg = readRDS(report.config$tsg)## cgc[get('Is Tumor Suppressor Gene') == 'Yes', get('Hugo Symbol')]
                onc = readRDS(report.config$onc)## cgc[get('Is Oncogene') == 'Yes', get('Hugo Symbol')]
                
                ## parsing drivers
                ## this can be a text file that contains
                ## either a 1-d vector of gene names
                ## in this case, the gene.type later on is basename(opt$drivers)
                ## or it can be a 2-column table
                ## 1st column is gene name,
                ## and 2nd col is annotation for gene.type
                ## header can be commented out via "#"
                txt.ptrn = "(.tsv|.txt|.csv|.tab)(.xz|.bz2|.gz){0,}$"
                no_ext <- function (x, compression = FALSE) 
                {
                    if (compression) 
                        x <- sub("[.](gz|bz2|xz)$", "", x)
                    sub("([^.]+)\\.[[:alnum:]]+$", "\\1", x)
                }
                if (check_file(opt$drivers)) {
                    if (grepl(".rds$", opt$drivers, ignore.case = TRUE)) {
                        drivers = readRDS(opt$drivers)
                        if (!(is.character(drivers) ||
                              inherits(drivers, c("matrix", "data.frame"))))
                            drivers = NULL
                            
                    } else if (grepl(txt.ptrn, opt$drivers)) {
                        drivers = fread(opt$driver)
                    }

                    if (!is.null(drivers)) {

                        if (NCOL(drivers) == 1) {
                            ## annotating drivers by file name
                            drivers = as.data.frame(drivers)
                            drivers[[2]] = no_ext(basename(opt$drivers))  
                        } 
                        ## presuming the formatting above
                        drivers = drivers[,c(1,2), drop = F]
                        colnames(drivers) = c("gene", "driver.type")
                        driver.mutations.dt$ord345987234 = 1:NROW(driver.mutations.dt)
                        driver.mutations.dt = merge(driver.mutations.dt,
                                                    drivers,
                                                    by = "gene",
                                                    all.x = TRUE)
                        driver.mutations.dt = driver.mutations.dt[order(ord345987234)]
                        driver.mutations.dt$ord345987234 = NULL
                    }

                }


                driver.mutations.dt[gene %in% tsg, gene.type := 'TSG']
                driver.mutations.dt[gene %in% onc, gene.type := 'ONC']
                                            # some genes are annotated in CGC as both
                driver.mutations.dt[gene %in% onc & gene %in% tsg, gene.type := 'ONC|TSG']
                cols = c("gene", "gene.type", "driver.type", "Tier", "seqnames", "pos", "impact", "REF", "ALT", "variant.p", "vartype", "annotation")
                driver.mutations.dt = driver.mutations.dt[,cols[cols %in% colnames(driver.mutations.dt)],
                                    with = FALSE]
            }

            

            fwrite(driver.mutations.dt, report.config$driver_mutations)
        }

        ## ##################
        ## Fusions
        ## ##################
        if (check_file(report.config$driver_fusions, opt$overwrite, opt$verbose) &
            check_file(report.config$other_fusions, opt$overwrite, opt$verbose)) {
            message("Fusions analysis already exists!")
        } else {
            
            ## ## if opt$fusions not available, run it
            ## if (check_file(report.config$fusions, overwrite = FALSE, verbose = opt$verbose)) {
            ##     fu = readRDS(opt$fusions)
            ## } else {
                
            ##     warning('Fusions were not supplied! Computing fusions can be time/memory intensive')

            ##     if (!exists("gff")){
                    
            ##         gff = skidb::read_gencode(fn = opt$gencode)
            ##     }
                
            ##     fu = fusions(gg, gff)
            ##     saveRDS(fu, paste0(opt$outdir, "/fusions.rds"))
            ##     opt$fusions = paste0(opt$outdir, "/fusions.rds")
            ##     saveRDS(opt, "cmd.args.rds")

            ## }
                
            message("Preparing fusion genes report")
            fusions.slickr.dt = fusion.wrapper(fusions.fname = opt$fusions,
                                               complex.fname = report.config$complex,
                                               cvgt.fname = report.config$coverage_gtrack,
                                               agt.fname = report.config$allele_gtrack,
                                               cgc.fname = report.config$cgc,
                                               gngt.fname = report.config$gencode_gtrack,
                                               pad = 0.5,
                                               height = 900,
                                               width = 1000,
                                               server = opt$server,
                                               pair = opt$pair,
                                               outdir = opt$outdir)

            ## make sure data table is not empty
            if (nrow(fusions.slickr.dt) == 0) {

                ## write empty data table, rmd file will deal with this.
                fwrite(fusions.slickr.dt, report.config$driver_fusions)
                fwrite(fusions.slickr.dt, report.config$other_fusions)
                
            } else {

                ## save data table for drivers and non-drivers separately
                fwrite(fusions.slickr.dt[driver == TRUE,], report.config$driver_fusions)
                fwrite(fusions.slickr.dt[driver == FALSE,], report.config$other_fusions)

            }
        }
        
        
        ## ##################
        ## SV gallery code
        ## ##################
        if (check_file(report.config$sv_gtracks, opt$overwrite, opt$verbose)) {
            message("SV gallery files already exist")
        } else {
            message("Preparing SV gallery")

            sv.slickr.dt = gallery.wrapper(complex.fname = report.config$complex,
                                           background.fname = system.file("extdata", "sv.burden.txt", package = "casereport"),
                                           cvgt.fname = report.config$coverage_gtrack,
                                           gngt.fname = report.config$gencode_gtrack,
                                           cgcgt.fname = report.config$cgc_gtrack,
                                           agt.fname = report.config$allele_gtrack,
                                           server = opt$server,
                                           pair = opt$pair,
                                           pad = 0.5,
                                           height = 900, ## png image height
                                           width = 1000, ## png image width
                                           outdir = opt$outdir)
            
            fwrite(sv.slickr.dt, report.config$sv_gtracks)
        } 

        ## #################
        ## CN gallery
        ## #################
        if (!check_file(report.config$scna_gtracks, overwrite = opt$overwrite)) {
            message("preparing CN gallery")
            ## grab ploidy
            pl = readRDS(opt$jabba_rds)$ploidy
            cn.slickr.dt = cn.plot(drivers.fname = report.config$driver_scna,
                                   report.config$complex,
                                   cvgt.fname = report.config$coverage_gtrack,
                                   gngt.fname = report.config$gencode_gtrack,
                                   cgcgt.fname = report.config$cgc_gtrack,
                                   agt.fname = report.config$allele_gtrack,
                                   server = opt$server,
                                   pair = opt$pair,
                                   amp.thresh = opt$amp_thresh,
                                   ploidy = pl,
                                   pad = 0.5,
                                   height = 900,
                                   width = 1000,
                                   outdir = opt$outdir)
            fwrite(cn.slickr.dt, report.config$scna_gtracks)
        } else {
            message("CN gallery files already exist")
        }

        ## histograms of over and under-expressed genes
        if (check_file(report.config$expression_histograms, opt$overwrite, opt$verbose)) {
            message("expression histograms already exist")
        } else {
            expr.histograms.dt = plot_expression_histograms(report.config$rna_change,
                                                            report.config$tpm_quantiles,
                                                            pair = opt$pair,
                                                            outdir = opt$outdir,
                                                            res = 150)
            fwrite(expr.histograms.dt, report.config$expression_histograms)
        }
        
        ## gTrack of over and under-expressed genes
        if (check_file(report.config$expression_gtracks, overwrite = opt$overwrite, verbose = opt$verbose)) {
            message("gTracks of over/underexpressed genes already exist")
        } else {
            message("preparing expression gallery")
            expr.slickr.dt = cn.plot(drivers.fname = report.config$rna_change,
                                     complex.fname = report.config$complex,
                                     cvgt.fname = report.config$coverage_gtrack,
                                     gngt.fname = report.config$gencode_gtrack,
                                     cgcgt.fname = report.config$cgc_gtrack,
                                     agt.fname = report.config$allele_gtrack,
                                     server = opt$server,
                                     pair = opt$pair,
                                     amp.thresh = opt$amp_thresh,
                                     ploidy = report.config$ploidy,
                                     pad = 0.5,
                                     height = 1600,
                                     width = 1000,
                                     outdir = opt$outdir)
            fwrite(expr.slickr.dt, report.config$expression_gtracks)
        }

        ## expression waterfall plot
        if (check_file(report.config$waterfall_plot, overwrite = opt$overwrite, verbose = opt$verbose)) {
            message("Waterfall plot already exists")
        } else {
            if (file.good(opt$tpm) & file.good(opt$tpm_cohort)) {
                message("Generating waterfall plot")
                pt = rna.waterfall.plot(melted.expr.fn = report.config$tpm_quantiles,
                                        rna.change.fn = report.config$rna_change,
                                        pair = opt$pair)
                ppng(print(pt), filename = report.config$waterfall_plot,
                     width = 1600, height = 1200, res = 150)
            } else {
                message("Skipping waterfall plot - RNA expression not supplied")
            }
        } 


        ## ##################
        ## deconstructSigs composition plot
        ## ##################
        if (check_file(report.config$sig_composition, opt$overwrite, opt$verbose)) {
            message("Signatures composition plot already exists")
        } else {
            if (file.good(opt$deconstruct_sigs)) {
                message("Creating signatures composition plot from deconstructSigs input")
                sig = readRDS(opt$deconstruct_sigs)
                png(filename = report.config$sig_composition, width = 1000, height = 1000)
                deconstructSigs::plotSignatures(sig)
                dev.off()
            } else {
                message("deconstructSigs not supplied. skipping signatures composition")
            }
        }

        ######################
        ## deconstructSigs histogram
        ######################
        if (check_file(report.config$sig_histogram, overwrite = opt$overwrite, verbose = opt$verbose)) {
            message("Signatures histogram plot already exists")
        } else {

            if (file.good(opt$deconstruct_variants)) {

                ## identify correct background type
                sig.fn = system.file("extdata", "all.signatures.txt", package = "casereport")
                background.type = "Cell cohort"
                if (file.good(opt$sigs_cohort)){
                    sig.fn = opt$sigs_cohort
                    background.type = "supplied" ## what is the background dist, used for plot title
                } else {
                    if (!is.null(opt$tumor_type) & !is.na(opt$tumor_type)) {

                        ## tumor type specific signature
                        tumor.sig.fn = system.file("extdata", paste0(opt$tumor_type, ".signatures.txt"), package = "casereport")

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
            sigMet=fread(system.file("extdata", "sig.metadata.txt", package = "casereport"), sep="\t")
                sigbar = deconstructsigs_histogram(sigs.fn = opt$deconstruct_variants,
                                                   sigs.cohort.fn = sig.fn,
                                                   id = opt$pair,
                                                   cohort.type = background.type,
                            sigMet=sigMet,
                            outdir=opt$outdir)
                ppng(print(sigbar),
                     filename = report.config$sig_histogram,
                     height = 800, width = 800, res = 150)

            presentSigs=fread(file.path(opt$outdir,"Sig.csv"))
                presentSigs$pair=NULL
                thisMet=sigMet[sigMet$Signature %in% presentSigs$Signature,]
                metTable=merge(presentSigs,thisMet,by="Signature")
            metTable=metTable[,c("Signature","Mutational.process","sig_count","perc")]
            colnames(metTable)=c("Signature","Mutational.process","sig_count","quantile")
                metTable=metTable[order(-metTable$quantile), ]
                fwrite(metTable, file.path(opt$outdir,"signatureMetadata.csv"))

                
            } else {
                message("deconstructSigs output not supplied.")
            }
        }
                

        ## ##################
        ## HRDetect results
        ## 
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

                hrd = merge.data.table(hrd.out, hrd.dat[, .(variable, data = value)], by = "variable")

                ## original training data
                hrd_cohort = fread(system.file("extdata", "/data/hrdetect.og.txt", package = "casereport"))
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
                    theme_minimal() +
                    theme(legend.position = "bottom",
                          legend.text = element_text(size=20),
                          title = element_text(size = 20, family = "sans"),
                          axis.title = element_text(size = 20, family = "sans"),
                          axis.text.x = element_text(size = 15, family = "sans"),
                          axis.text.y = element_text(size = 20, family = "sans"))

                ## finally draw the output probability in a waterfall plot
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
        ## Oneness / Twoness
        ## 
        ## ##################
        all.ot = check_file(report.config$ot_log, opt$overwrite) &
            check_file(report.config$ot_prob, opt$overwrite) &
            check_file(report.config$oneness, opt$overwrite) &
            check_file(report.config$twoness, opt$overwrite)
        if (!all.ot){
            if (file.good(opt$ot_results)) {
                ot.res = readRDS(opt$ot_results)
                ot = merge(ot.res$expl_variables,
                           ot.res$ot_scores)
                ot$`.` = NULL
                ## ot = data.table::melt(ot)
                ot$pair = opt$pair
                tx = paste0("ot[,c(", paste(!names(ot) %in% c("qrdup", "qrdel", "tib"), collapse = ","), "),drop=F]")
                ot = eval(parse(text = tx)) ## this is just for robustness
                ## ## melt the input data matrix
                ## ot.dat = as.data.table(ot.res$data_matrix) %>% data.table::melt()
                ## ot.dat[, pair := opt$pair]
                ## ## melt the output results
                ## ot.out = as.data.table(ot.res$expl_variables)
                ## ot.out$`.` = NULL
                ## ot.out = data.table::melt(ot.out)
                ## ot.out[, pair := opt$pair]
                ## ot.out = ot.out[!variable %in% c("qrdup", "qrdel", "tib")]
                saveRDS(ot, paste0(opt$outdir, "/ot.rds"))

                ## ot = merge(ot.out, ot.dat[, .(variable, data = value)], by = "variable")

                ## original training data
                ot_cohort = readRDS(system.file("extdata", "/data/ot_scores_cohort.rds", package = "casereport"))

                ## if (opt$pair %in% hrd_cohort$pair){
                ot_cohort = rbind(ot_cohort[!pair %in% opt$pair], ot, fill = T)
                ## }

                ## ot_cohort[, range(value), by = .(variable, is.hrd)]

                ot_cohort[, ot.status := ifelse(BRCA1 > 0.5, "Oneness",
                                         ifelse(BRCA2 > 0.5, "Twoness",
                                         ifelse(SUM12 > 0.9, "OT+", "OT-")))]

                ot.stat = ot_cohort[ot_cohort$pair %in% opt$pair]$ot.status

                ot.stat = c(
                    "Oneness" = "#C93312",
                    "Twoness" = "#3B9AB2",
                    "OT+" = "#0B775E",
                    ## "OT-" = "#899DA4"
                    "OT-" = "#899DA480"
                )

                ot_dat = data.table::melt(ot_cohort, id.vars = c("pair", "ot.status",
                                                                 "fmut_bi"))

                

                ## these are the dimensions that need to be logged
                ldat = ot_dat[!variable %in% c("del.mh.prop", "BRCA1", "BRCA2", "SUM12", "OTHER")]
                ldat[, variable := factor(variable, levels = setdiff(levels(variable), "del.mh.prop"))]
                ## plot the distribution of raw count
                ot.plot = ggplot(ldat, aes(x = value, y = variable)) +
                    geom_density_ridges(aes(point_colour = ot.status,
                                            fill = ot.status),
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
                    scale_discrete_manual("point_colour", values = ot.stat) +
                    scale_fill_manual(values = ot.stat) +
                    geom_segment(data = ldat[pair==opt$pair],
                                 aes(x = value, xend = value,
                                     y = as.numeric(variable), yend = as.numeric(variable) + 0.9, colour = ot.status),
                                 alpha = 1, lwd = 3) +
                    scale_colour_manual(values = ot.stat) + 
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
                rdat = ot_dat[variable %in% "del.mh.prop"]
                rdat[, variable := factor(variable, levels = "del.mh.prop")]
                ## plot the distribution of raw count
                ot.plot.2 = ggplot(rdat, aes(x = value, y = variable)) +
                    geom_density_ridges(aes(point_colour = ot.status,
                                            fill = ot.status),
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
                    scale_discrete_manual("point_colour", values = ot.stat, labels = NULL) +
                    scale_fill_manual(values = ot.stat,
                                      labels = c("Non HR deficient", "HR deficient")) +
                    ## scale_y_discrete(labels = hrdetect.dims) +
                    labs(x = "Proportion") +
                    guides(point_colour = FALSE) +
                    geom_segment(data = rdat[pair==opt$pair],
                                 aes(x = value, xend = value,
                                     y = as.numeric(variable), yend = as.numeric(variable) + 3, colour = ot.status),
                                 alpha = 1, lwd = 3) +
                    scale_colour_manual(values = ot.stat) + 
                    theme_minimal() +
                    theme(legend.position = "bottom",
                          legend.text = element_text(size=20),
                          title = element_text(size = 20, family = "sans"),
                          axis.title = element_text(size = 20, family = "sans"),
                          axis.text.x = element_text(size = 15, family = "sans"),
                          axis.text.y = element_text(size = 20, family = "sans"))

                ## finally draw the output probability in a waterfall plot

                prob_dat1 = ot_dat[variable %in% c("BRCA1")]
                prob_dat1$marked_pair = prob_dat1$pair %in% opt$pair
                prob_dat1[, fpair := reorder(pair, value)]


                prob.cols = binary.cols
                prob.cols["TRUE"] = ifelse(prob_dat1[pair == opt$pair]$value > 0.5,
                       binary.cols["TRUE"], "black")
                
                prob.brca1 = ggplot(prob_dat1, aes(x = fpair, y = value)) +
                    geom_bar(aes(color = marked_pair), stat = "identity") +
                    scale_colour_manual(values = prob.cols) +
                    ylab("Oneness") +
                    theme_minimal() +
                    theme(legend.position = "none",
                          legend.text = element_text(size=20),
                          title = element_text(size = 20, family = "sans"),
                          axis.title = element_text(size = 20, family = "sans"),
                          axis.text.x = element_blank(),
                          axis.title.x = element_blank(),
                          axis.text.y = element_text(size = 20, family = "sans"))

                
                prob_dat2 = ot_dat[variable %in% c("BRCA2")]
                prob_dat2$marked_pair = prob_dat2$pair %in% opt$pair
                prob_dat2[, fpair := reorder(pair, value)]


                prob.cols = binary.cols
                prob.cols["TRUE"] = ifelse(prob_dat2[pair == opt$pair]$value > 0.5,
                       binary.cols["TRUE"], "black")
                
                prob.brca2 = ggplot(prob_dat2, aes(x = fpair, y = value)) +
                    geom_bar(aes(color = marked_pair), stat = "identity") +
                    scale_colour_manual(values = prob.cols) +
                    ylab("Twoness") + 
                    theme_minimal() +
                    theme(legend.position = "none",
                          legend.text = element_text(size=20),
                          title = element_text(size = 20, family = "sans"),
                          axis.title = element_text(size = 20, family = "sans"),
                          axis.text.x = element_blank(),
                          axis.title.x = element_blank(),
                          axis.text.y = element_text(size = 20, family = "sans"))

                
                
                ## draw the plots
                png(report.config$ot_log, width = 800, height = 800)
                print(ot.plot)
                dev.off()

                png(report.config$ot_prob, width = 800, height = 200)
                print(ot.plot.2)
                dev.off()

                png(report.config$oneness, width = 800, height = 200)
                print(prob.brca1)
                dev.off()

                png(report.config$twoness, width = 800, height = 200)
                print(prob.brca2)
                dev.off()
            } else {
                message("Oneness/twoness data not provided")
            }
        } else {
            message("Oneness/twoness analysis already exists.")
        }

        ## ##################
        ## Enhancer/gene proximity
        ## ##################
        proximity.gallery.fname = paste0(opt$outdir, "/", "proximity.gallery.txt")
        if (check_file(proximity.gallery.fname, overwrite = opt$overwrite, verbose = opt$verbose)) {
            message("Proximity gallery already exists")
        } else {
            if (!file.good(report.config$rna_change)) {
                message("No RNA analysis, skipping proximity")
                proximity.gallery.dt = data.table(gene = as.character(),
                                                  plot.fname = as.character())
            } else if (!file.good(opt$proximity)) {
                message("No proximity.rds supplied - skipping!")
                proximity.gallery.dt = data.table(gene = as.character(),
                                                  plot.fname = as.character())
            } else {
        
                prox = readRDS(opt$proximity)

                if (length(prox)) {
                    pdt = prox$dt
                    pdt = pdt[reldist<0.25 & refdist>5e6]
                } else {
                    pdt = data.table(gene_name = c())
                }
                cool.exp = fread(report.config$rna_change) %>% setkey("gene")

                if (any(cool.exp[direction=="over", gene] %in% pdt$gene_name)) {
                    
                    ## filter by overexpressed genes
                    pdt = merge.data.table(pdt, cool.exp[direction=="over"], by.x = "gene_name", by.y = "gene", all.x = TRUE)
                    pgs = pdt[(direction=="over"), unique(gene_name)]
                    gt.cgc = readRDS(report.config$cgc_gtrack)
                    gt.cgc$name = "CGC"
            gt.cgc$height=15
                    gt.cgc$cex.label=1.2
            enh = readRDS(opt$enhancer)

                    ## read some gTracks
                    cvgt = readRDS(report.config$coverage_gtrack)

                    ## read gGraph
                    gg = readRDS(report.config$complex)

                    proximity.fns = lapply(setNames(nm = pgs),
                                           function(g) {
                                               this.png = paste0(opt$outdir, "/", g, ".proximity.png")
                                               w = prox[pdt[gene_name==g, walk.id]]
                                               this.enh = copy(enh)
                                               this.enh$col = binary.cols[as.character(seq_along(enh) %in% w$dt$qid)]
                                               gt.enh = gTrack(this.enh, name = "enhancer")
                           cvgt$xaxis.cex.label=2.0
                           cvgt$xaxis.cex.tick=2.0
                                               pgt = c(cvgt, gg$gtrack(height = 30),
                                                       w$gtrack(name = "shortest walk", height=15), gt.cgc, gt.enh)
                                               png(this.png, height = 1200, width = 1800, pointsize = 20)
                                               plot(pgt, w$footprint + 1e6, legend.params = list(plot = FALSE), y0 = 0, mai=c(0,0,0,0))
                                               dev.off()
                                               return(this.png)
                                           }) %>% unlist
                    proximity.gallery.dt = data.table(gene = names(proximity.fns),
                                                      plot.fname = proximity.fns)

                } else {
                    proximity.gallery.dt = data.table(gene = as.character(),
                                                      plot.fname = as.character())
                }
            }
            fwrite(proximity.gallery.dt, proximity.gallery.fname)
        }


        ## #################
        ## summarize
        ## #################
        if (check_file(report.config$summary_stats, opt$overwrite, opt$verbose)) {
            message("Summary already exists, skipping")
        } else {
            message("Computing CN/mutation summary")
            summary.list = create.summary(jabba_rds = opt$jabba_rds,
                                          snv_vcf = opt$snv_vcf,
                                          indel_vcf = opt$indel_vcf,
                                          verbose = TRUE,
                                          amp.thresh = opt$amp_thresh,
                                          del.thresh = opt$del_thresh)
            saveRDS(summary.list, report.config$summary_stats)
        }

        
        ## #################
        ## deconv
        ## #################
        if (check_file(report.config$deconv, opt$overwrite, opt$verbose)) {
          message("Deconvolution data already exists, skipping")
        } else {
          if (file.good(opt$tpm)){
            message("Running Deconvolution algorithm")
            tpm_raw = as.character(opt$tpm)
            tpm_read <- data.table::fread(tpm_raw, header = TRUE)
            tpm_read_new <- tpm_read[,-1]
            tpm_read_new_name <- as.matrix(tpm_read[,1])
            rownames(tpm_read_new) <- tpm_read_new_name[,1] 
            deconv_results = immunedeconv::deconvolute(tpm_read_new, opt$deconv)
            data.table::fwrite(deconv_results, file.path(opt$outdir,"deconv_results.txt"), sep = '\t', quote = F, row.names = F)
          }
        }
        
        
        ## ################
        ## create oncotable
        ## ################

        if (check_file(report.config$oncotable, opt$overwrite, opt$verbose)) {
            message("Oncotable already exists. Skipping!")
        } else {
            message("Preparing oncotable...")
            oncotable.input = data.table(pair = opt$pair,
                                         jabba_rds = opt$jabba_rds,
                                         fusions = opt$fusions,
                                         complex = report.config$complex,
                                         scna = report.config$gene_cn,
                                         annotated_snv_bcf = opt$snpeff_snv_bcf,
                                         annotated_indel_bcf = opt$snpeff_indel_bcf,
                                         rna = report.config$rna_change_all,
                                         proximity = opt$proximity,
                                         deconstruct_variants = opt$deconstruct_variants,
                         hrd_results = opt$hrd_results,
                                         key = "pair")
            oncotable = oncotable(oncotable.input,
                                  gencode = opt$genes,
                                  verbose = TRUE)
            saveRDS(oncotable, report.config$oncotable)
        } 

        ## ################
        ## create summaryTable
        ## ################


       if (check_file(report.config$summaryTable, opt$overwrite, opt$verbose)) {
            message("Summary Table already exists. Skipping!")
        } else {
            message("Generating summary table")
        wol=makeSummaryTable(report.config$driver_scna,report.config$driver_fusions, report.config$rna_change_with_cn,report.config$driver_mutations,report.config$oncotable,opt$libdir)
        wol$tier=as.character(wol$tier)
        wol[is.na(wol$tier),]$tier="Undefined"
        fwrite(wol,report.config$summaryTable)
        }


    }

    message("Optimizing PNGs")
    cm = paste("python",
               paste0(opt$libdir, "/", "optimize_png.py"),
               "--dirname",
               opt$outdir)
    system(cm)


    message("Start knitting")
    rmarkdown::render(
        input = normalizePath(paste0(opt$libdir, "/wgs.report.rmd")),
        output_format = "html_document",
        output_file = normalizePath(paste0(opt$outdir, "/", opt$pair,".wgs.report.html")),
        knit_root_dir = normalizePath(opt$outdir),
        ## params = report.config,
        params = list(set_title = paste0(opt$pair),
                      pair = opt$pair,
                      jabba_rds = normalizePath(opt$jabba_rds),
                      outdir = normalizePath(opt$outdir),
                      tumor_type = opt$tumor_type,
                      server = opt$server),
        quiet = FALSE)

    message("WGS case report completed.")
}

casereport = function(thisp){

    message(paste0("Beginning Case Analysis for ", thisp$pair))

    onepline = "/gpfs/commons/groups/imielinski_lab/git/mskilab/flows/modules/SnpEff/source/snpEff/scripts/vcfEffOnePerLine.pl"

    library(gGnome)
    dt.patho = gr2dt(readRDS("/gpfs/commons/groups/imielinski_lab/DB/ClinVar/pathogenic.rds"))

    filt = "java -Xmx1g -Xms1g -jar /gpfs/commons/groups/imielinski_lab/git/mskilab/flows/modules/SnpEff/source/snpEff/SnpSift.jar filter \"( ANN =~ 'missense|splice|stop_gained|frame' )\""

    gff = readRDS(gzcon(url('http://mskilab.com/gGnome/hg19/gencode.v19.annotation.gtf.gr.rds')))
    ge = read_gencode('gene'); message('loaded GENCODE')
    pge = ge %Q% (gene_type == 'protein_coding' & transcript_status == "KNOWN");
    cds = read_gencode('cds'); message('loaded CDS')
    path = paste0(getwd(),"/casereports/",thisp$pair,"/")
    pge = genes %Q% (gene_type == 'protein_coding') %Q% (level <3 & transcript_status == 'KNOWN')
    cgc = fready('/gpfs/commons/groups/imielinski_lab/DB/COSMIC/CGC.txt')
    cgc$onc = cgc$Significant_across_any_tumor_type_lineage_if_gain_ == 'Yes'
    cgc$tsg = cgc$Significant_across_any_tumor_type_lineage_if_lost_ == 'Yes'
    pge.cgc = pge %Q% (gene_name %in% cgc[onc == TRUE | tsg == TRUE, Gene_Symbol])
    cgc_tsg = pge %Q% (gene_name %in% cgc[tsg == TRUE, Gene_Symbol])
    cgc_onc = pge %Q% (gene_name %in% cgc[onc == TRUE, Gene_Symbol])
    cgc_fus = pge %Q% (gene_name %in% cgc[nchar(Translocation_Partner) > 0, Gene_Symbol])
    Sys.setenv(DEFAULT_BSGENOME = "/gpfs/commons/groups/imielinski_lab/DB/UCSC/hg19.broad.chrom.sizes")

    message(paste0("Returning Outputs to Directory: ", path))
    system(paste0("mkdir -p ",path))

    ##Load and Print COSMIC Signatures
    if(!is.na(thisp$signature_counts)){
        message("Returning top 10 COSMIC Signature Counts")
        sig = try(fread(thisp$signature_counts))
        cols = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a") ##Make a color map real quick to match to the names that I will want to color by

        names(cols) = sig$Signature[1:10]
        barplot = ggplot(sig[1:10],
                         aes(x = Signature, y = num_events)) +
            geom_col(aes(fill = Signature)) +
            ylim(0,max(sig$Signature)) +
            scale_fill_manual(values = cols,aesthetics = "fill") +
            labs(fill = " ",
                 x = "Cosmic SNV Signature",
                 y = "SNV Contributing Counts",
                 title = paste0("COSMIC SNV Signature Contributions to ", thisp$pair[1]),
                 size = 20) +
            theme(legend.position = "none",
                  axis.text.x = element_text(size = 10,
                                             angle = 45,
                                             vjust = 0.5),
                  axis.text.y = element_text(size = 15)
                  )
        ppng(print(barplot), filename = paste0(path,"COSMIC_signatures.png"))
    } else {
        message(paste0("No Signatures generated for ",thisp$pair))
    }

    ##Load and print CBS Ratio Coverage
    kar = karyogram()
    message("Defining covcbs and covgt functions")
    covcbs = function(x){gTrack(readRDS(x),
                                y.field = "ratio",
                                name = " ",
                                max.ranges = 1e4,
                                circles = TRUE,
                                lwd.border = 0.2)}
    covgts = function(x){gTrack(readRDS(x),
                                y.field = "seg.mean",
                                name = " ",
                                max.ranges = 1e4,
                                lwd.border = 0.2)}
    covgtn = function(x){gTrack(readRDS(x),
                                y.field = "mean",
                                name = " ",
                                max.ranges = 1e4,
                                lwd.border = 0.2)}

    if(!is.na(thisp$cbs_cov_rds)){
        message("Loading gTracks")
        cbs.gt = covcbs(thisp$cbs_cov_rds)
        cbs.gt$name = "CBS Ratio"
        cov.t.gt = covgts(thisp$cbs_seg_rds)
        cov.t.gt$name = "Segmentation "
        cov.n.gt = covgtn(thisp$cbs_nseg_rds)
        cov.n.gt$name = "Normal Segment Mean"
        win = si2gr(hg_seqlengths())[c(-25,-26)]
        message("Printing Coverage Plots")
        ppdf(
            lapply(1:length(win),function(x){
                plot.new()
                plot(c(kar,cov.n.gt,cov.t.gt,cbs.gt),
                     win[x])}),
            width = 20,
            title = paste0("Genome Scale Coverage for ", thisp$pair),
            filename = paste0(path,"Coverage_profile.pdf"))
        message("Done Printing Coverage")
    } else {
        message(paste0("No CBS data for ", thisp$pair))
    }

    ##Load and print JaBbA findings
    if(!is.na(thisp$jabba_rds)){
        message("Printing JaBbA Genome Scale View")
        gg = gG(jabba = thisp$jabba_rds)
        message("Loaded JaBbA Object")
        gt = gg$gt
        cbs.gt = covcbs(thisp$cbs_cov_rds)
        cbs.gt$name = "CBS Ratio"
        ppdf(lapply(1:length(win),function(x){
            plot.new()
            plot(c(kar,cbs.gt,gt),
                 win[x])}),
            width = 20,
            title = paste0("JaBbA Genome Scale Model: ", thisp$pair),
            filename = paste0(path,"JaBbA.global.pdf"))
        message("Done printing JaBbA genome scale")
        message("Now printing complex events")
        message("Identifying and printing chromothripsis events")
        gg = chromothripsis(gg)
        ct = gg[chromothripsis > 0]$copy$subgraph(d = 1e4)
        if(length(ct) != 0){
            ppdf(plot(c(kar,cbs.gt,gt),ct$footprint+1e4),
                 width = 15,
                 title = paste0("Chromothripsis Events in: ", thisp$pair),
                 filename = paste0(path,"JaBbA.chromothripsis.pdf"))
        } else {message("No Chromothripsis found")}
        message("Identifying and printing bfb events")
        gg = bfb(gg)
        bfb = gg[bfb > 0]$copy$subgraph(d = 1e4)
        if(length(bfb) != 0){
            ppdf(plot(c(kar,gt.ge,cbs.gt,gt),gg$nodes[bfb>0]$edges$footprint+1e4),
                 width = 15,
                 title = paste0("Breakage Fusion Bridges in: ", thisp$pair),
                 filename = paste0(path,"JaBbA.bfb.pdf"))
        } else {message("No BFBs found")}
        message("Identifying and printing Templated insertion Chains")
        gg = tic(gg)
        tic = gg[tic > 0]$copy$subgraph(d = 1e4)
        if(length(tic) != 0){
            ppdf(plot(c(kar,cbs.gt,gt),
                      gg$nodes[tic > 0]$edges$footprint +1e4),
                 width = 15,
                 title = paste0("Templated Insertion Chains in: ", thisp$pair),
                 filename = paste0(path,"JaBbA.tic.pdf"))
        } else {message("No TICs found")}
        message("Identifying and printing Pyrgos")
        gg = pyrgos(gg)
        pyrgo = gg[pyrgo > 0]$copy$subgraph(d = 1e4)
        if(length(pyrgo) != 0){
            ppdf(plot(c(kar,cbs.gt,gt),pyrgo$footprint+1e4),
                 width = 15,
                 title = paste0("Pyrgos in: ", thisp$pair),
                 filename = paste0(path,"JaBbA.pyrgos.pdf"))
        } else {message("No Pyrgos found")}
        message("Identifying and printing Rigmas")
        gg = rigma(gg)
        rigma = gg[rigma > 0]$copy$subgraph(d = 1e4)
        if(length(rigma) != 0){
            ppdf(plot(c(kar,gt.ge,cbs.gt,gt),gg$nodes[rigma > 0]$edges$footprint+1e4),
                 width = 15,
                 title = paste0("Rigmas in: ", thisp$pair),
                 filename = paste0(path,"JaBbA.rigma.pdf"))
        } else {message("No Rigmas found")}
        message("Calling Fusions")
        fus = gGnome::fusions(gg,gff)
        if(length(fus) != 0){
            dt = fus$dt[in.frame == TRUE][silent == FALSE]
            dt = dt[ ,c(1:19,22:23)]
            write.csv(dt, paste0(path, "fusions.summary.csv"))
            fus = gGnome::fus[!duplicated(fus$dt$name)]
            message("Printing Fusions")
            ppdf(lapply(1:length(fus),function(x){
                plot.new()
                plot(c(cbs.gt,gt,fus$gt),fus[x]$footprint+1e5)}),
                width = 15,
                title = paste0("Fusions in: ", thisp$pair),
                filename = paste0(path,"Fusion.plots.pdf")
                )
        } else {message("No fusions found")}
        message("Calling CNAs from JaBbA Output")
        seg = gg$nodes$dt[, pair := thisp$pair]
        seg[, ploidy := sum(width*cn, na.rm = TRUE)/sum(width[!is.na(cn)], na.rm = TRUE), by = pair]
        seg[, cn.rel := cn/ploidy]
        seg[, amp := cn.rel>2]
        seg[, homdel := cn==0]
        seg[, hetdel := cn==1]
        seg.gr = dt2gr(seg[!is.na(cn), ])
        seg.gr.onc = cgc_onc[, "gene_name"] %*% (seg.gr %Q% (amp == TRUE)) %Q% (!duplicated(cbind(gene_name, pair)))
        seg.gr.onc$type = 'amp'
        seg.gr.tsg = cgc_tsg[, "gene_name"] %*% (seg.gr %Q% (hetdel | homdel)) %Q% (!duplicated(cbind(gene_name, pair)))
        seg.gr.tsg$type = ifelse(seg.gr.tsg$homdel, 'homdel', ifelse(seg.gr.tsg$hetdel, 'hetdel', NA))
        cna = seg2gr(seg[, alt := ifelse(amp, 'amp', ifelse(hetdel | homdel, 'del', 'wt'))][alt != 'wt', ])[, c('cn', 'cn.rel', 'alt', 'pair')]
        cnag = cna %*% (pge.cgc[, 'gene_name'])
        dt.cna = gr2dt(cnag)
        write.csv(dt.cna, paste0(path, "cna.genes.csv"))
    } else {
        message("No JaBbA complete for this sample")
    }

                                        #Parse and return SNPEFF calls
    if(!is.na(thisp$annotated_vcf_somatic)){
        message("Parsing SNPEFF Somatic Calls")
        vcf_path = thisp$annotated_vcf_somatic  ##Fix your input
        out.name = paste0("tmp_", rand.string(), ".vcf")
        tmp.path = paste0(tempdir(), "/", out.name)
        cmd = sprintf("cat %s | %s | %s > %s", vcf_path, onepline, filt, tmp.path)
        system(cmd)
        fn = gsub('\\.vcf$', 'coding.somatic.rds', vcf_path)  ##Choose where to write it
        vcf = grok_vcf(tmp.path, long = TRUE)
        vcf = vcf %Q% (!impact %in% c("LOW", "MODIFIER"))
        unlink(tmp.path)
        dt = gr2dt(vcf)
        dt = dt[,-c("RU","LOF","NMD")]
        write.csv(dt, paste0(path,"Somatic.mutations.csv"))
        message("Parsing SNPEFF Germline Calls")
        vcf_path = thisp$annotated_vcf_germline  ##Fix your input
        out.name = paste0("tmp_", rand.string(), ".vcf")
        tmp.path = paste0(tempdir(), "/", out.name)
        cmd = sprintf("cat %s | %s | %s > %s", vcf_path, onepline, filt, tmp.path)
        system(cmd)
        fn = gsub('\\.vcf$', 'coding.germline.rds', vcf_path)  ##Choose where to write it
        vcf = grok_vcf(tmp.path, long = TRUE)
        vcf = vcf %Q% (!impact %in% c("LOW", "MODIFIER"))
        unlink(tmp.path)
        dt.g = gr2dt(vcf)
        dt.g = merge(dt.g, copy(dt.patho)[, pathogenic := TRUE],
                     by = c("seqnames","start","end","REF"),
                     all.x = TRUE,
                     allow.cartesian = TRUE)[, pathogenic := replace_na(pathogenic, FALSE)]
        dt.g = dt.g[pathogenic == TRUE]
        dt.g = dt.g[,-c("RU","LOF","NMD")]
        write.csv(dt.g, paste0(path,"Germline.mutations.csv"))
        dt.union = rbind(dt, dt.g, fill = TRUE)
        dt.union = dt.union[,-c("RU","LOF","NMD")]
        write.csv(dt.union, paste0(path, "All.mutations.csv"))
    } else{
        message("No SNPEFF Annotations Available")}
}



