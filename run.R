library(optparse)
options(bitmapType='cairo')

if (!exists('opt'))
{
    option_list = list(
        make_option(c("--libdir"), type = "character", help = "libdir"),
        make_option(c("--pair"), type = "character", help = "pair.id"),
        make_option(c("--tumor_bam"), type = "character", help = "tumor bam"),
        make_option(c("--jabba_rds"), type = "character", help = "jabba output"),
        make_option(c("--complex"), type = "character", help = "complex event caller"),
        make_option(c("--fusions"), type = "character", help = "fusions module output"),
        make_option(c("--proximity"), type = "character", help = "proximity module output"),
        make_option(c("--cbs_cov_rds"), type = "character", help = "cbs_cov output"),
        make_option(c("--deconstruct_sigs"), type = "character", help = "deconstruct_sigs"),
        make_option(c("--strelka2_somatic_variants"), type = "character", help = "strelka2_somatic_variants"),
        make_option(c("--svaba_somatic_vcf"), type = "character", help = "svaba_somatic_vcf"),
        make_option(c("--annotated_vcf"), type = "character", help = "snpeff annotated somatic vcf file"),
        make_option(c("--annotated_vcf_germline"), type = "character", help = "snpeff annotated germline vcf file"),
        make_option(c("--gencode"), type = "character", default = "~/DB/GENCODE/hg19/gencode.v19.annotation.gtf",
                    help = "GENCODE gene models in GTF/GFF3 formats"),
        make_option(c("--chrom_sizes"), type = "character", default = "/gpfs/commons/groups/imielinski_lab/DB/UCSC/hg19.broad.chrom.sizes",
                    help = "chrom.sizes file of the reference genome"),
        make_option(c("--chr"), type = "logical", default = FALSE, action = "store_true", help = "Force the reference genome to have chr prefix"),
        ## make_option(c("--germline_muts"), type = "character", help = "alterations module output"),
        ## make_option(c("--germline_truncs"), type = "character", help = "alterations module output"),
        ## make_option(c("--somatic_muts"), type = "character", help = "alterations module output"),
        ## make_option(c("--somatic_truncs"), type = "character", help = "alterations module output"),
        ## make_option(c("--homdels"), type = "character", help = "alterations module output"),
        ## make_option(c("--hetdels"), type = "character", help = "alterations module output"),
        ## make_option(c("--amplifications"), type = "character", help = "alterations module output"),
        ## make_option(c("--mutations"), type = "character", help = "alterations module output"),
        ## make_option(c("--proteins"), type = "character", help = "alterations module output"),
        make_option(c("--knit_only"), type = "logical", default = FALSE, action = "store_true",
                    help = "if true, skip module and just knit"),
        make_option(c("--outdir"), type = "character", default = './', help = "Directory to dump output into"),
        make_option(c("--amp_thresh"), type = "numeric", default = 4,
                    help = "Threshold over ploidy to call amplification"),
        make_option(c("--del_thresh"), type = "numeric", default = 0.5,
                    help = "Threshold over ploidy to call deletion"),
        make_option(c("--server"), type = "character", default = "https://mskilab.com/gGraph/", help = "URL of the gGnome.js browser"),
        make_option(c("--tumor_type"), type = "character", default = "", help = "tumor type"),
        make_option(c("--somatic_snv_cn"), type = "character", help = "MLE somatic SNV CN estimate from snv_multiplicity2 module"),
        make_option(c("--oncokb_token"), default = '~/keys/oncokb.token', type = "character", help = "a token to use when querying OncoKB. If you don't have a token, you must first obtain one from: https://www.oncokb.org/apiAccess. By default looking for the key in '~/keys/oncokb.token'. If no valid key is provided then OncoKB annotations will be skipped.")
        make_option(c("--germline_snv_cn"), type = "character", help = "MLE germline SNV CN estimate from snv_multiplicity2 module"),
        make_option(c("--overwrite"), type = "logical", default = FALSE, action = "store_true", help = "overwrite existing data in the output dir")
    )
    parseobj = OptionParser(option_list=option_list)
    opt = parse_args(parseobj)
    opt$outdir = normalizePath(opt$outdir)
    saveRDS(opt, paste(opt$outdir, 'cmd.args.rds', sep = '/'))
}
message(normalizePath(opt$outdir))

message("Loading Libraries -- Please wait...")
suppressMessages(expr = {
    suppressPackageStartupMessages(expr = {
        library(skidb)
        library(stringr)
        ## library(skitools)
        library(gGnome)
        library(gTrack)
        ## library(naturalsort)
        library(wesanderson)
        library(gUtils)
        library(tidyr)
        library(dplyr)
        library(ggforce)
        library(httr)
        library(jsonlite)
        library(knitr)
        library(rmarkdown)
        message("Loading critical dependencies from KevUtils")
        source(paste0(opt$libdir, "/utils.R"))
        source(paste0(opt$libdir, "/config.R"))
    })
})
options(warn=-1)
opt$outdir = normalizePath(opt$outdir)

## source("/gpfs/commons/groups/imielinski_lab/home/khadi/DThelp.R")
## source("/gpfs/commons/groups/imielinski_lab/home/khadi/hkev_utils.R")
## source('/gpfs/commons/groups/imielinski_lab/home/khadi/safety_func.R')

## no regenerating the data
if(opt$knit_only == FALSE){
    ## message("knit.only set to true, attempting knit")
    ## out = normalizePath(opt$outdir)
    ## opt$outdir = out
    ## rmarkdown::render(input = "/gpfs/commons/home/jrosiene/notebooks/report.rmd",
    ##                   output_file = paste0(opt$outdir, "/", opt$pair,".report.html"),
    ##                   knit_root_dir = normalizePath(opt$outdir),
    ##                   params = list(set_title = paste0(opt$pair),
    ##                                 pair = opt$pair,
    ##                                 outpath = out,
    ##                                 bam = opt$tumor_bam),
    ##                   quiet = FALSE)
    ## quit()
    ## all the reference data are temporarily stored in module directory
    ## TODO: make it customizable
    message("Setting system DEFAULT_BSGENOME variable to hg19")
    Sys.setenv(DEFAULT_BSGENOME = opt$chrom_sizes)

    message("Loading reference genome")
    sl_fn = paste0(opt$outdir, "/sl.rds")
    if (file.exists(sl_fn) & file.size(sl_fn) > 0 & !opt$overwrite){
        sl = readRDS(sl_fn)
    } else {
        sl = hg_seqlengths(chr = opt$chr)
        saveRDS(sl, sl_fn)
    }
    hg_fn = paste0(opt$outdir, "/hg.rds")
    if (file.exists(hg_fn) & file.size(hg_fn) > 0 & !opt$overwrite){
        hg = readRDS(hg_fn)
    } else {
        hg = si2gr(sl) %>% gr.stripstrand %>% unname
        saveRDS(hg, hg_fn)
    }

    message("Loading Gencode genes")
    ## gff = readRDS(gzcon(url('http://mskilab.com/gGnome/hg19/gencode.v19.annotation.gtf.gr.rds')))
    gff = NA # we assign NA here so that we can later only load it when needed (for example for fusions calls below)
    pge_fn = paste0(opt$outdir, "/pge.rds")
    if (check_file(pge_fn, overwrite = opt$overwrite)){
        pge = readRDS(pge_fn)
    } else {
        gff = skidb::read_gencode(fn = opt$gencode)
        pge = gff %Q% (type=="gene") %Q% (gene_type=="protein_coding") %Q% (level<3)
        saveRDS(pge, pge_fn)
    }
    ## ge = readRDS(paste0(opt$libdir,"/db/ge.rds"))
    ## pge = readRDS(paste0(opt$libdir,"/db/pge.rds"))
    ## ug = readRDS(paste0(opt$libdir,"/db/ug.rds"))
    ## gco2 = split(gff, gff$gene_name)
    ## gco = readRDS("~/DB/GENCODE/hg19//gencode.composite.collapsed.rds")
    gt.ge.fn = paste0(opt$outdir, "/gt.ge.rds")
    if (file.exists(gt.ge.fn) & !opt$overwrite){
        gt.ge = readRDS(gt.ge.fn)
    } else {
        gt.ge = track.gencode(gencode = opt$gencode, cached = FALSE)
        saveRDS(gt.ge, gt.ge.fn)
    }
    ## gt.ge = readRDS(paste0(opt$libdir,"/db/gt.ge.rds"))

    message("Loading gene lists")
    cgc = fready(paste0(opt$libdir,"/db/CGC.txt"))
    tsg = readRDS(paste0(opt$libdir, "/db/tsg.rds"))
    onc = readRDS(paste0(opt$libdir, "/db/onc.rds"))
    ddr = readRDS(paste0(opt$libdir, "/db/prl2015.rds"))
    
    ## cgc$onc = cgc$Significant_across_any_tumor_type_lineage_if_gain_ == 'Yes'
    ## cgc$tsg = cgc$Significant_across_any_tumor_type_lineage_if_lost_ == 'Yes'
    ## pge.cgc = pge %Q% (gene_name %in% cgc[onc == TRUE | tsg == TRUE, Gene_Symbol])
    ## cgc_tsg = pge %Q% (gene_name %in% cgc[tsg == TRUE, Gene_Symbol])
    ## cgc_onc = pge %Q% (gene_name %in% cgc[onc == TRUE, Gene_Symbol])
    ## cgc_fus = pge %Q% (gene_name %in% cgc[nchar(Translocation_Partner) > 0, Gene_Symbol])

    message("Returning Purity, Ploidy, and Coverage Variance \nAlso returning bam qc stats")
    jabba = readRDS(opt$jabba_rds)
    gg_fn = paste0(opt$outdir, "/complex.rds")
    if (file.exists(opt$complex) & file.size(opt$complex)>0){
        file.copy(opt$complex, paste0(opt$outdir, "/complex.rds"))
        gg = readRDS(opt$complex)
    } else {
        if (!(file.exists(gg_fn) & file.size(gg_fn) > 0) | opt$overwrite){
            gg = events(gG(jabba = jabba))
            # TODO: there is another instance in which we override this file later. Seems redundant
            saveRDS(gg, gg_fn)
        } else {
            gg = readRDS(gg_fn)
        }
    }

    cvgt_fn = paste0(opt$outdir, "/coverage.gtrack.rds")
    if (!file.exists(cvgt_fn) | opt$overwrite){
        cvgt = covcbs(opt$cbs_cov_rds, purity = jabba$purity, ploidy = jabba$ploidy, rebin = 5e3)
        saveRDS(cvgt, cvgt_fn)
    } else {
        cvgt = readRDS(cvgt_fn)
    }

    if (!file.exists(paste0(opt$outdir, "/oncotable.rds")) | opt$overwrite){
        ## container for all variants and events
        out = data.table()

        ## eli = rtracklayer::import.bed(paste0(opt$libdir, "/db/eligible_wgs.hg19.bed"))
        ## eli = gr.fix(eli, bcf)
        ## bcf = parsesnpeff(opt$annotated_vcf)
        ## bcf = bcf %&% eli

        ## ##########################
        ## TMB
        ## ##########################
        message("Calculating TMB")
        mut.density.fn = paste0(opt$outdir, '/mut.density.rds')
        if (check_file(mut.density.fn, overwrite = opt$overwrite)){
            mut.density = readRDS(mut.density.fn)
        } else {
            bcf = khtools::parsesnpeff(opt$annotated_vcf, filterpass = TRUE, coding_alt_only = FALSE, geno = "GT", altpipe = FALSE)
            ## bcf = grok_vcf(opt$annotated_vcf, label = opt$pair, long = TRUE, snpeff.ontology = paste0(opt$libdir, "/db/snpeff_ontology.rds"))
            ## bcf = bcf %Q% (FILTER=="PASS")
            genome.size = sum(seqlengths(bcf), na.rm = T)/1e6
            nmut = data.table(as.character(seqnames(bcf)), start(bcf), end(bcf), bcf$REF, bcf$ALT) %>% unique %>% nrow
            mut.density = data.table(
                id = opt$pair, value = c(nmut, nmut/genome.size), type = c('count', 'density'),
                track = 'tmb', source = 'annotated_bcf')
            saveRDS(mut.density, mut.density.fn)
        }
        out = rbind(out, mut.density, fill = TRUE, use.names = TRUE)

        ## TODO add overall TMB background distribution

        ## ##########################
        ## TODO: come back to SNV signatures later
        ## SNV PATTERNS
        ## ##########################
        if (file.exists(opt$deconstruct_sigs) & file.size(opt$deconstruct_sigs)>0){
            message("Returning SNV Signature stats from deconstruct_sigs pipeline")
            sig = readRDS(opt$deconstruct_sigs)
            weights = data.table::melt(as.data.table(sig$weights))
            weights = weights[order(value,decreasing=TRUE)]
            hits = as.character(weights[value >= 0.2, variable])
            weights[ ,label := gsub("Signature.","",variable)]
            ## TO BE CONTINUED
        }

        ## sigplot = ggplot(weights, aes(x = reorder(weights$label, as.numeric(weights$label)), y = value)) +
        ##     geom_col(aes(fill = variable)) +
        ##     scale_fill_manual(values = all.pal) +
        ##     theme_classic() +
        ##     theme(axis.text.x = element_text(angle = 0,vjust = 0.5,size = 18)) +
        ##     labs(title = "COSMIC Version 2 SNV Signature Contributions",
        ##          subtitle = paste0("Top hits in case ",opt$pair, ": ",
        ##                            paste0(hits,collapse = ", ")),
        ##          caption = "Generated via deconstructSigs pipeline\nReference: https://cancer.sanger.ac.uk/cosmic/signatures",
        ##          x = "",
        ##          y = "Proportional Signature Contribution",
        ##          fill = "COSMIC SNV Signature")
        ## ppng(print(sigplot),height = 5, width = 15,
        ##      filename = paste0(opt$outdir,"/SNV_Signatures.png"))

        ## message("Returning mutational burdens within TCGA context")
        ## counts = readRDS(paste0(opt$libdir,"/db/tumor_mutational_burdens.rds"))
        ## these.counts = data.table(pair = opt$pair,
        ##                           snv_count = nmut,
        ##                           sv_count = nrow(gr2dt(grok_vcf(opt$svaba_somatic_vcf))),
        ##                           tumor_type = "This_Sample")
        ## counts = rbind(counts, these.counts)
        ## counts[ ,snv_density := snv_count * 1e6/2897310462]
        ## counts[ ,sv_density := sv_count * 1e6/2897310462]
        ## cm = melt(counts, id.vars = c("pair","tumor_type"))

        ## g.g = ggplot(cm[variable == "snv_density" & !is.na(tumor_type)],
        ##             aes(x = reorder(tumor_type, value, median), value)) +
        ##     geom_sina(aes(color = tumor_type),size = 1) +
        ##     geom_point(data=cm[tumor_type == "This_Sample" & variable == "snv_density", ], colour="green", size=5) +
        ##     ##  geom_sina(data = cm[tumor_type == "This_Sample", ],color = "blue3", size = 20) +
        ##     scale_y_log10() +
        ##     scale_color_manual(values = all.pal, name = "Tumor Type") +
        ##     theme_classic() +
        ##     theme(plot.title = element_text(size = 25),
        ##           plot.subtitle = element_text(size = 20),
        ##           axis.title.x = element_text(size = 20),
        ##           axis.title.y = element_text(size = 20),
        ##           axis.text.x = element_text(size = 20,angle = 45,vjust = 0.5),
        ##           legend.title = element_text(size = 20),
        ##           legend.text = element_text(size = 15)) +
        ##     labs(title = "Tumor Mutational Burden as Compared to TCGA + IPM Whole Genomes",
        ##          subtitle = paste0("This sample: ",counts[tumor_type == "This_Sample", snv_density], " mutations per megabase"),
        ##          x = "TCGA Tumor Type",
        ##          y = "SNV + Indel Density, per Megabase")
        ## ppng(print(gg),width = 20,filename = paste0(opt$outdir,"/",opt$pair,"TMB_Context.png"))
        ## message("Returning Burden Weighted Signatures")
        ## snv.total = cm[tumor_type == "This_Sample" & variable == "snv_count", value]
        ## weights$burden = weights$value * snv.total
        ## weights2 = weights[burden != 0]
        ## sigplot = ggplot(weights2, aes(x = reorder(weights2$label, as.numeric(weights2$label)), y = burden)) +  geom_col(aes(fill = variable)) +
        ##     scale_fill_manual(values = all.pal) +
        ##     theme_classic() +
        ##     theme(axis.text.x = element_text(angle = 0,vjust = 0.5,size = 18)) +
        ##     labs(title = "SNV Signature Burden",
        ##          subtitle = paste0("Top hits in case ",opt$pair, ": ",
        ##                            paste0(hits,collapse = ", ")),
        ##          caption = "Generated via deconstructSigs pipeline\nReference: https://cancer.sanger.ac.uk/cosmic/signatures",
        ##          x = "",
        ##          y = "Proportional Signature Contribution",
        ##          fill = "COSMIC SNV Signature") +
        ##     geom_text(aes(y = burden + 40, label = floor(burden)), size = 6)

        ## ppng(print(sigplot),height = 5, width = 15, units = "in", res = 1,
        ##      filename = paste0(opt$outdir,"/SNV_Sig_burden.png"))



        ## #####################
        ## SV EVENT PATTERNS
        ## #####################
        message("Returning Windows of complex variants")
        ## SV events
        ev_fn = paste0(opt$outdir, "/events.rds")
        if (file.exists(ev_fn) & file.size(ev_fn) >0 & !opt$overwrite){
            ev = readRDS(ev_fn)
        } else {
            ev = gg$meta$event
            if (nrow(ev)>0){
                ev[, ev.ix := seq_len(.N), by = type]
                ev[, ev.id := paste(type, ev.ix, sep = "_")]
                ev.fp = ev[, gr.stripstrand(parse.gr(footprint)) + 1e4]
                ev.fp$ev.id = ev$ev.id[ev.fp$grl.ix]
                ev.fp = split(ev.fp, ev.fp$ev.id)

                ## ev[, ev.png := paste0(opt$outdir,"/sv_events/",opt$pair,"_",ev.id,".png")]
                ## ev[, ev.js.range := str_replace(gsub("+|-", "", footprint), ",", "%20")]
                ev[, ev.js.range := {
                    gr = gr.reduce(gr.stripstrand(parse.gr(footprint)) + 1e4)
                    paste(gr.string(gr), collapse = "%20|%20")
                }, by = ev.id]
                ev[, ev.js := paste0(
                         opt$server,
                         "index.html?file=",
                         opt$pair,
                         ".json&location=",
                         ev.js.range,
                         "&view=")]
                saveRDS(ev, ev_fn)

            }
        }
        if (nrow(ev) > 0){
            ## TODO: also add burden
            ## only save the count in oncotable
            sv = ev[, .(value = .N), by = type][, id := opt$pair][
              , track := ifelse(type %in% c('del', 'dup', 'invdup', 'tra', 'inv'),
                                'simple sv', 'complex sv')][
              , source := 'complex']
            out = rbind(out, sv, fill = TRUE, use.names = TRUE)
        }

        ## generate gTracks
        ## for (i in seq_len(nrow(ev))){
        ##     if (!file.exists(ev[i, ev.png]) | overwrite){
        ##         ppng({
        ##             fp = ev.fp %Q% (ev.id==ev[i, ev.id])
        ##             fpw = median(width(fp))
        ##             plot(c(cvgt, gg$gt), fp + fpw/2)
        ##         }, ev[i, ev.png],
        ##         width = 1600, height = 900)
        ##     }
        ## }

        ## #####################
        ## SNV & INDELs
        ## #####################
        ## message("Filtering the input VCF by quality and annotation")
        ## out.name = paste0(opt$pair, ".filtered.vcf.gz")
        ## ## tmp.path = paste0(tempdir(), "/", out.name)
        ## tmp.path = paste0(opt$outdir, "/", out.name)
        ## catcmd = if (grepl("(.gz)$", opt$annotated_vcf)) "zcat" else "cat"
        ## onepline = paste0(opt$libdir,"/db/vcfEffOnePerLine.pl")
        ## ## dt.patho = gr2dt(readRDS(paste0(opt$libdir,"/db/pathogenic.rds")))
        ## ## filt = paste0("java -Xmx1g -Xms1g -jar ",opt$libdir,"/SnpEff/source/snpEff/SnpSift.jar filter \"( ANN =~ 'missense|splice|stop_gained|frame' )\"")
        ## sift = "java -Xmx20m -Xms20m -XX:ParallelGCThreads=1 -jar /gpfs/commons/groups/imielinski_lab/git/mskilab/flows/modules/SnpEff/source/snpEff/SnpSift.jar filter \"( ANN =~ 'chromosome_number_variation|exon_loss_variant|rare_amino_acid|stop_lost|transcript_ablation|coding_sequence|regulatory_region_ablation|TFBS|exon_loss|truncation|start_lost|missense|splice|stop_gained|frame' )\""
        ## cmd = sprintf(
        ##     paste(catcmd, "%s | %s | %s | bcftools view -i 'FILTER==\"PASS\"' | bgzip -c > %s"),
        ##     opt$annotated_vcf, onepline, sift, tmp.path)

        ## if (!file.exists(tmp.path)){
        ##     system(cmd)
        ## }
        
        message("Loading coding mutations only.")
        ## SOMATIC first
        ## mutations with some effect
        ## WARNING: not working with KH update Thursday, Mar 25, 2021 09:44:55 AM
        ## som = khtools::parsesnpeff(tmp.path, filterpass = TRUE, coding_alt_only = TRUE, geno = "GT", altpipe = FALSE)
        som.dt.fn = paste0(opt$outdir, "/somatic.mutations.rds")
        if (file.exists(som.dt.fn) & file.size(som.dt.fn) > 0 & !opt$overwrite){
            som.dt = readRDS(som.dt.fn)
        } else {
            som = bcf %Q% (grepl('chromosome_number_variation|exon_loss_variant|rare_amino_acid|stop_lost|transcript_ablation|coding_sequence|regulatory_region_ablation|TFBS|exon_loss|truncation|start_lost|missense|splice|stop_gained|frame', annotation))

            ## NOTE: mutations can be in multiple lines as it affects different transcripts!!

            ## annotate copy numbers
            som.cn = readRDS(opt$somatic_snv_cn)
            som = som %$% dt2gr(som.cn[!is.na(cn)])[, c("cn", "est_cn_ll")]
            som$variant.g = paste0(seqnames(som), ':', start(som), '-', end(som), ' ', som$REF, '>', som$ALT)
            ## deduplicate!!!!!!!!!!
            som = som %Q% (!duplicated(variant.g))
            ## variant.g = paste(gr.string(som), paste0(som$REF, ">", som$ALT))
            som$type = "short"
            som$track = "variants"
            som$source = "annotated_bcf"
            som.dt = gr2dt(som)
            ## som.dt[, factor(impact, levels = rev(c("LOW", "MODERATE", "HIGH")))]
            ## som.dt = som.dt[order(impact)]
            saveRDS(som.dt, som.dt.fn)
        }
        vars = som.dt[, .(id = opt$pair, gene, vartype, variant.g, variant.p, distance, annotation, type = "short", track = 'variants', source = 'annotated_bcf', vcn = est_cn_ll)] %>% unique
        out = rbind(out, vars, fill = TRUE, use.names = TRUE)

        ## annotate cool.ge with overlapping functional somatic short variants
        ## cool.ge %*%

        ## GERMLINE mutations
        message("Loading germline pathogenic mutations")
        ger.dt.fn = paste0(opt$outdir, "/germline.mutations.rds")
        if (file.exists(ger.dt.fn) & file.size(ger.dt.fn) > 0 & !opt$overwrite){
            ger.dt = readRDS(ger.dt.fn)
        } else {
            ger = khtools::parsesnpeff(opt$annotated_vcf_germline, filterpass = TRUE, coding_alt_only = TRUE, geno = "GT", altpipe = FALSE)
            ## TOO many, only select for really impactful ones
            ger = ger %Q% (impact=="HIGH")
            ger.cn = readRDS(opt$germline_snv_cn)
            ger.cn = dt2gr(ger.cn[, .(seqnames, start, end, REF, ALT, cn, est_cn_ll)])
            ger = ger %$% ger.cn[, c("cn", "est_cn_ll")]
            ger$variant.g = paste0(seqnames(ger), ':', start(ger), '-', end(ger), ' ', ger$REF, '>', ger$ALT)
            ger = ger %Q% (!duplicated(variant.g))
            ## variant.g = paste(gr.string(ger), paste0(ger$REF, ">", ger$ALT))
            ger$type = "short_germline"
            ger$track = "variants"
            ger$source = "annotated_vcf_germline"
            ger.dt = gr2dt(ger)
            ## ger.dt[, factor(impact, levels = rev(c("LOW", "MODERATE", "HIGH")))]
            ## ger.dt = ger.dt[order(impact)]
            saveRDS(ger.dt, ger.dt.fn)
        }
        gvars = ger.dt[, .(id = opt$pair, gene, vartype, variant.g, variant.p, distance, annotation, type, track = 'variants', source, vcn = est_cn_ll)] %>% unique
        out = rbind(out, gvars, fill = TRUE, use.names = TRUE)

        ## ##################
        ## FUNCTIONAL CNA EVENTS
        ## ##################
        message("Calling CNA events involving any cancer gene")
        pl = gg$nodes$dt[seqnames %in% names(sl)][, sum(cn * width, na.rm = TRUE)/sum(width, na.rm = TRUE)]
        gg$set(ploidy = pl)
        gg$set(purity = jabba$purity)
        ## FIXME: we already save gg to this file earlier so let's get rid of one of these.
        saveRDS(gg, paste0(opt$outdir, "/complex.rds"))

      # get the ncn data from jabba
      kag = readRDS(gsub("jabba.simple.rds", "karyograph.rds", opt$jabba_rds))
      nseg = NULL
      if ('ncn' %in% names(mcols(kag$segstats))){
          nseg = kag$segstats[,c('ncn')]
      }

      scna.fn = paste0(opt$outdir, '/scna.rds')
      if (check_file(scna.fn, overwrite = opt$overwrite)){
          scna = readRDS(scna.fn)
      } else {
          scna = skitools::get_gene_ampdels_from_jabba(jabba, amp.thresh = opt$amp_thresh,
                                         del.thresh = opt$del_thresh, pge = pge, nseg = nseg)

            if (nrow(scna))
            {
                scna$col = case_when(
                    scna$type == "homdel" ~ "#2166ac",
                    scna$type == "hetdel" ~ "#92c5de",
                    scna$type == "amp" ~ "#b2182b",
                    TRUE ~ NA_character_
                )
                saveRDS(scna, scna.fn)
            }
      }
      if (scna[, .N] > 0){
            cool.scna = scna[gene_name %in% c(onc, tsg, ddr)]
            ## cool.scna = rbind(
            ##     scna[(gene_name %in% union(tsg, ddr)) & grepl("del", type)],
            ##     scna[(gene_name %in% onc) & grepl("amp", type)]
            ## )
            if (nrow(cool.scna))
            {
                cool.scna[, track := 'variants'][, source := 'jabba_rds'][, vartype := 'scna']
                # FIXME: why only save the oncogenes and TSGs to the oncotable? let's save all gene CNVs and wherever we need, we will subset according to the type of gene
                out = rbind(out,
                            cool.scna[, .(id = opt$pair, value = type, type, track, gene = gene_name)],
                            fill = TRUE, use.names = TRUE)
            }
        }

        ## ##################
        ## FUSION GENES
        ## ##################
        message("Returning Fusion Calls, only multi gene, in-frame fusions")
        ## gg = gg ## gG(jabba = opt$jabba_rds)
        fu.fn = paste0(opt$outdir, "/fusions.rds")
        if (check_file(fu.fn, overwrite = opt$overwrite)){
            if (file.exists(opt$fusions) & file.size(opt$fusions)>0){
                fu = readRDS(opt$fusions)
            } else {
                if (is.na(gff)){
                    gff = skidb::read_gencode(fn = opt$gencode)
                }
                fu = fusions(gg, gff)
            }

            fu$set(og.wid = fu$dt$walk.id)
            fus <- fu$meta
            fix = fus[silent == FALSE, ][!duplicated(genes), walk.id]
            if (length(fix)>0){
                fu = fu[fix]
                saveRDS(fu, fu.fn)
                ## TODO: Let's pull out the DNA junctions supporting reads
                ## TODO: Let's pull out the RNA fusion transcripts
                ## if (file.exists(opt$star_bam)){}
            }
        } else {
            fu = readRDS(fu.fn)
        }
        fus = fu$meta
        ## save into the oncotable
        fus[, vartype := ifelse(in.frame == TRUE, 'fusion', 'outframe_fusion')] # annotate out of frame fusions
        fus = fus[, .(gene = strsplit(genes, ',') %>% unlist,
                      vartype = rep(vartype, sapply(strsplit(genes, ','), length)))][
          , id := opt$pair][, track := 'variants'][, type := vartype][, source := 'fusions']
        out = rbind(out, fus, fill = TRUE, use.names = TRUE)

        ## ###################
        ## OVEREXPRESSION of ONCs
        ## ###################
        message("Loading RNA quantities matrix")
        ## TODO

        message("Save oncotable content")
        ## annotate gene roles
        out[, role := case_when(
                  out$gene %in% onc ~ "ONC",
                  out$gene %in% tsg & out$gene %in% ddr ~ "TSG;DDR",
                  out$gene %in% tsg & !out$gene %in% ddr ~ "TSG",
                  !out$gene %in% tsg & out$gene %in% ddr ~ "DDR",
                  TRUE ~ "none"
              )]

        ## annotate which rows are LOF of TSG/DDR
        ## annotate which rows are GOF of ONC
        out[grepl("DDR|TSG", role) & type=="homdel", lof := TRUE]
        out[grepl("DDR|TSG", role) & grepl("short", type) &
            grepl("stop|frameshift|missense|disruptive|trunc", annotation), lof := TRUE]
        out[grepl("DDR|TSG", role) & type=="hetdel", lof := TRUE]
        out[grepl("ONC", role) & type=="amp", gof := TRUE]
        out[grepl("ONC", role) & grepl("short", type) & grepl("missense", annotation), gof := TRUE]
        saveRDS(out, paste0(opt$outdir, "/oncotable.rds"))

        ## prepare the gene-centric dosage table
        gs = out[, .(role = role[1],
                     nlof = sum(lof==TRUE, na.rm = T),
                     ngof = sum(gof==TRUE, na.rm = T),
                     homdel = "homdel" %in% type), by = gene]
        gs = rbind(gs[grepl("TSG|DDR", role) & (homdel==TRUE | nlof>1)],
                   gs[grepl("ONC", role) & (ngof>0)])
        if (nrow(gs)>0){
            ## cool.ge = pge[, "gene_name"] %Q% (gene_name %in% gs$gene)
            ## cool.ge$role = case_when(
            ##     cool.ge$gene_name %in% onc ~ "ONC",
            ##     cool.ge$gene_name %in% tsg & cool.ge$gene_name %in% ddr ~ "TSG;DDR",
            ##     cool.ge$gene_name %in% tsg & !cool.ge$gene_name %in% ddr ~ "TSG",
            ##     !cool.ge$gene_name %in% tsg & cool.ge$gene_name %in% ddr ~ "DDR",
            ##     TRUE ~ "none"
            ## )
            saveRDS(gs, paste0(opt$outdir, "cool.ge.rds"))
        }
    }
    message("Everything ready, now knit")

    ## if(length(fus)>0){
    ##     lapply(unique(fus$dt[numgenes > 1, genes]),function(x){
    ##         ppng(plot(c(cov,gg$gt,gt.ge),
    ##                   fus[name == x]$footprint + gr2dt(fus[name == x]$footprint)$width,
    ##                   legend.params = list(plot = FALSE)),
    ##              filename = paste0(opt$outdir,"/fusions/",opt$pair,"-",x,".png"))})
    ##     saveRDS(fus, paste0(opt$outdir,"/fusion.table.rds"))
    ## }

    ## ######################
    ## match up with CIVIC
    ## ######################
    ## message("Accessing CIVIC database for detected mutations, data updated as of sept 1, 2019")
    ## ## civ = fread("/gpfs/commons/groups/imielinski_lab/git/mskilab/flows/modules/CaseReport/01-Sep-2019-VariantSummaries.tsv", fill = TRUE)
    ## civ = fread("~/DB/civic/nightly-VariantSummaries.tsv", sep = "\t", fill = TRUE)[, 1:27]

    ## civ.vartypes = civ[, unique(variant_types)]
    ## ## civ = fread("~/DB/civic/01-Mar-2021-VariantSummaries.tsv", sep = "\t")

    ## message("An EXTREMELY hypermutated sample may overwhelm the CIVIC API at this step and fail to return variants")
    ## hits = mt[mt$gene %in% civ$gene, ]
    ## ## civic = "https://civicdb.org/api/genes/"
    ## civic = "https://civicdb.org/api/variants"

    ## request = GET(url = paste0(civic,paste0(hits$gene, collapse = ",")),
    ##               query = list(identifier_type = "entrez_symbol"))

    ## res =  data.frame(fromJSON(content(request, as = "text", encoding = "UTF-8"),
    ##                            flatten = TRUE)$records)
########OOF we need protein variants here - gotta go get them from Alterations
    ## if(nrow(res) >0){
    ##     vars = lapply(1:nrow(res), function(n){
    ##         var = res[n, "variants"][[1]]
    ##         var$gene = res$name[n]
    ##         return(as.data.table(var))})
    ##     vars = do.call(rbind, vars)
    ##     vars[ ,variantp := paste0(vars$gene, ":", name)]
    ##     vars$variantp = gsub("OVEREXPRESSION", "amp", vars$variantp)
    ##     vars$variantp = gsub("AMP.*", "amp", vars$variantp)
    ##     hits$variantp = gsub("hetdel|homdel","UNDEREXPRESSION", hits$variantp)
    ##                                         #vars$variantp = gsub("UNDEREXPRESSION", "homdel", vars$variantp)
    ##     vars$variantp = gsub("LOSS.*", "UNDEREXPRESSION", vars$variantp)
    ##     vars$variantp = gsub("DEL.*", "homdel", vars$variantp)
    ##     hits2 = hits[hits$variantp %in% vars$variantp, ]
    ##     hits3= merge(hits2, vars[ ,c("variantp","id","evidence_items.accepted_count")], by = "variantp")
    ##     setnames(hits3, old = "id", new = "CIVIC.id")
    ## } else {
    ##     hits3 = data.table()
    ## }

    ## if(nrow(hits3) >0){
    ##     evidence = lapply(1:nrow(hits3),function(n){
    ##         civic2 = "https://civicdb.org/api/variants/"
    ##         request = GET(url = paste0(civic2,paste0(hits3$CIVIC.id[n], collapse = ",")))
    ##         civic.var =  as.data.table(data.frame(fromJSON(content(request,
    ##                                                                as = "text",
    ##                                                                encoding = "UTF-8"),
    ##                                                        flatten = TRUE)$evidence_items))
    ##         civic.var$gene = hits3$gene[n]
    ##         return(civic.var[ , .(gene, id, name, description, drugs, rating, evidence_level, evidence_type, clinical_significance, evidence_direction, variant_origin, phenotypes, disease.name, source.citation, source.source_url, source.publication_date.year)])})
    ##     evidence = do.call(rbind, evidence)
    ##     saveRDS(evidence, paste0(opt$outdir,"/","CIVIC.report.rds"))
    ##     saveRDS(hits3, paste0(opt$outdir,"/","CIVIC.hits.rds"))
    ## }

    ## ######################
    ## match up with OncoKB
    ## ######################
    oncokb.token = opt$oncokb_token
    if (file.exists(oncokb.token)){
        message('Querying OncoKB to annotate genomic alterations')
        oncokb = get_oncokb_response(som.dt, oncokb.token = oncokb.token)
        oncokb_annotations = get_oncokb_annotations(oncokb)
        oncokb_annotations[, key_for_merging := .I]
        som.dt[, key_for_merging := .I]
        som.dt = merge(som.dt, oncokb_annotations, by = 'key_for_merging')
        som.dt$key_for_merging = NULL
        som.dt$onco_kb_entry_url = get_oncokb_gene_entry_url(oncokb)
    } else {
        message('No OncoKB token was provided so skipping OncoKB annotations')
    }

    ## build one query
    ## oncokb = "https://www.oncokb.org/api/v1/annotate"
    ## oncokb.token = readLines("~/keys/oncokb.token")

    ## wtf = som %Q% (gene=="HCFC1") %>% head(1)
    ## wtf = som[1]
    ## all.res = lapply(
    ##     1:nrow(som.dt), 
    ##     function(i){
    ##         res = GET(
    ##             url = paste0(
    ##                 oncokb,
    ##                 "/mutations/byGenomicChange?",
    ##                 "genomicLocation=", som.dt[i, asc(seqnames)], "%2C",
    ##                 som.dt[i, start], "%2C",
    ##                 som.dt[i, end], "%2C",
    ##                 som.dt[i, REF], "%2C",
    ##                 som.dt[i, ALT],
    ##                 "&referenceGenome=GRCh37"),
    ##             add_headers(authorization = paste("Bearer", oncokb.token)),
    ##             add_headers(accept = "application/json")
    ##         )
    ##     })

    ## som.dt$oncokb.status = sapply(all.res, function(res) res$status_code)
    ## success.ix = which(between(som.dt$oncokb.status, 200, 299))
    ## som.dt$oncokb.variantExist = sapply(all.res, function(res) res$)

    ## message("Coverage data QC")
    ## cbs = readRDS(opt$cbs_cov_rds)
    ## variance = dlrs(log(cbs$ratio))
    ## qc = data.table(purity = jabba$purity, ploidy = jabba$ploidy, variance = variance)
    ## qcm = data.table::melt(qc)
    ## qcm[ ,remainder := 1 - value]

    ## ## a pie chart of one number: purity
    ## pur = ggplot(qcm[variable == "purity"], aes(x = 2, y = value, fill = variable)) +
    ##     geom_bar(stat = "identity", color = "white") +
    ##     coord_polar(theta = "y", start = 0)+
    ##     scale_fill_manual(values = pal1) +
    ##     theme_void()+
    ##     theme(legend.position = "none") +
    ##     xlim(0.5, 2.5) +
    ##     ylim(0,1) +
    ##     labs(title = "Sequenza estimated sample purity",
    ##          subtitle = paste0("Sample ",opt$pair),
    ##          caption = "Note: 30% purity is the rough cutoff for \ndecent quality SV/JaBbA modeling") +
    ##     annotate(geom = 'text', x = 0.5, y = 0, label = qc$purity,size = 15)

    ## ppng(print(pur),height = 300, width = 350,
    ##      filename = paste0(opt$outdir,"/purity_donut.png"))

    ## ## TODO: add way more details to the quality control info
    ## ## 1) PP fit contour
    ## ## 2) variant supporting reads of SNVs (IGV)
    ## ## 3) variant supporting reads of SVs (samplot)
    ## write.csv(qc, paste0(opt$outdir,"/",opt$pair,'_quality_metrics.csv'))


}


####################################
## KNITING THE PAGE
####################################
message("now we will try to knit")
rmarkdown::render(
    ## input = "~/git/CaseReport/report.rmd",
    input = paste0(opt$libdir, "/report.rmd"),
    ## TODO: temporary sandbox!!!!!!!!!
    ## input = "/gpfs/commons/groups/imielinski_lab/projects/IPM/Flow.test/casereport/PM1087-Z1/report.rmd",
    output_format = "html_document",
    output_file = paste0(opt$outdir, "/", opt$pair,".report.html"),
    ## knit_root_dir = normalizePath(opt$outdir),
    knit_root_dir = opt$outdir,
    ## "/gpfs/commons/groups/imielinski_lab/projects/IPM/Flow.test/casereport/PM1087-Z1",
    params = list(set_title = paste0(opt$pair),
                  pair = opt$pair,
                  ## pair_short = opt$pair_short,
                  outdir = normalizePath(opt$outdir),
                  jabba = normalizePath(opt$jabba_rds),
                  bam = opt$tumor_bam,
                  tumor_type = opt$tumor_type,
                  server = opt$server),
    quiet = FALSE)

