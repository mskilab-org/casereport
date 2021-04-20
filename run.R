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
        make_option(c("--somatic_snv_cn"), type = "character", help = "MLE somatic SNV CN estimate from snv_multiplicity2 module"),
        make_option(c("--germline_snv_cn"), type = "character", help = "MLE germline SNV CN estimate from snv_multiplicity2 module")
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
    sl = hg_seqlengths(chr = opt$chr)
    saveRDS(sl, paste0(opt$outdir, "/sl.rds"))
    hg = si2gr(sl) %>% gr.stripstrand %>% unname
    saveRDS(hg, paste0(opt$outdir, "/hg.rds"))
    ## bands.gt = gTrack::karyogram(height = 2.5)
    ## bands = grl.unlist(bands.gt@data[[1]])
    ## bands = bands[, "band"]
    ## values(bands)$arm = gsub("^(.*)([pq])(.*)$", "\\1\\2", values(bands)$band)
    ## arm = gr.reduce(bands, by = "arm")

    message("Loading Gencode genes")
    ## gff = readRDS(gzcon(url('http://mskilab.com/gGnome/hg19/gencode.v19.annotation.gtf.gr.rds')))
    gff = skidb::read_gencode(fn = opt$gencode)
    pge = gff %Q% (type=="gene") %Q% (gene_type=="protein_coding")
    saveRDS(pge, paste0(opt$outdir, "/pge.rds"))
    ## ge = readRDS(paste0(opt$libdir,"/db/ge.rds"))
    ## pge = readRDS(paste0(opt$libdir,"/db/pge.rds"))
    ## ug = readRDS(paste0(opt$libdir,"/db/ug.rds"))
    ## gco2 = split(gff, gff$gene_name)
    ## gco = readRDS("~/DB/GENCODE/hg19//gencode.composite.collapsed.rds")
    gt.ge = track.gencode(gencode = opt$gencode, cached = FALSE)
    saveRDS(gt.ge, paste0(opt$outdir, "/gt.ge.rds"))
    ## gt.ge = readRDS(paste0(opt$libdir,"/db/gt.ge.rds"))

    message("Loading gene lists")
    cgc = fready(paste0(opt$libdir,"/db/CGC.txt"))
    tsg = readRDS(paste0(opt$libdir, "db/tsg.rds"))
    onc = readRDS(paste0(opt$libdir, "db/onc.rds"))
    ddr = readRDS(paste0(opt$libdir, "db/prl2015.rds"))
    ## cgc$onc = cgc$Significant_across_any_tumor_type_lineage_if_gain_ == 'Yes'
    ## cgc$tsg = cgc$Significant_across_any_tumor_type_lineage_if_lost_ == 'Yes'
    ## pge.cgc = pge %Q% (gene_name %in% cgc[onc == TRUE | tsg == TRUE, Gene_Symbol])
    ## cgc_tsg = pge %Q% (gene_name %in% cgc[tsg == TRUE, Gene_Symbol])
    ## cgc_onc = pge %Q% (gene_name %in% cgc[onc == TRUE, Gene_Symbol])
    ## cgc_fus = pge %Q% (gene_name %in% cgc[nchar(Translocation_Partner) > 0, Gene_Symbol])

    message("Returning Purity, Ploidy, and Coverage Variance \nAlso returning bam qc stats")
    jabba = readRDS(opt$jabba_rds)
    gg = readRDS(opt$complex)
    file.symlink(opt$complex, paste0(opt$outdir, "/complex.rds"))

    if (!file.exists(paste0(opt$outdir, "/coverage.gtrack.rds"))){
        cvgt = covcbs(opt$cbs_cov_rds, purity = jabba$purity, ploidy = jabba$ploidy, rebin = 5e3)
        saveRDS(cvgt, paste0(opt$outdir, "/coverage.gtrack.rds"))
    } else {
        cvgt = readRDS(paste0(opt$outdir, "/coverage.gtrack.rds"))
    }
    

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
    bcf = khtools::parsesnpeff(opt$annotated_vcf, filterpass = TRUE, coding_alt_only = FALSE, geno = "GT", altpipe = FALSE)
    ## bcf = grok_vcf(opt$annotated_vcf, label = opt$pair, long = TRUE, snpeff.ontology = paste0(opt$libdir, "/db/snpeff_ontology.rds"))
    ## bcf = bcf %Q% (FILTER=="PASS")
    genome.size = sum(seqlengths(bcf), na.rm = T)/1e6
    nmut = data.table(as.character(seqnames(bcf)), start(bcf), end(bcf), bcf$REF, bcf$ALT) %>% unique %>% nrow
    mut.density = data.table(id = opt$pair, value = c(nmut, nmut/genome.size), type = c('count', 'density'),
                             track = 'tmb', source = 'annotated_bcf')
    out = rbind(out, mut.density, fill = TRUE, use.names = TRUE)

    ## ##########################
    ## TODO: come back to SNV signatures later
    ## SNV PATTERNS
    ## ##########################
    message("Returning SNV Signature stats from deconstruct_sigs pipeline")
    sig = readRDS(opt$deconstruct_sigs)
    weights = data.table::melt(as.data.table(sig$weights))
    weights = weights[order(value,decreasing=TRUE)]
    hits = as.character(weights[value >= 0.2, variable])
    weights[ ,label := gsub("Signature.","",variable)]

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



#######################
    ## SV EVENT PATTERNS
#######################
    message("Returning Windows of complex variants")
    ## SV events
    ev = gg$meta$event
    ev[, ev.ix := seq_len(.N), by = type]
    ev[, ev.id := paste(type, ev.ix, sep = "_")]
    ev.fp = ev[, parse.gr(footprint)]
    ev.fp$ev.id = ev$ev.id[ev.fp$grl.ix]

    ## ev[, ev.png := paste0(opt$outdir,"/sv_events/",opt$pair,"_",ev.id,".png")]
    ev[, ev.js.range := str_replace(footprint, ",", "%20")]
    ev[, ev.js := paste0(
             opt$server,
             "index.html?file=",
             opt$pair,
             ".json&location=",
             ev.js.range,
             "&view=")]
    saveRDS(ev, paste0(opt$outdir, "/events.rds"))

    ## TODO: also add burden
    ## only save the count in oncotable
    sv = ev[, .(value = .N), by = type][, id := opt$pair][
      , track := ifelse(type %in% c('del', 'dup', 'invdup', 'tra', 'inv'),
                        'simple sv', 'complex sv')][
      , source := 'complex']
    out = rbind(out, sv, fill = TRUE, use.names = TRUE)

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
    som = bcf %Q% (grepl('chromosome_number_variation|exon_loss_variant|rare_amino_acid|stop_lost|transcript_ablation|coding_sequence|regulatory_region_ablation|TFBS|exon_loss|truncation|start_lost|missense|splice|stop_gained|frame', annotation))

    ## NOTE: mutations can be in multiple lines as it affects different transcripts!!

    ## annotate copy numbers
    som.cn = readRDS(opt$somatic_snv_cn)
    som = som %$% dt2gr(som.cn[!is.na(cn)])[, c("cn", "est_cn_ll")]
    som$variant.g = paste0(seqnames(som), ':', start(som), '-', end(som), ' ', som$REF, '>', som$ALT)
    ## variant.g = paste(gr.string(som), paste0(som$REF, ">", som$ALT))
    som$type = "short"
    som$track = "variants"
    som$source = "annotated_bcf"

    som.dt = gr2dt(som)
    ## som.dt[, factor(impact, levels = rev(c("LOW", "MODERATE", "HIGH")))]
    ## som.dt = som.dt[order(impact)]

    vars = som.dt[, .(id = opt$pair, gene, vartype, variant.g, variant.p, distance, annotation, type = "short", track = 'variants', source = 'annotated_bcf')] %>% unique
    out = rbind(out, vars, fill = TRUE, use.names = TRUE)

    saveRDS(som.dt, paste0(opt$outdir, "/somatic.mutations.rds"))

    ## ##################
    ## FUNCTIONAL CNA EVENTS
    ## ##################
    message("Calling CNA events involving any cancer gene")
    pl = gg$nodes$dt[seqnames %in% names(sl)][, sum(cn * width, na.rm = TRUE)/sum(width, na.rm = TRUE)]
    gg$set(ploidy = pl)

    ## TODO:  ## infer CNAs
    scna = rbind(
        gg$nodes$dt[
            cn>=opt$amp_thresh * gg$meta$ploidy, ][, type := 'amp'],
        gg$nodes$dt[
            !seqnames %in% c("X", "Y")][
            (cn == 1) | (cn < opt$del_thresh * gg$meta$ploidy), ][
          , type := 'hetdel'],
        gg$nodes$dt[
            !seqnames %in% c("X", "Y")][cn == 0, ][, type := 'homdel']
    )

    if (nrow(scna))
    {
        ## TODO: ignore X/Y chromosome for now
        scna = dt2gr(scna, seqlengths = seqlengths(gg)) %*% pge[, 'gene_name'] %>% gr2dt
        cool.scna = rbind(
            scna[(gene_name %in% union(tsg, ddr)) & grepl("del", type)],
            scna[(gene_name %in% onc) & grepl("amp", type)]
        )
        if (nrow(cool.scna))
        {
            cool.scna[, track := 'variants'][, source := 'jabba_rds'][, vartype := 'scna']
            out = rbind(out,
                        cool.scna[, .(id = opt$pair, value = cn, type, track, gene = gene_name)],
                        fill = TRUE, use.names = TRUE)
        }
    }

    ## ##################
    ## FUSION GENES
    ## ##################
    message("Returning Fusion Calls, only multi gene, in frame fusions")
    ## gg = gg ## gG(jabba = opt$jabba_rds)
    if (file.exists(opt$fusions) & file.size(opt$fusions)>0){
        fu = readRDS(opt$fusions)
    } else {
        fu = fusions(gg, gff)
    }

    fu$set(og.wid = fu$dt$walk.id)
    fus <- fu$meta
    fix = fus[silent == FALSE, ][!duplicated(genes), walk.id]
    if (length(fix)>0){
        fu = fu[fix]
        fus = fu$meta
        ## save into the oncotable
        fus[, vartype := ifelse(in.frame == TRUE, 'fusion', 'outframe_fusion')] # annotate out of frame fusions
        fus = fus[, .(gene = strsplit(genes, ',') %>% unlist,
                      vartype = rep(vartype, sapply(strsplit(genes, ','), length)))][
          , id := opt$pair][, track := 'variants'][, type := vartype][, source := 'fusions']
        out = rbind(out, fus, fill = TRUE, use.names = TRUE)
        fu.fn = paste0(opt$outdir, "/fusions.rds")
        saveRDS(fu, fu.fn)
        ## TODO: Let's pull out the RNA fusion transcripts
        ## if (file.exists(opt$star_bam)){}
    }


    message("Save oncotable content")
    ## annotate gene roles
    out[, role := case_when(
              out$gene %in% onc ~ "ONC",
              out$gene %in% tsg & out$gene %in% ddr ~ "TSG;DDR",
              out$gene %in% tsg & !out$gene %in% ddr ~ "TSG",
              !out$gene %in% tsg & out$gene %in% ddr ~ "DDR",
              TRUE ~ "none"
          )]
    saveRDS(out, paste0(opt$outdir, "/oncotable.rds"))


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
    ##                 oncokb, "/mutations/byGenomicChange?",
    ##                 "genomicLocation=", som.dt[i, asc(seqnames)], "%2C",
    ##                 som.dt[i, start], "%2C",
    ##                 som.dt[i, end], "%2C",
    ##                 som.dt[i, REF], "%2C", som.dt[i, ALT],
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
                  bam = opt$tumor_bam,
                  tumor_type = opt$tumor_type,
                  server = opt$server),
    quiet = FALSE)
