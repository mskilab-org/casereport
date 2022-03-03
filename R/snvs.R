##### Utility functions to process SNVs and indels #####

#' @name process_germline_muts
#' @title process_germline_muts
#'
#' @description
#'
#' Read and process the germline mutations in coding genes
#'
#' @param germline_coding the "annotated.germline.coding.rds" output from the Alterations module
#' @param driver_germline_mutations path to output file in which to save germline mutations of driver genes
#' @param germline_genes list of germline genes that should be included in the output (if not provided then any gene annotated as "Germline" in CGC in addition to all DDR genes will be included
#' @param ... other parameters accepted by casereport:::annotate_germline_mutations
#' @param return.type (character) one of "data.table" or "GRanges"
#' 
#' @return data.table with node.fp, amp.fp, and ev.fp for node, event, and amplicon footprints
process_germline_muts = function(germline_coding, driver_germline_mutations, name = '', germline_snpeff_snv_bcf = NULL, germline_genes = NULL, ...){
    if (!file.good(germline_coding)){
        if (file.good(germline_snpeff_snv_bcf)){
            annotate_germline_mutations(germline_snpeff_snv_bcf, name = name,
                                        germline_coding = germline_coding,
                                        ...)
        } else {
            empty_germline_dt = data.table(gene = as.character(),
                                     variant.c = as.character(),
                                     variant.p = as.numeric(),
                                     vartype = as.numeric(),
                                     CLNDN = as.character(),
                                     CLNSIG = as.character(),
                                     CLNVC = as.character(),
                                     annotation = as.character(),
                                     impact = as.character())
            fwrite(empty_empty_germline_dt, driver_germline_mutations)
            return(empty_empty_germline_dt)
        }
    }

    germline_dt = readRDS(germline_coding)

    germline_dt_dedup = germline_dt
    germline_dt_dedup[, uid := paste(gene, variant.c)]

    # in case the AA variation is different in different transcripts we will take all versions
    germline_dt_dedup[, variant.p := paste(unique(variant.p), collapse = ','), by = uid]

    # deduping
    germline_dt_dedup = germline_dt_dedup[!duplicated(uid)]
    if (is.null(germline_genes) || length(germline_genes) == 0){
        # keeping only genes annotated as "Germline" in CGC
        # in addition to all DDRs from Pearl et al.
        cgc = fread(system.file('extdata', 'cgc.tsv', package = 'casereport'))
        germline_genes = cgc[Germline == 'yes', get("Gene Symbol")]
        pearl = fread(system.file('extdata', 'ddr_pearl2015.tsv', package = 'casereport'))
        germline_genes = unique(c(germline_genes, pearl$V1))
    }
    germline_dt_dedup = germline_dt_dedup[gene %in% germline_genes]
    # order so that pathogenic mutations appear first, then germline, and lastly VUS
    germline_dt_dedup = germline_dt_dedup[order(germline_pathogenic, germline_truncating, germline_vus)]

    # keeping only some samples
    gcols = c('gene', 'variant.c', 'variant.p', 'vartype', 'CLNDN', 'CLNSIG', 'CLNVC', 'annotation', 'impact')
    fwrite(germline_dt_dedup[,..gcols], driver_germline_mutations)
    return(germline_dt_dedup)
}


#' @name annotate_germline_mutations
#' @title annotate_germline_mutations
#'
#' @description
#'
#' Read and process the germline mutations in coding genes
#'
#' @param germline_coding the "annotated.germline.coding.rds" output from the Alterations module
#' @param driver_germline_mutations path to output file in which to save germline mutations of driver genes
#' @param germline_genes list of germline genes that should be included in the output (if not provided then any gene annotated as "Germline" in CGC in addition to all DDR genes will be included
#' @param return.type (character) one of "data.table" or "GRanges"
#' 
#' @return data.table with node.fp, amp.fp, and ev.fp for node, event, and amplicon footprints
annotate_germline_mutations = function(germline_snpeff_snv_bcf, pathogenic = NULL, comvar = NULL,
                                       mask = NULL, name ='',
                                       germline_coding = 'annotated.germline.coding.rds'){
    if (file.exists(germline_snpeff_snv_bcf) & (file.info(germline_snpeff_snv_bcf)$size > 0)) {
        if (is.null(comvar)){
            # FIXME: hardcoded local path.
            comvar = '/gpfs/commons/groups/imielinski_lab/DB/modules/Alterations/db/comvar2.rds'
        }
        if (is.null(pathogenic)){
            # FIXME: hardcoded local path.
            pathogenic = '/gpfs/commons/groups/imielinski_lab/DB/modules/Alterations/db/pathogenic.rds'
        }
        if (is.null(mask)){
            # FIXME: hardcoded local path.
            mask = '/gpfs/commons/groups/imielinski_lab/DB/Broad/um75-hs37d5.bed.gz'
        }
        message('Reading common variants file: ', comvar)
        comvar.dt = readRDS(comvar)
        message('Reading pathogenic variants file: ', pathogenic)
        dt.patho = gr2dt(readRDS(pathogenic))
        # loading this list of genes to use to set the factor value for gene names. Notice that this is legacy from the original Alterations module. This is probably not necesary but I included it here for now to preseve compatibility with what was done in the Alterations module.
        ug = readLines(system.file('extdata', 'ug.txt', package = 'casereport'))

        message("Loading and Filtering germline call set")
        message("Filtering on heng-li mask, pathogenic variants (ClinVar), common variants")
        ## lofre = '(stop)|(frameshift)'
        lofre = paste0(paste0("(",
                              c("stop", "splice_acceptor", "splice_donor", "rare_amino_acid", "frameshift", "start_lost", "exon_loss", "transcript_ablation"),
                              ")"),
                       collapse = "|")

        temp_vcf = NA
        if (grepl('.bcf$', germline_snpeff_snv_bcf)){
            # to avoid collisions we include a random string in the temporary file name
            temp_vcf = paste0(tempdir(), '/', paste0(rand.string(), collapse = ''), '-germline_snpeff_snv.vcf')
            message('Making temporary VCF copy of the input BCF: ', temp_vcf)
            system(paste('bcftools index', germline_snpeff_snv_bcf))
            system(paste('bcftools view', germline_snpeff_snv_bcf, '-o', temp_vcf))
            germline_snpeff_snv_bcf = temp_vcf
        }
        som = gr.nochr(parsesnpeff(germline_snpeff_snv_bcf, geno = "GT"))
        if (NROW(som) == 0) {
            germ.mut = data.table(pair = character(0),
                                  gene = character(0),
                                  fpair = factor(levels = name),
                                  fgene = factor(levels = ug))
        } else {
            rm_col = sapply(mcols(som), function(x) inherits(x, c("list", "List")))
            colnames(mcols(som))[grepl("^GT", colnames(mcols(som)))] = "GT"
            mcols(som) = mcols(som)[, colnames(mcols(som))[!rm_col]]
            if (!is.null(mask)){
                message('Reading and marking SNVs according to the provided mask')
                mask.gr = import(mask)
                mcols(som)$masked.region = som %^% mask.gr
            }
            mcols(som)$file = NULL
            dt = as.data.table(som)[, pair := name]
            message('Marking pathogenic variants')
            dt = merge.data.table(dt,
                                  select(dt.patho, -strand, -width),
                                  by = c("seqnames", "start", "end", "REF", "ALT", "gene"),
                                  all.x = TRUE,
                                  allow.cartesian = TRUE)
            message('Marking common variants')
            dt = merge.data.table(dt,
                                  comvar.dt,
                                  by = c("seqnames", "start", "end", "REF", "ALT"),
                                  all.x = TRUE,
                                  allow.cartesian = TRUE)[, common_variant := replace_na(common_variant, FALSE)]
            if (!is.null(dt$QUAL) && !all(is.na(dt$QUAL))) {
                dt = dt[QUAL >= 5]
            }
            dt = dt[common_variant == FALSE][FILTER == "PASS"]
            if (NROW(dt) == 0) {
                ## germ.mut = NULL
                germ.mut = data.table(pair = character(0),
                                      gene = character(0),
                                      fpair = factor(levels = name),
                                      fgene = factor(levels = ug))
            } else {
                dt.patho = dt[grepl("[Pp]athogenic$", CLNSIG)]
                dt.patho = distinct(dt.patho, pair, seqnames, start, end, REF, ALT, gene, .keep_all = TRUE)
                dt.gtrunc = dt[!grepl("[Pp]athogenic$", CLNSIG) & (grepl(lofre, annotation))]
                if (!is.null(mask)){
                    dt.gtrunc = dt.gtrunc[masked.region == FALSE]
                }
                dt.gtrunc = distinct(dt.gtrunc, pair, seqnames, start, end, REF, ALT, gene, .keep_all = TRUE)
                dt.vus = dt[!grepl("[Pp]athogenic$", CLNSIG) | is.na(CLNSIG)] ## non-Pathogenic
                dt = rbind(dt.patho[, germline_pathogenic := TRUE], dt.gtrunc[, germline_truncating := TRUE], fill = T,
                           dt.vus[, germline_vus := TRUE])
                message(paste0(nrow(dt)," pathogenic germline variants found"))
                germ.mut = dt
                germ.mut[, fpair := factor(name)]
            }
        }
    } else {
        germ.mut = data.table(pair = character(0),
                              gene = character(0),
                              fpair = factor(levels = name),
                              fgene = factor(levels = ug))
    }
    saveRDS(germ.mut, germline_coding)
}

#' @name parsesnpeff
#' @title parse snpeff output into granges
#'
#'
#' @param vcf path to snpeff vcf
#' @param pad Exposed argument to skitools::ra.overlaps()
#' @return GRangesList of breakpoint pairs with junctions that overlap removed
#' @export
parsesnpeff = function (vcf, id = NULL, filterpass = TRUE, coding_alt_only = TRUE, 
    geno = NULL, gr = NULL, keepfile = FALSE, altpipe = FALSE, 
    debug = FALSE) 
{
    if (debug) 
        browser()
    out.name = paste0("tmp_", rand.string(), ".vcf.gz")
    tmp.path = paste0(tempdir(), "/", out.name)
    if (!keepfile) 
        on.exit(unlink(tmp.path))
    try2({
        catcmd = if (grepl("(.gz)$", vcf)) "zcat" else "cat"
        onepline = "/gpfs/commons/groups/imielinski_lab/git/mskilab/flows/modules/SnpEff/source/snpEff/scripts/vcfEffOnePerLine.pl"
        if (coding_alt_only) {
            filt = "java -Xmx20m -Xms20m -XX:ParallelGCThreads=1 -jar /gpfs/commons/groups/imielinski_lab/git/mskilab/flows/modules/SnpEff/source/snpEff/SnpSift.jar filter \"( ANN =~ 'chromosome_number_variation|exon_loss_variant|rare_amino_acid|stop_lost|transcript_ablation|coding_sequence|regulatory_region_ablation|TFBS|exon_loss|truncation|start_lost|missense|splice|stop_gained|frame' )\""
            if (filterpass)
                cmd = sprintf(paste(catcmd, "%s | %s | %s | bcftools view -i 'FILTER==\"PASS\"' | bgzip -c > %s"), 
                  vcf, onepline, filt, tmp.path)
            else cmd = sprintf("cat %s | %s | %s | bcftools norm -Ov -m-any | bgzip -c > %s", 
                vcf, onepline, filt, tmp.path)
        }
        else {
            filt = ""
            if (filterpass) 
                cmd = sprintf(paste(catcmd, "%s | %s | bcftools view -i 'FILTER==\"PASS\"' | bgzip -c > %s"), 
                  vcf, onepline, tmp.path)
            else cmd = sprintf(paste(catcmd, "%s | %s | bcftools norm -Ov -m-any | bgzip -c > %s"), 
                vcf, onepline, tmp.path)
        }
        system(cmd)
    })
    if (!altpipe) 
        out = grok_vcf(tmp.path, long = TRUE, geno = geno, gr = gr)
    else {
        vcf = readVcf(tmp.path)
        vcf = S4Vectors::expand(vcf)
        rr = within(rowRanges(vcf), {
            REF = as.character(REF)
            ALT = as.character(ALT)
        })
        ann = as.data.table(tstrsplit(unlist(info(vcf)$ANN), 
            "\\|"))[, 1:15, with = FALSE, drop = FALSE]
        fn = c("allele", "annotation", "impact", "gene", "gene_id", 
            "feature_type", "feature_id", "transcript_type", 
            "rank", "variant.c", "variant.p", "cdna_pos", "cds_pos", 
            "protein_pos", "distance")
        data.table::setnames(ann, fn)
        if ("AD" %in% names(geno(vcf))) {
            adep = setnames(as.data.table(geno(vcf)$AD[, , 1:2]), 
                c("ref", "alt"))
            gt = geno(vcf)$GT
        }
        else if (all(c("AU", "GU", "CU", "TU", "TAR", "TIR") %in% 
            c(names(geno(vcf))))) {
            this.col = dim(geno(vcf)[["AU"]])[2]
            d.a = geno(vcf)[["AU"]][, , 1, drop = F][, this.col, 
                1]
            d.g = geno(vcf)[["GU"]][, , 1, drop = F][, this.col, 
                1]
            d.t = geno(vcf)[["TU"]][, , 1, drop = F][, this.col, 
                1]
            d.c = geno(vcf)[["CU"]][, , 1, drop = F][, this.col, 
                1]
            mat = cbind(A = d.a, G = d.g, T = d.t, C = d.c)
            rm("d.a", "d.g", "d.t", "d.c")
            refid = match(as.character(VariantAnnotation::fixed(vcf)$REF), colnames(mat))
            refid = ifelse(!isSNV(vcf), NA_integer_, refid)
            altid = match(as.character(VariantAnnotation::fixed(vcf)$ALT), colnames(mat))
            altid = ifelse(!isSNV(vcf), NA_integer_, altid)
            refsnv = mat[cbind(seq_len(nrow(mat)), refid)]
            altsnv = mat[cbind(seq_len(nrow(mat)), altid)]
            this.icol = dim(geno(vcf)[["TAR"]])[2]
            refindel = d.tar = geno(vcf)[["TAR"]][, , 1, drop = F][, 
                this.icol, 1]
            altindel = d.tir = geno(vcf)[["TIR"]][, , 1, drop = F][, 
                this.icol, 1]
            adep = data.table(ref = coalesce(refsnv, refindel), 
                alt = coalesce(altsnv, altindel))
            gt = NULL
        }
        else {
            message("ref and alt count columns not recognized")
            adep = NULL
            gt = NULL
        }
        mcols(rr) = BiocGenerics::cbind(mcols(rr), ann, adep, 
            gt = gt[, 1])
        out = rr
    }
    this.env = environment()
    return(this.env$out)
}
