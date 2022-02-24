##### Utility functions to process germline mutations #####
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
#' @param return.type (character) one of "data.table" or "GRanges"
#' 
#' @return data.table with node.fp, amp.fp, and ev.fp for node, event, and amplicon footprints
process_germline_muts = function(germline_coding, driver_germline_mutations, germline_genes = NULL){
    if (!file.good(germline_coding)){
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
