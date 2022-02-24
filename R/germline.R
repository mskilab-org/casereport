##### Utility functions to process germline mutations #####

read_germline_muts = function(report.config){
    if (!file.good(report.config$germline.coding.rds)){
        empty_germline_dt = data.table(gene = as.character(),
                                 variant.c = as.character(),
                                 variant.p = as.numeric(),
                                 vartype = as.numeric(),
                                 CLNDN = as.character(),
                                 CLNSIG = as.character(),
                                 CLNVC = as.character(),
                                 annotation = as.character(),
                                 impact = as.character())
        return(empty_empty_germline_dt)
    }

    germline_dt = readRDS(report.config$germline.coding.rds)

    germline_dt_dedup[, uid := paste(gene, variant.c)]

    # in case the AA variation is different in different transcripts we will take all versions
    germline_dt_dedup[, variant.p := paste(unique(variant.p), collapse = ','), by = uid]

    # deduping
    germline_dt_dedup = germline_dt_dedup[!duplicated(uid)]
    # keeping only genes annotated 
    cgc = fread(system.file('extdata', 'cgc.tsv', package = 'casereport'))
    germline_genes = cgc[Germline == 'yes', get("Gene Symbol")]
    germline_dt_dedup = germline_dt_dedup[gene %in% germline_genes]
    # order so that pathogenic mutations appear first, then germline, and lastly VUS
    germline_dt_dedup = germline_dt_dedup[order(germline_pathogenic, germline_truncating, germline_vus)]

    # keeping only some samples
    gcols = c('gene', 'variant.c', 'variant.p', 'vartype', 'CLNDN', 'CLNSIG', 'CLNVC', 'annotation', 'impact')
    fwrite(germline_dt_dedup[,..gcols], report.config$germline.coding.summary)
    return(germline_dt_dedup)
}
