

#' @name get_pmkb_dt
#' @title get_pmkb_dt
#'
#' @description
#'
#' Internal function to read the pmkb TSV (this is the simplified version that we produced from the PMKB interpretations CSV)
#'
#' @param pmkb_tsv (character) path to PMKB TSV (if not provided then the default table supplied with the package will be used
get_pmkb_dt = function(pmkb_tsv = NA){
    if (is.na(pmkb_tsv)){
        pmkb_tsv = system.file("extdata", "pmkb-tier.csv", package = "casereport")
        message('Using default PMKB annotations: ', pmkb_tsv)
    }
        pmkb.dt = fread(pmkb_tsv)
        return(pmkb.dt)
}

#' @name get_pmkb_tier_table
#' @title get_pmkb_tier_table
#'
#' @description
#'
#' Internal function to generate a table with a single tier value per gene
#' 
#' The output table includes columns: gene, type, and tier
#' the type takes one of the following values: 'ONC', 'TSG', or 'ONC|TSG'. As of writing this the only gene with 'ONC|TSG' annotation was CDH1
#'
#' @param pmkb_tsv (character) path to PMKB TSV (if not provided then the default table supplied with the package will be used
get_pmkb_tier_table = function(pmkb_tsv = NA){
    pmkb.dt = get_pmkb_dt(pmkb_tsv)
    return(pmkb.dt[, .(type = paste(sort(unique(gene.type)), collapse = '|'), tier = min(Tier)), by = gene])
}



#' @name annotate_with_pmkb
#' @title annotate_with_pmkb
#'
#' @description
#'
#' annotate driver mutations with tier and gene.type information from the PMKB table
#'
#' @param pmkb_tsv (character) path to PMKB TSV (if not provided then the default table supplied with the package will be used
#' @param driver.mutations.dt (data.table) must contain column "gene"
annotate_with_pmkb = function(driver.mutations.dt, pmkb_tsv = NA){
    pmkb.dt = get_pmkb_dt(pmkb_tsv)
    driver.mutations.dt.annotated.with.tier = merge.data.table(pmkb.dt[, .(gene.type = unique(gene.type), Tier = min(Tier)), by = gene], driver.mutations.dt, by = 'gene', all.y = TRUE)
    return(driver.mutations.dt.annotated.with.tier[order(Tier)])
}

#' @name extract_variant_types_from_pmkb
#' @title extract_variant_types_from_pmkb
#'
#' @description
#'
#' Internal function to parse the PMKB interpretations CSV and extract unique variant types
#'
#' We used this function to generate pmkb_tsv and it is here for the sake of documentation.
#' pmkb_tsv is a data.table with the following columns: gene, variant, type, position.type, any.variant, gene.type.
#' gene - gene name (this column could contain duplicates since there will be a row for each variant value per gene
#' variant - a single variant description (the original PMKB interpretations CSV could have multiple variant description per line, separated by "|", we split these so that each variant appears in a separate line).
#' type - one of the following: missense, nonsense, frameshift, any mutation, insertion, deletion, copy number loss, stop gain, indel, copy number gain, rearrangement
#' position.type - either NA, codon, or exon. This is column is intended to identify these types of variants: "T618I" (type = codon), "codon(s) 515 missense" (type = codon), "exon(s) 34 frameshift" (type = exon)
#' any.variant - TRUE/FALSE. intended to identify variants such as: "any nonsense", "codon(s) 515 missense", "exon(s) 26-27 deletion", "any mutation"
#' gene.type - "TSG"/"ONC". We used a simple heuristic to annotate each gene as either ONC or TSG: see annotate_pmkb_entries_as_tsg_or_onc below for details. 
#'
#' @param pmkb.csv (character) path to PMKB interpretations CSV (if not provided then the default table supplied with the package will be used). The table was downloaded from https://pmkb.weill.cornell.edu/about on 6/11/2021.
#' @param output_file (optional character) if supplied then the PMKB TSV will be saved to this path
#' @return data.table 
extract_variant_types_from_pmkb = function(pmkb.csv = NA, output_file = NA){
    if (is.na(pmkb.csv)){
        message('No PMKB interpretations CSV input specified.')
        pmkb.csv = system.file("extdata", "pmkb-interpretations-06-11-2021.csv", package = "casereport")
        message('Using default PMKB interpretations CSV: ', pmkb.csv)
    }
    pmkb.interpretations = fread(pmkb.csv)

    bla = unique(unlist(strsplit(pmkb.interpretations[, get('Variant(s)')], '\\|')))
    # build a data.table with gene name and mutation description
    bli = strsplit(bla, ' ')
    genes = sapply(bli, '[', 1)

    # remove the gene name from the variant description
    vars = sapply(lapply(bli, '[', -1), function(x) paste(x, collapse = ' '))

    pmkb.dt = data.table(gene = genes, variant = vars)

    # Add variant type for each variant
    # we do so by a series of conditions to deal with all the PMKB entry types
    pmkb.dt[grepl('rearrangement', variant), type := variant]

    pmkb.dt[, any.variant := F] # column to deal with "any" anntations
    pmkb.dt[grepl('any', variant), any.variant := T]
    # extract the mutation type for "any" annotations (e.g. "any missense" is of type = missense)
    pmkb.dt[grepl('^any', variant), type := sapply(strsplit(variant, ' '), function(x){x[2]})]
    pmkb.dt[grepl('any$', variant), type := sapply(strsplit(variant, ' '), function(x){x[2]})]

    # some variants are limited to specific codons or exons
    # we add the position.type column to identify these
    pmkb.dt[grepl('codon', variant), position.type := 'codon']
    # these are never specific (so these are all any.variant == TRUE
    pmkb.dt[grepl('codon', variant), any.variant := TRUE]
    pmkb.dt[grepl('exon', variant), unique(variant)]
    pmkb.dt[grepl('exon', variant), position.type := 'exon']
    pmkb.dt[grepl('exon', variant), any.variant := TRUE]
    # this is a bit tricky: we remove all numeric and commas and the position name (i.e. codon(s) or exons(s)) from the variant and then we get the variant type
    pmkb.dt[position.type %in% c('codon', 'exon'), type := gsub('codon\\(s\\)|exon\\(s\\)|\\,|[0-9]|\\ |\\-', '', variant)]

    # deal with single amino acid variants (e.g. "G67C")
    pmkb.dt[grepl('[A-Z][0-9]+[A-Z]', variant), type := 'missense']
    pmkb.dt[grepl('[A-Z][0-9]+[A-Z]', variant), position.type := 'codon']

    pmkb.dt[grepl('\\*', variant), type := 'stop gain']

    # deal with CNVs
    pmkb.dt[grepl('copy', variant), type := variant]

    pmkb.dt[grepl('del', variant), type := 'deletion']
    pmkb.dt[grepl('ins', variant), type := 'insertion']
    pmkb.dt[grepl('delins|indel', variant), type := 'indel']
    pmkb.dt[grepl('fs', variant), type := 'frameshift']

    # "any mutation" is in some cases listed as "any" and somtimes listed as "any mutation"
    # here we catch both options an annotate as type "any mutation"
    pmkb.dt[grepl('any|mutation', type), type := 'any mutation']

    # add Tier data to the pmkb.dt
    pmkb.dt = merge.data.table(pmkb.dt, pmkb.interpretations[, .(Tier = min(Tier)), by = Gene], by.x = 'gene', by.y = 'Gene', all.x = T)

    # lastely, we add the gene.type column ("TSG" or "ONC")
    pmkb.dt = annotate_pmkb_entries_as_tsg_or_onc(pmkb.dt)

    if (!is.na(output_file)){
        fwrite(pmkb.dt, output_file, sep = '\t')
    }
    return(pmkb.dt)
}

#' @name annotate_pmkb_entries_as_tsg_or_onc
#' @title annotate_pmkb_entries_as_tsg_or_onc
#'
#' @description
#'
#' Internal function to annotate genes in the PMKB table as either TSG or ONC
#'
#' We used this function to generate pmkb_tsv and it is here for the sake of documentation.
#'
#' @param pmkb (data.table) the pmkb data.table
#' @return data.table with added TSG/ONC annotation for each gene
annotate_pmkb_entries_as_tsg_or_onc = function(pmkb){
    # Important notice: the heuristic implemented here is very naive
    # and hence some ONCs might get a "TSG" label.
    # for example I noticed that KIT gets labeled as a TSG because in at least one cancer type "any mutation" qualifies in the PMKB table
    # this is not so bad, since any way if a gene is already listed in our onc and tsg lists then those label prevail (see the code for "ONC and TSG annotations" in wgs.run.R)
    if (!inherits(pmkb, 'data.table')){
        tryCatch({
                 pmkb = fread(pmkb)
        },
        error = function(e){
            message('The PMKB interpretations must be provided as a data.table or a path to a tabular file. While trying to read the PMKB interpretations we got this error:')
            stop(e)
        }) 
    }

    tsg_variant_types = c('any mutation', 'copy number loss', 'frameshift', 'nonsense', 'stop gain')
    tsgs = unique(pmkb[is.na(position.type)][type %in% tsg_variant_types, gene])

    pmkb[, gene.type := 'ONC'] # anything that is not a TSG is an ONC
    pmkb[gene %in% tsgs, gene.type := 'TSG']
    # in addition, if a gene accepts any missense or indel then it is a TSG
    pmkb[gene.type != 'TSG' & is.na(position.type) & any.variant == TRUE, gene.type := 'TSG']
}

# 
# in the future we can expend this with more parsers
mutation.tier.annotators = list('PMKB' = annotate_with_pmkb)
