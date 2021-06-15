

# take the 
get_pmkb_dt = function(pmkb_tsv){
    if (is.na(pmkb_tsv)){
        pmkb_tsv = file.path(opt$libdir, "data", "pmkb-tier.tsv")
        message('Using default PMKB annotations: ', pmkb_tsv)
    }
        pmkb.dt = fread(pmkb_tsv)
        return(pmkb.dt)
}


pmkb_parser = function(driver.mutations.dt, pmkb_tsv = NA){

    pmkb.dt = get_pmkb_dt(pmkb_tsv)

    driver.mutations.dt.annotated.with.tier = merge.data.table(driver.mutations.dt, pmkb.dt[, .(Tier = min(Tier)), by = gene], by = 'gene')

    return(driver.mutations.dt.annotated.with.tier[order(Tier)])
}

# in the future we can expend this with more parsers
mutation.tier.parsers = list('PMKB' = pmkb_parser)
