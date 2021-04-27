#' zchoo Tuesday, Apr 27, 2021 10:49:13 AM
#' this is to generate data tables and plots for fusions

#' @name fusion.table
#' @title fusion.table
#'
#' @description
#'
#' Prepare table of fusion events for displaying
#'
#' @param fusions.fname (character) file name containing fusion gWalks
#' @param complex.fname (character) file name to events
#' @param ccg.fname (character) cancer gene census file name
#' @param outdir (character) output directory
#'
#' @return data.table with columns:
#' - walk.id
#' - name (e.g. involved genes)
#' - num.aa (total number of aminos)
#' - gene.pc (protein coordinates)
#' - driver (logical, is a gene a driver)
#' - ev.id (complex event IDs overlapping with walk)
#' - ev.type (complex event type overlapping with walk)
fusion.table = function(fusions.fname = NULL,
                        complex.fname = NULL,
                        ccg.fname = "/data/ccg.txt",
                        outdir = "./") {

    if (!file.exists(fusions.fname)) {
        stop("fusions.fname does not exist")
    }
    if (!file.exists(complex.fname)) {
        stop("complex.fname does not exist")
    }

    ## filter to include only in-frame non-silent
    this.fusions = readRDS(fusions.fname) ## gWalk object
    filtered.fusions = this.fusions[in.frame == TRUE & silent == FALSE]

    ## compute total number of amino acids and mark
    n.aa = sapply(filtered.fusions$dt$gene.pc, function(pc) { sum(width(parse.grl(pc))) })
    filtered.fusions$set(total.aa = n.aa)

    ## filter so that gene pairs are unique (if multiple choose the one with most AA's)
    filtered.fusions = filtered.fusions[order(n.aa, decreasing = TRUE)]
    filtered.fusions = filtered.fusions[which(!duplicated(filtered.fusions$dt$name))]

    ## then 
}

#' @name fusion.plot
#' @title fusion.plot
#'
#' @description
#'
#' Create .png files for each fusion
#'
#' @param fusions.fname (character) file name containing fusion gWalks
#' @param cvgt.fname (character) coverage gTrack file name
#' @param pad (numeric) gWalk pad for plotting default 1e5
#' @param height (numeric) plot height default 1e3
#' @param width (numeric) plot width default 1e3
#' @param outdir (character) output directory
#'
#' @return data.table with columns walk.id, plot.fname (input to slickR)
fusion.plot = function(fusions.fname = NULL,
                       cvgt.fname = NULL,
                       pad = 1e5,
                       height = 1e3,
                       width = 1e3,
                       outdir = "./") {
    return()
}
