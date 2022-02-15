
nm = 'G32831'
case.dir = paste0(tempdir(), '/casereport-sv-tests')
system(paste('rm -rf', case.dir))
dir.create(case.dir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
opt = list(
         pair = nm,
         outdir = case.dir,
         jabba_rds = system.file('extdata', 'G32831/jabba.simple.rds', package = "casereport"),
         cbs_cov_rds = system.file('extdata', 'G32831/G32831.Detroit_562.chr11.window.rds', package = "casereport"),
         fusions = system.file('extdata', 'G32831/fusions.rds', package = "casereport"),
         complex = system.file('extdata', 'G32831/complex.rds', package = "casereport"),
         het_pileups_wgs = system.file('extdata', 'G32831/G32831.Detroit_562.chr11.window.hets.txt', package = "casereport"),
         deconstruct_sigs = system.file('extdata', 'G32831/deconstructSigs.out.rds', package = "casereport"),
         deconstruct_variants = system.file('extdata', 'G32831/variants.out.txt.gz', package = "casereport"),
         tpm = system.file('extdata', 'G32831/tpm.txt.gz', package = "casereport"),
         snv_vcf = system.file('extdata', 'G32831/G32831.Detroit_562.vcf', package = "casereport"),
         snpeff_snv_bcf = system.file('extdata', 'G32831/annotated.bcf', package = "casereport"))

gencode_fn = system.file('extdata', 'mock.gencode.gtf', package = 'casereport')
test_that("create_gene_gtrack", {
    gngt = create_genes_gtrack('TP53',
                        cached.dir = opt$outdir,
                        gencode.fname = gencode_fn)
    gn.gt.fname = paste0(opt$outdir, '/gencode.composite.collapsed.rds')
    expect_true(file.exists(gn.gt.fname))
    expect_true(inherits(gngt, 'gTrack'))
    expect_true('gene_label' %in% names(values(gngt@data[[1]])))
})
