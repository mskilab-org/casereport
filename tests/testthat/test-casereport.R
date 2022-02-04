test_that("Run case report for G32831", {
    gdt = data.table(
             pair = 'G32831',
             jabba_rds = system.file('extdata', 'G32831/jabba.simple.rds', package = "casereport"),
             cbs_cov_rds = system.file('extdata', 'G32831/G32831.Detroit_562.chr11.window.rds', package = "casereport"),
             fusions = system.file('extdata', 'G32831/fusions.rds', package = "casereport"),
             complex = system.file('extdata', 'G32831/complex.rds', package = "casereport"),
             het_pileups_wgs = system.file('extdata', 'G32831/G32831.Detroit_562.chr11.window.hets.txt', package = "casereport"),
             deconstruct_sigs = system.file('extdata', 'G32831/deconstructSigs.out.rds', package = "casereport"),
             deconstruct_variants =  system.file('extdata', 'G32831/variants.out.txt.gz', package = "casereport"),
             snv_vcf = system.file('extdata', 'G32831/G32831.Detroit_562.vcf', package = "casereport"),
             snpeff_snv_bcf = system.file('extdata', 'G32831/annotated.bcf', package = "casereport"))
})
