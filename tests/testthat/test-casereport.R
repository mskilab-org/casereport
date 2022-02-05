test_that("Run case report for G32831", {

    opt = list(
             pair = 'G32831',
             outdir = tempdir(),
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

    wgs.report(opt)



#    library(Flow)
#    task.file = system.file('extdata', 'case_report_module/wgs.casereport.task', package = 'casereport')
#    task.content = readLines(task.file)
#    module_path = system.file('extdata', 'case_report_module', package = 'casereport')
#    task.content = gsub('MODULE_PATH', module_path, task.content)
#    new.task.file = paste0(tmpdir, '/wgs.casereport.task')
#    writeLines(task.content, new.task.file)
#
#    # FIXME: we need to find 
#    jb = Flow::Job(new.task.file,
#                   gdt,
#                   rootdir = tmpdir)
#
#    srun(jb)

})
