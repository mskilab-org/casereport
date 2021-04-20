casereport = function(thisp){

    message(paste0("Beginning Case Analysis for ", thisp$pair))

    onepline = "/gpfs/commons/groups/imielinski_lab/git/mskilab/flows/modules/SnpEff/source/snpEff/scripts/vcfEffOnePerLine.pl"

    library(gGnome)
    dt.patho = gr2dt(readRDS("/gpfs/commons/groups/imielinski_lab/DB/ClinVar/pathogenic.rds"))

    filt = "java -Xmx1g -Xms1g -jar /gpfs/commons/groups/imielinski_lab/git/mskilab/flows/modules/SnpEff/source/snpEff/SnpSift.jar filter \"( ANN =~ 'missense|splice|stop_gained|frame' )\""

    gff = readRDS(gzcon(url('http://mskilab.com/gGnome/hg19/gencode.v19.annotation.gtf.gr.rds')))
    ge = read_gencode('gene'); message('loaded GENCODE')
    pge = ge %Q% (gene_type == 'protein_coding' & transcript_status == "KNOWN");
    cds = read_gencode('cds'); message('loaded CDS')
    path = paste0(getwd(),"/casereports/",thisp$pair,"/")
    pge = genes %Q% (gene_type == 'protein_coding') %Q% (level <3 & transcript_status == 'KNOWN')
    cgc = fready('/gpfs/commons/groups/imielinski_lab/DB/COSMIC/CGC.txt')
    cgc$onc = cgc$Significant_across_any_tumor_type_lineage_if_gain_ == 'Yes'
    cgc$tsg = cgc$Significant_across_any_tumor_type_lineage_if_lost_ == 'Yes'
    pge.cgc = pge %Q% (gene_name %in% cgc[onc == TRUE | tsg == TRUE, Gene_Symbol])
    cgc_tsg = pge %Q% (gene_name %in% cgc[tsg == TRUE, Gene_Symbol])
    cgc_onc = pge %Q% (gene_name %in% cgc[onc == TRUE, Gene_Symbol])
    cgc_fus = pge %Q% (gene_name %in% cgc[nchar(Translocation_Partner) > 0, Gene_Symbol])
    Sys.setenv(DEFAULT_BSGENOME = "/gpfs/commons/groups/imielinski_lab/DB/UCSC/hg19.broad.chrom.sizes")

    message(paste0("Returning Outputs to Directory: ", path))
    system(paste0("mkdir -p ",path))

    ##Load and Print COSMIC Signatures
    if(!is.na(thisp$signature_counts)){
        message("Returning top 10 COSMIC Signature Counts")
        sig = try(fread(thisp$signature_counts))
        cols = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a") ##Make a color map real quick to match to the names that I will want to color by

        names(cols) = sig$Signature[1:10]
        barplot = ggplot(sig[1:10],
                         aes(x = Signature, y = num_events)) +
            geom_col(aes(fill = Signature)) +
            ylim(0,max(sig$Signature)) +
            scale_fill_manual(values = cols,aesthetics = "fill") +
            labs(fill = " ",
                 x = "Cosmic SNV Signature",
                 y = "SNV Contributing Counts",
                 title = paste0("COSMIC SNV Signature Contributions to ", thisp$pair[1]),
                 size = 20) +
            theme(legend.position = "none",
                  axis.text.x = element_text(size = 10,
                                             angle = 45,
                                             vjust = 0.5),
                  axis.text.y = element_text(size = 15)
                  )
        ppng(print(barplot), filename = paste0(path,"COSMIC_signatures.png"))
    } else {
        message(paste0("No Signatures generated for ",thisp$pair))
    }

    ##Load and print CBS Ratio Coverage
    kar = karyogram()
    message("Defining covcbs and covgt functions")
    covcbs = function(x){gTrack(readRDS(x),
                                y.field = "ratio",
                                name = " ",
                                max.ranges = 1e4,
                                circles = TRUE,
                                lwd.border = 0.2)}
    covgts = function(x){gTrack(readRDS(x),
                                y.field = "seg.mean",
                                name = " ",
                                max.ranges = 1e4,
                                lwd.border = 0.2)}
    covgtn = function(x){gTrack(readRDS(x),
                                y.field = "mean",
                                name = " ",
                                max.ranges = 1e4,
                                lwd.border = 0.2)}

    if(!is.na(thisp$cbs_cov_rds)){
        message("Loading gTracks")
        cbs.gt = covcbs(thisp$cbs_cov_rds)
        cbs.gt$name = "CBS Ratio"
        cov.t.gt = covgts(thisp$cbs_seg_rds)
        cov.t.gt$name = "Segmentation "
        cov.n.gt = covgtn(thisp$cbs_nseg_rds)
        cov.n.gt$name = "Normal Segment Mean"
        win = si2gr(hg_seqlengths())[c(-25,-26)]
        message("Printing Coverage Plots")
        ppdf(
            lapply(1:length(win),function(x){
                plot.new()
                plot(c(kar,cov.n.gt,cov.t.gt,cbs.gt),
                     win[x])}),
            width = 20,
            title = paste0("Genome Scale Coverage for ", thisp$pair),
            filename = paste0(path,"Coverage_profile.pdf"))
        message("Done Printing Coverage")
    } else {
        message(paste0("No CBS data for ", thisp$pair))
    }

    ##Load and print JaBbA findings
    if(!is.na(thisp$jabba_rds)){
        message("Printing JaBbA Genome Scale View")
        gg = gG(jabba = thisp$jabba_rds)
        message("Loaded JaBbA Object")
        gt = gg$gt
        cbs.gt = covcbs(thisp$cbs_cov_rds)
        cbs.gt$name = "CBS Ratio"
        ppdf(lapply(1:length(win),function(x){
            plot.new()
            plot(c(kar,cbs.gt,gt),
                 win[x])}),
            width = 20,
            title = paste0("JaBbA Genome Scale Model: ", thisp$pair),
            filename = paste0(path,"JaBbA.global.pdf"))
        message("Done printing JaBbA genome scale")
        message("Now printing complex events")
        message("Identifying and printing chromothripsis events")
        gg = chromothripsis(gg)
        ct = gg[chromothripsis > 0]$copy$subgraph(d = 1e4)
        if(length(ct) != 0){
            ppdf(plot(c(kar,cbs.gt,gt),ct$footprint+1e4),
                 width = 15,
                 title = paste0("Chromothripsis Events in: ", thisp$pair),
                 filename = paste0(path,"JaBbA.chromothripsis.pdf"))
        } else {message("No Chromothripsis found")}
        message("Identifying and printing bfb events")
        gg = bfb(gg)
        bfb = gg[bfb > 0]$copy$subgraph(d = 1e4)
        if(length(bfb) != 0){
            ppdf(plot(c(kar,gt.ge,cbs.gt,gt),gg$nodes[bfb>0]$edges$footprint+1e4),
                 width = 15,
                 title = paste0("Breakage Fusion Bridges in: ", thisp$pair),
                 filename = paste0(path,"JaBbA.bfb.pdf"))
        } else {message("No BFBs found")}
        message("Identifying and printing Templated insertion Chains")
        gg = tic(gg)
        tic = gg[tic > 0]$copy$subgraph(d = 1e4)
        if(length(tic) != 0){
            ppdf(plot(c(kar,cbs.gt,gt),
                      gg$nodes[tic > 0]$edges$footprint +1e4),
                 width = 15,
                 title = paste0("Templated Insertion Chains in: ", thisp$pair),
                 filename = paste0(path,"JaBbA.tic.pdf"))
        } else {message("No TICs found")}
        message("Identifying and printing Pyrgos")
        gg = pyrgos(gg)
        pyrgo = gg[pyrgo > 0]$copy$subgraph(d = 1e4)
        if(length(pyrgo) != 0){
            ppdf(plot(c(kar,cbs.gt,gt),pyrgo$footprint+1e4),
                 width = 15,
                 title = paste0("Pyrgos in: ", thisp$pair),
                 filename = paste0(path,"JaBbA.pyrgos.pdf"))
        } else {message("No Pyrgos found")}
        message("Identifying and printing Rigmas")
        gg = rigma(gg)
        rigma = gg[rigma > 0]$copy$subgraph(d = 1e4)
        if(length(rigma) != 0){
            ppdf(plot(c(kar,gt.ge,cbs.gt,gt),gg$nodes[rigma > 0]$edges$footprint+1e4),
                 width = 15,
                 title = paste0("Rigmas in: ", thisp$pair),
                 filename = paste0(path,"JaBbA.rigma.pdf"))
        } else {message("No Rigmas found")}
        message("Calling Fusions")
        fus = gGnome::fusions(gg,gff)
        if(length(fus) != 0){
            dt = fus$dt[in.frame == TRUE][silent == FALSE]
            dt = dt[ ,c(1:19,22:23)]
            write.csv(dt, paste0(path, "fusions.summary.csv"))
            fus = gGnome::fus[!duplicated(fus$dt$name)]
            message("Printing Fusions")
            ppdf(lapply(1:length(fus),function(x){
                plot.new()
                plot(c(cbs.gt,gt,fus$gt),fus[x]$footprint+1e5)}),
                width = 15,
                title = paste0("Fusions in: ", thisp$pair),
                filename = paste0(path,"Fusion.plots.pdf")
                )
        } else {message("No fusions found")}
        message("Calling CNAs from JaBbA Output")
        seg = gg$nodes$dt[, pair := thisp$pair]
        seg[, ploidy := sum(width*cn, na.rm = TRUE)/sum(width[!is.na(cn)], na.rm = TRUE), by = pair]
        seg[, cn.rel := cn/ploidy]
        seg[, amp := cn.rel>2]
        seg[, homdel := cn==0]
        seg[, hetdel := cn==1]
        seg.gr = dt2gr(seg[!is.na(cn), ])
        seg.gr.onc = cgc_onc[, "gene_name"] %*% (seg.gr %Q% (amp == TRUE)) %Q% (!duplicated(cbind(gene_name, pair)))
        seg.gr.onc$type = 'amp'
        seg.gr.tsg = cgc_tsg[, "gene_name"] %*% (seg.gr %Q% (hetdel | homdel)) %Q% (!duplicated(cbind(gene_name, pair)))
        seg.gr.tsg$type = ifelse(seg.gr.tsg$homdel, 'homdel', ifelse(seg.gr.tsg$hetdel, 'hetdel', NA))
        cna = seg2gr(seg[, alt := ifelse(amp, 'amp', ifelse(hetdel | homdel, 'del', 'wt'))][alt != 'wt', ])[, c('cn', 'cn.rel', 'alt', 'pair')]
        cnag = cna %*% (pge.cgc[, 'gene_name'])
        dt.cna = gr2dt(cnag)
        write.csv(dt.cna, paste0(path, "cna.genes.csv"))
    } else {
        message("No JaBbA complete for this sample")
    }

                                        #Parse and return SNPEFF calls
    if(!is.na(thisp$annotated_vcf_somatic)){
        message("Parsing SNPEFF Somatic Calls")
        vcf_path = thisp$annotated_vcf_somatic  ##Fix your input
        out.name = paste0("tmp_", rand.string(), ".vcf")
        tmp.path = paste0(tempdir(), "/", out.name)
        cmd = sprintf("cat %s | %s | %s > %s", vcf_path, onepline, filt, tmp.path)
        system(cmd)
        fn = gsub('\\.vcf$', 'coding.somatic.rds', vcf_path)  ##Choose where to write it
        vcf = grok_vcf(tmp.path, long = TRUE)
        vcf = vcf %Q% (!impact %in% c("LOW", "MODIFIER"))
        unlink(tmp.path)
        dt = gr2dt(vcf)
        dt = dt[,-c("RU","LOF","NMD")]
        write.csv(dt, paste0(path,"Somatic.mutations.csv"))
        message("Parsing SNPEFF Germline Calls")
        vcf_path = thisp$annotated_vcf_germline  ##Fix your input
        out.name = paste0("tmp_", rand.string(), ".vcf")
        tmp.path = paste0(tempdir(), "/", out.name)
        cmd = sprintf("cat %s | %s | %s > %s", vcf_path, onepline, filt, tmp.path)
        system(cmd)
        fn = gsub('\\.vcf$', 'coding.germline.rds', vcf_path)  ##Choose where to write it
        vcf = grok_vcf(tmp.path, long = TRUE)
        vcf = vcf %Q% (!impact %in% c("LOW", "MODIFIER"))
        unlink(tmp.path)
        dt.g = gr2dt(vcf)
        dt.g = merge(dt.g, copy(dt.patho)[, pathogenic := TRUE],
                     by = c("seqnames","start","end","REF"),
                     all.x = TRUE,
                     allow.cartesian = TRUE)[, pathogenic := replace_na(pathogenic, FALSE)]
        dt.g = dt.g[pathogenic == TRUE]
        dt.g = dt.g[,-c("RU","LOF","NMD")]
        write.csv(dt.g, paste0(path,"Germline.mutations.csv"))
        dt.union = rbind(dt, dt.g, fill = TRUE)
        dt.union = dt.union[,-c("RU","LOF","NMD")]
        write.csv(dt.union, paste0(path, "All.mutations.csv"))
    } else{
        message("No SNPEFF Annotations Available")}
}



