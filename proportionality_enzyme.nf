#!/usr/bin/env nextflow


log.info """\
         PROPORTIONALITY - TEST ON ENZYME GENES
         ======================================="
         Input dataset ID                               : ${params.id}
         Input dataset file (DATASET)                   : ${params.data}
         Input genes (FILE)                             : ${params.genes}
         Metric (str)                                   : ${params.metric}
         Transformation (str)                           : ${params.ivar}
         Normalization (str)                            : ${params.norm}
         If specific tissue                             : ${params.tis}
         If get n random samples (bool)                 : ${params.usensample}
         Number samples (int)                           : ${params.nsamples}
         Random sample seeds                            : ${params.seeds}
         If use all samples (bool)                      : ${params.allsample}
         Number permutation (int)                       : ${params.permutation}
         Output directory (DIRECTORY)                   : ${params.outdir}
         Trace report (DIRECTORY)                       : ${params.tracedir}
         If check file already exists (bool)            : ${params.check}
         """
         .stripIndent()

/*
 *  Manage channels
 */


// nsamp channel

// set nsample seeds channel
Channel
  .fromPath(params.seeds)
  .map { item -> [ item.baseName.tokenize('_')[1], item ] } 
  .combine( Channel.from(params.nsamples) )
  .into { ch_seeds1; ch_seeds2 }     // Eg. [seeds2, /path/to/nsamp_seeds2_singularity, 100]



/* 1. Proportionality-enzyme with all samples */
process propr_allsamp {

    memory '84.GB'

    tag "allsamp"
    publishDir "${params.outdir}/allsamp", mode: 'copy', overwrite: true

    input:
      path(data) from params.data
      path(genes) from params.genes
      val(permutation) from params.permutation
      val(metric) from params.metric
      val(ivar) from params.ivar
      val(norm) from params.norm
       
    output:
      file("enzyme_*") into ch_out1
    
    when:
      if (params.check){
        if (!file("${params.outdir}/allsamp/enzyme_results.Rdata").exists() && params.allsample && params.tis==false) {true} else {false}
      }else{
        if (params.allsample && params.tis==false) {true} else {false}
      }

    script:
      """
      # print cluster working node's name
      hostname
      
      # print R version
      R --version
      Rscript --version
      
      # run proportionality analysis
      Rscript ${baseDir}/bin/propr/compute_proportionality_alldata_fromgenes.R \
                                              --data ${data} \
                                              --genes ${genes} \
                                              --outdir . \
                                              --permutation ${permutation} \
                                              --metric ${metric} \
                                              --ivar ${ivar} \
                                              --norm ${norm}
      """
}

/* 2. Proportionality-enzyme with all samples from a specific tissue */
process propr_allsamp_t {

    memory '24.GB'

    tag "allsamp"
    publishDir "${params.outdir}/allsamp", mode: 'copy', overwrite: true

    input:
      path(data) from params.data
      path(genes) from params.genes
      path(tissue) from params.tissue
      val(permutation) from params.permutation
      val(metric) from params.metric
      val(ivar) from params.ivar
      val(norm) from params.norm
       
    output:
      file("enzyme_*") into ch_out2
    
    when:
      if (params.check){
        if (!file("${params.outdir}/allsamp/enzyme_results.Rdata").exists() && params.allsample && params.tis) {true} else {false}
      }else{
        if (params.allsample && params.tis) {true} else {false}
      }

    script:
      """
      # print cluster working node's name
      hostname
      
      # print R version
      R --version
      Rscript --version
      
      # run proportionality analysis
      Rscript ${baseDir}/bin/propr/compute_proportionality_alldata_fromgenes.R \
                                              --data ${data} \
                                              --genes ${genes} \
                                              --outdir . \
                                              --permutation ${permutation} \
                                              --tissue ${tissue} \
                                              --metric ${metric} \
                                              --ivar ${ivar} \
                                              --norm ${norm}
      """
}

/* 3. Proportionality-enzyme with n randomly selected samples */
process propr_nsamp {

    // memory { nsamp < 1000 ? '16.GB' : 
    //          nsamp < 5000 ? '48.GB' : 
    //          '64.GB' }
    memory { nsamp < 1000 ? '24.GB' : 
             nsamp < 5000 ? '64.GB' :
             '84.GB'} 

    tag "${nsamp}_${nseed}"
    publishDir "${params.outdir}/${nseed}_singularity/n${nsamp}", mode: 'copy', overwrite: true

    input:
      path(data) from params.data
      path(genes) from params.genes
      val(permutation) from params.permutation
      val(metric) from params.metric
      val(ivar) from params.ivar
      val(norm) from params.norm
      tuple val(nseed), file(seed), val(nsamp) from ch_seeds1
       
    output:
      file("enzyme_*") into ch_out3
    
    when:
      if (params.check){
        if (!file("${params.outdir}/${nseed}_singularity/n${nsamp}/enzyme_results.Rdata").exists() && params.usensample && params.tis==false) {true} else {false}
      }else{
        if (params.usensample && params.tis==false) {true} else {false}
      }

    script:
      """
      # print cluster working node's name
      hostname
      
      # print R version
      R --version
      Rscript --version
      
      # run proportionality analysis
      Rscript ${baseDir}/bin/propr/compute_proportionality_alldata_fromgenes.R \
                                              --data ${data} \
                                              --genes ${genes} \
                                              --outdir . \
                                              --permutation ${permutation} \
                                              --nsamp ${nsamp} \
                                              --sampseed ${seed} \
                                              --metric ${metric} \
                                              --ivar ${ivar} \
                                              --norm ${norm}
      """
}

/* 4. Proportionality-enzyme with n randomly selected samples from a specific tissue */
process propr_nsamp_t {

    memory { nsamp < 1000 ? '16.GB' : 
             nsamp < 5000 ? '48.GB' : 
             '64.GB' }

    tag "${nsamp}_${nseed}"
    publishDir "${params.outdir}/${nseed}_singularity/n${nsamp}", mode: 'copy', overwrite: true

    input:
      path(data) from params.data
      path(genes) from params.genes
      path(tissue) from params.tissue
      val(permutation) from params.permutation
      val(metric) from params.metric
      val(ivar) from params.ivar
      val(norm) from params.norm
      tuple val(nseed), file(seed), val(nsamp) from ch_seeds2
       
    output:
      file("enzyme_*") into ch_out4
    
    when:
      if (params.check){
        if (!file("${params.outdir}/${nseed}_singularity/n${nsamp}/enzyme_results.Rdata").exists() && params.usensample && params.tis) {true} else {false}
      }else{
        if (params.usensample && params.tis) {true} else {false}
      }

    script:
      """
      # print cluster working node's name
      hostname
      
      # print R version
      R --version
      Rscript --version
      
      # run proportionality analysis
      Rscript ${baseDir}/bin/propr/compute_proportionality_alldata_fromgenes.R \
                                              --data ${data} \
                                              --genes ${genes} \
                                              --outdir . \
                                              --permutation ${permutation} \
                                              --nsamp ${nsamp} \
                                              --sampseed ${seed} \
                                              --tissue ${tissue} \
                                              --metric ${metric} \
                                              --ivar ${ivar} \
                                              --norm ${norm}
      """
}


workflow.onComplete {
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}"
}
