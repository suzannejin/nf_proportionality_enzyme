#!/usr/bin/env nextflow


log.info """\
         COMPARE PROPORTIONALITY VS PEARSON
         ======================================="
         Input dataset ID                               : ${params.id}
         Input dataset file (DATASET)                   : ${params.data}
         Input genes (FILE)                             : ${params.genes}
         Combs                                          : ${params.combs}
         Number permutation (int)                       : ${params.permutation}
         Output directory (DIRECTORY)                   : ${params.outdir}
         Trace report (DIRECTORY)                       : ${params.tracedir}
         If check file already exists (bool)            : ${params.check}
         """
         .stripIndent()

/*
 *  Manage channels
 */

// combs channel : [method, norm, ivar]
Channel
  .fromList(params.combs)
  .map { item -> [ item.tokenize(',')[0], item.tokenize(',')[1], item.tokenize(',')[2] ]}
  .into { ch_combs1; ch_combs2 }


process propr_rho {

    memory '84.GB'
    
    tag "${method}_${norm}_${ivar}"
    publishDir "${params.outdir}/${method}_${norm}_${ivar}", mode: 'copy', overwrite: true

    input:
      path(data) from params.data
      path(genes) from params.genes
      val(permutation) from params.permutation
      tuple val(method), val(norm), val(ivar) from ch_combs1
       
    output:
      file("enzyme_*") into ch_out1
    
    when:
      if (method == "rho")
      {
        if (params.check){
          if (!file("${params.outdir}/${method}_${norm}_${ivar}/enzyme_results.Rdata").exists()) {true} else {false}
        }else{
          true
        }
      }else{
        false
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
                                              --metric ${method} \
                                              --ivar ${ivar} \
                                              --norm ${norm}
      """
}


process propr_cor {

    memory '84.GB'

    tag "${method}_${norm}_${ivar}"
    publishDir "${params.outdir}/${method}_${norm}_${ivar}", mode: 'copy', overwrite: true

    input:
      path(data) from params.data
      path(genes) from params.genes
      val(permutation) from params.permutation
      tuple val(method), val(norm), val(ivar) from ch_combs2
       
    output:
      file("enzyme_*") into ch_out2
    
    when:
      if (method == "cor")
      {
        if (params.check){
          if (!file("${params.outdir}/${method}_${norm}_${ivar}/enzyme_results.Rdata").exists()) {true} else {false}
        }else{
          true
        }
      }else{
        false
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
                                              --metric ${method} \
                                              --ivar ${ivar} \
                                              --norm ${norm}
      """
}


workflow.onComplete {
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}"
}
