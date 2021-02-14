#!/usr/bin/env nextflow


log.info """\
         COMPARE PROPORTIONALITY VS PEARSON
         ======================================="
         Input dataset ID                               : ${params.id}
         Input dataset file (DATASET)                   : ${params.data}
         Input genes (FILE)                             : ${params.genes}
         Association+normalization+transformation combs : ${params.combs}
         Number permutation - propr (int)               : ${params.propr_permutation}
         Number permutation - graflex (int)             : ${params.graflex_permutation}
         Cutoff of propr/pearson                        : ${params.cutoff}
         Filter pathway-genes                           : ${params.filter}
         Main output directory (DIRECTORY)              : ${params.outdir}
         Trace report (DIRECTORY)                       : ${params.tracedir}
         If check file already exists (bool)            : ${params.check}
         """
         .stripIndent()



// combs channel : [comb, method, norm, ivar]
Channel
  .fromList(params.combs)
  .map { item -> [ item, item.tokenize('_')[0], item.tokenize('_')[1], item.tokenize('_')[2] ]}
  .set { ch_combs }


process propr {

    memory '84.GB'
    
    tag "${comb}"
    publishDir "${params.outdir}/network/${comb}", mode: 'copy', overwrite: true

    input:
      path(data) from params.data
      path(genes) from params.genes
      val(permutation) from params.propr_permutation
      tuple val(comb), val(method), val(norm), val(ivar) from ch_combs
       
    output:
      tuple val(comb), file("enzyme_*") into ch_out1
      tuple val(comb), file("enzyme_results.csv") into ch_proprOut
    
    when:
      if (params.check){
        if (!file("${params.outdir}/network/${comb}/enzyme_results.csv").exists()) {true} else {false}
      }else{
        true
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


// concat propr output with existing files
Channel
  .fromPath("${params.outdir}/network/*/enzyme_results.csv")
  .map { item -> [ item.getParent().baseName, item ] }  // [comb, propr result file]
  .concat( ch_proprOut )
  .unique()
  .set { ch_toGraflex }

// cutoff channel 
Channel
  .from(params.cutoff)
  .map { item -> [ item, (item*100).toInteger() ]}
  .set { ch_cutoff }


process graflex {

    memory '64.GB'

    tag "${cutoff}"
    publishDir "${params.output}/graflex${filter}/${comb}", mode: 'copy', overwrite: true

    input:
      tuple val(comb), file(input) from ch_toGraflex
      path(genes) from params.genes
      val(permutation) from params.graflex_permutation
      tuple val(cutoff), val(cutoff2) from ch_cutoff
      val(filter) from params.filter
    
    output:
      file("graflex_enzyme_c${cutoff2}_p${permutation}.tsv") into ch_out2

    when:
      if (params.check){
        if (!file("${params.outdir}/graflex${filter}/${comb}/graflex_enzyme_c${cutoff2}_p${permutation}.tsv").exists()) {true} else {false}
      }else{
        true
      }

    script:
      """
      
      # print cluster working node's name
      hostname
      
      # print R version
      R --version
      Rscript --version
      
      # run proportionality analysis
      Rscript ${baseDir}/bin/propr/calculateOR_enzyme.R \
                                           -i ${input} \
                                           -g ${genes} \
                                           -c ${cutoff} \
                                           -p ${permutation} \
                                           -o . \
                                           ${filter}
      
      """
}


workflow.onComplete {
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}"
}
