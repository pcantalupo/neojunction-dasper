#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

  // Pipeline configuration summary
  println "\n=== Pipeline Configuration Summary ==="
  println "Command Line       : ${workflow.commandLine}"
  println "Samplesheet        : ${params.samplesheet}"
  println "GTF File           : ${params.gtf}"
  println "Output Dir         : ${params.outdir}"
  println "Dasper Outdir      : ${params.dasper_outdir}"
  println "Min Intron         : ${params.min_intron}"
  println "Max Intron         : ${params.max_intron}"
  println "Scale Factor       : ${params.scale_factor}"
  println "Counts Column      : ${params.counts_col}"
  println "Project Dir        : ${workflow.projectDir}"
  println "Work Dir           : ${workflow.workDir}"
  println "PipelineInfo Dir   : ${params.pipeline_info}"
  println "Container Engine   : ${workflow.containerEngine ?: 'NONE'}"
  println "=======================================\n\n"

  // Samplesheet
  if (params.samplesheet) { file(params.samplesheet, checkIfExists: true) }
  else { exit 1, "Error: please specify --samplesheet with the path to your samplesheet file" }
  ch_samples = Channel.fromPath(params.samplesheet)
      .splitCsv(header: true, sep: ',')
      .map { row -> tuple(row.sample_id, file(row.sj_file, checkIfExists: true), row.library_size.toInteger()) }

  // GTF file
  if (params.gtf) { file(params.gtf, checkIfExists: true) }
  else { exit 1, "Error: please specify --gtf with the path to your gtf file" }
  ch_gtf = Channel.fromPath(params.gtf)


  // Run dasper on each sample
  DASPER(ch_samples.combine(ch_gtf))    // need 'combine' operator to attach it to each sample
  //DASPER.out.dasper.view()

  // Filter the dasper output
  FILTER_DASPER(DASPER.out.dasper)
  //FILTER_DASPER.out.filtered_dasper.view()

  /* Convert FILTER_DASPER output channel:
  [S18, /path/to/S18.dasper.tsv.filtered.tsv]
  [S19, /path/to/S19.dasper.tsv.filtered.tsv]
  Into a single tuple with 2 elements: first is a list of samples and second is a list of the corresponding files.
  [ [S18, S19], [/path/to/S18.dasper.tsv.filtered.tsv, /path/to/S19.dasper.tsv.filtered.tsv] ]
  */
  combined_filter_dasper_output_ch = FILTER_DASPER.out.filtered_dasper
    .collect(flat: false) // List of [sample_id, path] tuples
    .map { it.transpose() }

  // Run neojunction on all filtered dasper files
  NEOJUNCTION(combined_filter_dasper_output_ch)
  //NEOJUNCTION.out.results.view()

}

// Process: Run dasper.R on each SJ.out.tab file
process DASPER {
  tag "${sample_id}"

  container "${ workflow.containerEngine == 'singularity' ?
                  'docker://virushunter/dasper:v1.0.0' :
                  'virushunter/dasper:v1.0.0' }"

  publishDir "${params.outdir}/${params.dasper_outdir}", pattern: "*.dasper.tsv", mode: 'copy'
  publishDir "${params.outdir}/${params.dasper_outdir}/outs", pattern: "*.out", mode: 'copy'

  input:
  tuple val (sample_id), path (sj_file), val (library_size), path(gtf)

  output:
  tuple val (sample_id), path ("${sample_id}.dasper.tsv"), val (library_size),   emit: dasper
  path ("*.out"),                                                                emit: logs


  script:
  """
  dasper.R --prefix ${sample_id} --sj ${sj_file} --gtf ${gtf} --outdir "." &> ${sample_id}_dasper.out
  """
}

// Process: Run filter_dasper.R on each dasper output
// Creates a temporary librarysizefile per sample
process FILTER_DASPER {
  tag "${sample_id}"

  container "${ workflow.containerEngine == 'singularity' ?
                  'docker://virushunter/dasper:v1.0.0' :
                  'virushunter/dasper:v1.0.0' }"

  publishDir "${params.outdir}/${params.dasper_outdir}", pattern: "*.filtered.tsv", mode: 'copy'
  publishDir "${params.outdir}/${params.dasper_outdir}/outs", pattern: "*.out", mode: 'copy'

  input:
  tuple val (sample_id), path (dasper_file), val (library_size)

  output:
  tuple val (sample_id), path ("${sample_id}.dasper.tsv.filtered.tsv"), emit: filtered_dasper
  path ("*.out")                                             , emit: logs

  script:
  """
  filter_dasper.R --dasperannotfile ${dasper_file} --librarysize ${library_size} \
                  --min ${params.min_intron} --max ${params.max_intron} --sf ${params.scale_factor} \
                  --outdir "." &> ${sample_id}_filterdasper.out
  """
}

// Process: Run neojunction.R on all filtered.tsv files
process NEOJUNCTION {

  publishDir "${params.outdir}/", mode: 'copy'

  container "${ workflow.containerEngine == 'singularity' ?
                  'docker://virushunter/dasper:v1.0.0' :
                  'virushunter/dasper:v1.0.0' }"

  input:
  tuple val(sample_ids), path(filtered_files) // e.g. [ [S18, S19], [/path/to/S18.dasper.tsv.filtered.tsv, /path/to/S19.dasper.tsv.filtered.tsv] ]

  output:
  path ("neojunction_results.tsv"),  emit: results
  path ("*.out"),                    emit: logs

  script:
  /* Create a metadata string like "sample_id,filtered_file\nS19,S19.dasper.tsv.filtered.tsv\nS18,S18.dasper.tsv.filtered.tsv" so bash can output it as a file called metadata.csv for use by neojunction.R script
  */
  def metadataRows = sample_ids.withIndex().collect { sample, i -> "${sample},${filtered_files[i]}" }
  def metadataString = (["sample_id,filtered_file"] + metadataRows).join("\\n")
  """
  echo -e "${metadataString}" > metadata.csv
  neojunction.R --metadata metadata.csv --counts_col ${params.counts_col} &> neojunction.out
  """
}

