#!/usr/bin/env nextflow

// run CRAM files though bazam and kraken2 by sample ID

// sample IDs are one per line in text file
Channel
  file('params.input)
  .readlines()
  .set { samples_ch }

// run shell script with bazam and kraken2
process bazam_kraken {
  input:
  file(${sample}.final.cram)
  file(${sample}.final.cram.crai)

  output:
  file(${sample}_kraken2_report.txt)
  file(${sample}_unclassified_r1.fq)
  file(${sample}_unclassified_r2.fq)

  script:
  """
  bazam_kraken2.sh $sample
  """

}
