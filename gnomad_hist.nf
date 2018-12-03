#!/Usr/bin/env nextflow

/* USAGE: nextflow run gnomad_hist --in [target vcf file] 
 * make histograms of the gnomAD allele frequency for given vcf
 */

//vcf = file(params.in)

file_handle = Channel
  .fromPath(params.in)
  .map { file -> tuple(file.baseName, file) }


// run gnomAD_AF_hist vi nim
process collect_gnomad_af {
  module 'htslib/1.7'

  input:
    set handle, file(vcf) from file_handle

  output:
    file("${vcf}_high_impact_gnomAD_AF.txt") into counts1
    file("${vcf}_med_impact_gnomAD_AF.txt") into counts2

  script:
    """
    $SCRATCH/nextflow_tests/gnomAD_AF_hist $vcf
    """
}

// plot results with R
process plot_hist {
  publishDir '/scratch/ucgd/lustre/work/u0806040/nextflow_tests', mode: "copy"
  
  input:
    file x from counts1
    file y from counts2
  
  output:
    file("gnomAD_hist.pdf")

  script:
    """
    R CMD BATCH $SCRATCH/nextflow_tests/plot_gnomad_hist.R
    """
}
