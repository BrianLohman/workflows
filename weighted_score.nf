#!/usr/bin/env nextflow

//TODO: set up channel to handle the spilt half procedure

/* generated weighted score, 
 * assign group membership based on score, 
 * test for differences among groups
 * USAGE: nextflow run weighted_score.nf --in [target vcf]
 * REQUIRES:
 *  count_variant_effects (requries htslib)
 *    sfari_gene_score_dict.txt (dictonary of sfari gene sores for each gene)
 *  assign_groups_weighted_score.R
 *    simons_master.txt (text file of pedigree and PRS/traits)
 */

// name input vcf
vcf = file(params.in)

// generate weighted score with count_variant_effects (nim)
process generate_scores {
  module 'htslib/1.7'

  input:
    file(vcf)

  output:
    file("${vcf}_variant_scores_result.txt") into filestream

  script:
    """
    /scratch/ucgd/lustre/work/u0806040/KEGG_pathway_gene_lists/SFARI_genes/weighted-score/count_variant_effects $vcf
    """
}

// assign group membership based on score
process define_groups {
  module 'ucgd'

  input:
    file x from filestream
    file("./assign_groups_weighted_score.R")

  output:
    file("probands_half1_for_clinco.txt") into filestream_half1
    file("probands_half2_for_clinco.txt") into filestream_half2

  script:
    """
    R CMD BATCH /scratch/ucgd/lustre/work/u0806040/KEGG_pathway_gene_lists/SFARI_genes/weighted-score/assign_groups_weighted_score.R
    """
}

// test for differences among groups
process testing_half1 {
  publishDir "$baseDir/", mode: "copy"

  input:
    file("probands_half1_for_clinco.txt") from filestream_half1

  output:
    file("half1_model_testing_results.txt") 

  script:
  """
  python /scratch/ucgd/lustre/work/u0806040/clinco/clinical-components/clinco/clincorr.py probands_half1_for_clinco.txt --column group > half1_model_testing_results.txt
  """
}

process testing_half2 {
  publishDir "$baseDir/", mode: "copy"

  input:
    file("probands_half2_for_clinco.txt") from filestream_half2

  output:
    file("half2_model_testing_results.txt")

  script:
  """
  python /scratch/ucgd/lustre/work/u0806040/clinco/clinical-components/clinco/clincorr.py probands_half2_for_clinco.txt --column group > half2_model_testing_results.txt
  """
}
