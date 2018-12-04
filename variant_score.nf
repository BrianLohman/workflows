#!/usr/bin/env nextflow

/* generated VARIANT score (NO WEIGHTING SCHEME), 
 * assign group membership based on score, 
 * test for differences among groups
 * NO split half procedure
 *
 * USAGE: nextflow run weighted_score.nf --in [target vcf]
 * REQUIRES:
 *  count_variant_effects (requries htslib)
 *    sfari_gene_score_dict.txt (dictonary of sfari gene sores for each gene)
 *  assign_groups_variant_score.R
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
    file("${vcf}_variant_counts_result.txt") into filestream

  script:
    """
    /scratch/ucgd/lustre/work/u0806040/KEGG_pathway_gene_lists/SFARI_genes/variant-score/no_filter_count_variant_effects $vcf
    """
}

// assign group membership based on score
process define_groups {
  module 'ucgd'

  input:
    file x from filestream
    file("./assign_groups_variant_score.R")

  output:
    file("variant_scores.txt") into testing

  script:
    """
    R CMD BATCH /scratch/ucgd/lustre/work/u0806040/KEGG_pathway_gene_lists/SFARI_genes/variant-score/assign_groups_variant_score.R
    """
}

// test for differences among groups
process testing {
  publishDir "$baseDir/", mode: "copy"

  input:
    file("variant_scores.txt") from testing

  output:
    file("variant_model_testing_results.txt") 

  script:
  """
  python /scratch/ucgd/lustre/work/u0806040/clinco/clinical-components/clinco/clincorr.py variant_scores.txt --column group > variant_model_testing_results.txt
  """
}
