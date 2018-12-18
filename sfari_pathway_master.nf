#!/usr/bin/env nexflow

// Using bed files for each KEGG pathway do all downstream work

// set up channel to read in bed files
Channel.fromPath('*.bed').set{ bed_files }

// extract from original bcf
process extract {
  module 'bcftools/1.7'

  input:
    file bed from bed_files

  output:
    file("${bed}_pathway.vcf") into pathway_vcf_files

  script:
    """
    bcftools view /scratch/ucgd/lustre/work/u0806040/data/vcf_files/Utah_519.bcf -R $bed -O v -o ${bed}_pathway.vcf
    """
}

// VEP annotation
process vep {
  module 'ucgd'
  module 'vep/91.3'

  input:
    file pathway_vcf from pathway_vcf_files

  output:
    file("${pathway_vcf}_vep.vcf") into pathway_vep_vcf_files

  script:
    """
    vep --format vcf --input_file $pathway_vcf --vcf --output_file ${pathway_vcf}_vep.vcf \
        --cache --dir ~/VEP_Cache_2.0 --assembly GRCh37 --offline --fork 24 --everything --use_given_ref \
        --fasta /scratch/ucgd/lustre/work/u0806040/data/reference_genomes/g1k_v37_decoy.fa
    """
}

// sorting
process sort {
  module 'bcftools/1.7'

  input:
    file pathway_vep_vcf from pathway_vep_vcf_files

  output:
    file("${pathway_vep_vcf}_sorted.vcf") into pathway_vep_sorted_vcf_files

  script:
    """
    bcftools sort $pathway_vep_vcf -m 10G -T $TMPDIR -o ${pathway_vep_vcf}_sorted.vcf -O v
    """
}

// gnomAD annotation with vcfanno
process gnomad {
  module 'vcfanno/0.2.4'
 
  input:
    file pathway_vep_sorted_vcf from pathway_vep_sorted_vcf_files
    file("/scratch/ucgd/lustre/work/u0806040/data/gnomAD_hg37.conf")
 
  output:
    file("${pathway_vep_sorted_vcf}_gnomad.vcf") into pathway_vep_sorted_gnomad_vcf_files
 
  script:
    """
    vcfanno /scratch/ucgd/lustre/work/u0806040/data/gnomAD_hg37.conf  $pathway_vep_sorted_vcf > ${pathway_vep_sorted_vcf}_gnomad.vcf
    """
}
 
// generate variant count table
process generate_counts {
  module 'htslib/1.7'
  publishDir "$baseDir/", mode = "copy"
 
  input:
    file pathway_vep_sorted_gnomad_vcf from pathway_vep_sorted_gnomad_vcf_files 
    file("/scratch/ucgd/lustre/work/u0806040/data/sfari_gene_score_dict.txt")
 
  output:
    file("${pathway_vep_sorted_gnomad_vcf}_variant_table.txt")
    file("${pathway_vep_sorted_gnomad_vcf}_variant_metadata.txt")
 
  script:
  """
  /scratch/ucgd/lustre/work/u0806040/data/sfari_genes/split_by_pathway/neural/variant_table $pathway_vep_sorted_gnomad_vcf
 """
}
