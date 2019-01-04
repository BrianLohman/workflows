#!/USr/bin/env nexflow

// Using bed files for each KEGG pathway do all downstream work

// set up channel to read in bed files
Channel.fromPath('*.bed').set{ bed_files }

// extract from original bcf
process extract {
  module 'bcftools/1.7'
  cpus 1
  tag "$bed"
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
  cpus 16
  memory 120.GB
  tag "${pathway_vcf}_vep.vcf"
  
  input:
    file pathway_vcf from pathway_vcf_files

  output:
    file("${pathway_vcf}_vep.vcf") into pathway_vep_vcf_files

  script:
    """
    vep --format vcf --input_file $pathway_vcf --vcf --output_file ${pathway_vcf}_vep.vcf \
        --cache --dir ~/VEP_Cache_2.0 --assembly GRCh37 --offline --fork 16 --everything --use_given_ref \
        --fasta /scratch/ucgd/lustre/work/u0806040/data/reference_genomes/g1k_v37_decoy.fa
    """
}

// sorting
process sort {
  module 'bcftools/1.7'
  cpus 1
  memory 12.GB
  tag "${pathway_vep_vcf}_sorted.vcf"
  input:
    file pathway_vep_vcf from pathway_vep_vcf_files

  output:
    file("${pathway_vep_vcf}_sorted.vcf") into pathway_vep_sorted_vcf_files

  script:
    """
    bcftools sort $pathway_vep_vcf -m 5G -T $TMPDIR -o ${pathway_vep_vcf}_sorted.vcf -O v
    """
}

// gnomAD annotation with vcfanno
process gnomad {
  module 'vcfanno/0.2.4'
  cpus 2
  memory { 12.GB * task.attempt }
  errorStrategy { task.attempt == 1 ? 'retry' : 'finish' }
  tag "${pathway_vep_sorted_vcf}_gnomad.vcf"

  input:
    file pathway_vep_sorted_vcf from pathway_vep_sorted_vcf_files
    file("/scratch/ucgd/lustre/work/u0806040/data/gnomAD_hg37.conf")
 
  output:
    file("${pathway_vep_sorted_vcf}_gnomad.vcf") into pathway_vep_sorted_gnomad_vcf_files
 
  script:
    """
    vcfanno /scratch/ucgd/lustre/work/u0806040/data/gnomAD_hg37.conf $pathway_vep_sorted_vcf > ${pathway_vep_sorted_vcf}_gnomad.vcf
    """
}
 
// generate variant count table
process generate_counts {
  module 'htslib/1.7'
  publishDir "$baseDir", mode:"copy"
  cpus 3
  memory 40.GB
  tag "${pathway_vep_sorted_gnomad_vcf.baseName}_variant_table.txt"
  input:
    file pathway_vep_sorted_gnomad_vcf from pathway_vep_sorted_gnomad_vcf_files 
    file("/scratch/ucgd/lustre/work/u0806040/data/sfari_gene_score_dict.txt")
 
  output:
    file("${pathway_vep_sorted_gnomad_vcf.baseName}_variant_table.txt") into variant_tables
    file("${pathway_vep_sorted_gnomad_vcf.baseName}_variant_metadata.txt")
 
  script:
  """
  /scratch/ucgd/lustre/work/u0806040/sfari_genes/split_by_pathway/neural/variant_table $pathway_vep_sorted_gnomad_vcf
 """
}

// remove duplicates
process dedup {
  publishDir "$baseDir", mode:"copy"
  tag "${variant_table.baseName}_dedup.txt"
  input:
    file variant_table from variant_tables

  output:
    file("${variant_table.baseName}_dedup.txt")

  script:
  """
  sort ${variant_table} -u | sort -k1,1V -k 2,2n > ${variant_table.baseName}_dedup.txt
  """
}
