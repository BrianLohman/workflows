#!/usr/bin/env nextflow

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
    file("${bed.baseName}_pathway.vcf") into pathway_vcf_files

  script:
    """
    bcftools view /scratch/ucgd/lustre/work/u0806040/data/vcf_files/Utah_519.bcf -R $bed -O v -o ${bed.baseName}_pathway.vcf
    """
}

// VEP annotation
process vep {
  module 'ucgd'
  module 'vep/91.3'
  cpus 16
  memory 120.GB
  tag "${pathway_vcf.baseName}_vep.vcf"
  
  input:
    file pathway_vcf from pathway_vcf_files

  output:
    file("${pathway_vcf.baseName}_vep.vcf") into pathway_vep_vcf_files

  script:
    """
    vep --format vcf --input_file $pathway_vcf --vcf --output_file ${pathway_vcf.baseName}_vep.vcf \
        --cache --dir ~/VEP_Cache_2.0 --assembly GRCh37 --offline --fork 16 --everything --use_given_ref \
        --fasta /scratch/ucgd/lustre/work/u0806040/data/reference_genomes/g1k_v37_decoy.fa
    """
}

// sorting
process sort {
  module 'bcftools/1.7'
  cpus 1
  memory { 12.GB * task.attempt }
  errorStrategy { task.attempt == 1 ? 'retry' : 'finish' }
  tag "${pathway_vep_vcf.baseName}_sorted.vcf"
  input:
    file pathway_vep_vcf from pathway_vep_vcf_files

  output:
    file("${pathway_vep_vcf.baseName}_sorted.vcf") into pathway_vep_sorted_vcf_files

  script:
    """
    bcftools sort $pathway_vep_vcf -m 5G -T $TMPDIR -o ${pathway_vep_vcf.baseName}_sorted.vcf -O v
    """
}

// gnomAD annotation with vcfanno
process gnomad {
  module 'vcfanno/0.2.4'
  cpus 6
  memory { 20.GB * task.attempt }
  errorStrategy { task.attempt == 1 ? 'retry' : 'finish' }
  tag "${pathway_vep_sorted_vcf.baseName}_gnomad.vcf"
  publishDir "$baseDir/raw_vcf_files", mode:"copy"

  input:
    file pathway_vep_sorted_vcf from pathway_vep_sorted_vcf_files
    file("/scratch/ucgd/lustre/work/u0806040/data/gnomAD_hg37.conf")
 
  output:
    file("${pathway_vep_sorted_vcf.baseName}_gnomad.vcf") into pathway_vep_sorted_gnomad_vcf_files
 
  script:
    """
    vcfanno -p 6 /scratch/ucgd/lustre/work/u0806040/data/gnomAD_hg37.conf $pathway_vep_sorted_vcf > ${pathway_vep_sorted_vcf.baseName}_gnomad.vcf
    """
}
 
// generate variant count table
process generate_counts {
  module 'htslib/1.7'
  publishDir "$baseDir/variant_counts_files", mode:"copy"
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

// calculate the bins based on quartiles of variant counts, after applying filters
process generate_bins {
  publishDir "$baseDir/clinco_input", mode:"copy"
  tag "${variant_table.baseName}_bins"
  cpus 2
  memory { 16.GB * task.attempt }
  errorStrategy { task.attempt == 1 ? 'retry' : 'finish' }

  input:
    file variant_table from variant_tables

  output:
    file("${variant_table.baseName}_bins.txt") into bins
    file("dedup_$variant_table") into dedup_variant_tables 

  script:
  """
  python /scratch/ucgd/lustre/work/u0806040/GitHub/py_tools/bins_by_mutation_load.py -v ${variant_table} -o ${variant_table.baseName}_bins.txt 
  """
}

// run clinco to test for associations
process clinco {
  publishDir "$baseDir/clinco_results", mode:"copy"
  tag "${bins_table.baseName}_clinco"
  cpus 4
  memory { 12.GB * task.attempt }
  errorStrategy { task.attempt == 1 ? 'retry' : 'finish' }

  input:
    file bins_table from bins

  output:
    file("${bins_table.baseName}_clinco_results.txt")

  script:
  """
  python /scratch/ucgd/lustre/work/u0806040/clinco/clinical-components/clinco/clincorr.py ${bins_table} --column group > ${bins_table.baseName}_clinco_results.txt
  """
}

// make table of individuals and variant impact and plot
process plotting {
  publishDir "$baseDir/plots", mode:"copy"
  tag "${dedup_table}_plotting"
  memory { 30.GB * task.attempt }
  errorStrategy { task.attempt == 1 ? 'retry' : 'finish' }

  input:
    file dedup_table from dedup_variant_tables

  output:
    file("individuals_by_impact_${dedup_table}")
    file("${dedup_table}_variant_plot.html")

  script:
  """
  python /scratch/ucgd/lustre/work/u0806040/GitHub/py_tools/impact_table_summary.py ${dedup_table}
  """
}
