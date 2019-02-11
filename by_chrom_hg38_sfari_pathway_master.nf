#!/usr/bin/env nextflow

// Using bed files for each KEGG pathway do all downstream work

// set up channel to read in bed files
Channel.fromPath('*.bed').set{ bed_files }

// split bed file by chroms
process split_chroms {
  input:
    file bed from bed_files
  
  output:
    file("${bed.baseName}_*.bed") into bed_by_chroms

  script:
    """
    python $GIT/py_tools/split_bed_by_chrom.py $bed ${bed.baseName}
    """
}

// extract from original bcf
process extract {
  module 'bcftools/1.7'
  cpus 1
  errorStrategy { task.attempt == 1 ? 'retry' : 'finish' }
  tag "$chrom_bed"
  input:
    file chrom_bed from bed_by_chroms

  output:
    file("${chrom_bed.baseName}.vcf") into chrom_vcf_files

  script:
    """
    bcftools view /scratch/ucgd/lustre/work/u0806040/data/vcf_files/P231_hg38.bcf -R $chrom_bed \
        -O v -o ${chrom_bed.baseName}.vcf
    """
}

// VEP annotation
process vep {
  module 'ucgd'
  module 'vep/91.3'
  cpus 16
  errorStrategy { task.attempt == 1 ? 'retry' : 'finish' }
  memory 120.GB
  tag "${chrom_vcf.baseName}_vep.vcf"
  
  input:
    file chrom_vcf from chrom_vcf_files

  output:
    file("${chrom_vcf.baseName}_vep.vcf") into vep_chrom_vcf_files

  script:
    """
    vep --format vcf --input_file $chrom_vcf --vcf --output_file ${chrom_vcf.baseName}_vep.vcf \
        --cache --dir ~/VEP_Cache_2.0 --assembly GRCh38 --offline --fork 16 --everything --use_given_ref \
        --fasta /scratch/ucgd/lustre/work/u0806040/data/reference_genomes/hg38.fa
    """
}

// gnomAD annotation with vcfanno
process gnomad {
  module 'vcfanno/0.2.4'
  cpus 6
  memory { 30.GB * task.attempt }
  ERrorStrategy { task.attempt == 1 ? 'retry' : 'finish' }
  tag "${chrom_vep_vcf.baseName}_gnomad.vcf"
  publishDir "$baseDir/raw_vcf_files", mode:"copy"

  input:
    file vep_chrom_vcf from vep_chrom_vcf_files
    file("/scratch/ucgd/lustre/work/u0806040/data/gnomAD_hg37.conf")
 
  output:
    file("${vep_chrom_vcf.baseName}_gnomad.vcf") into gnomad_vcf_files
 
  script:
    """
    vcfanno -p 6 /scratch/ucgd/lustre/work/u0806040/data/gnomAD_hg38.conf $vep_chrom_vcf > ${vep_chrom_vcf.baseName}_gnomad.vcf
    """
}

// combine all chroms into single vcf
process combine_chroms {
  module 'ucgd'
  
  input:
    file x from gnomad_vcf_files.collect()

  output:
    file("${bed.baseName}_annotated_combined.vcf") into combined_vcf_files

  script:
    """
    bcftools concat $x -O v -o ${bed.baseName}_annotated_combined.vcf 
    """
}

// generate variant count table
process generate_counts {
  module 'htslib/1.7'
  publishDir "$baseDir/variant_counts_files", mode:"copy"
  errorStrategy { task.attempt == 1 ? 'retry' : 'finish' }
  cpus 3
  memory 120.GB
  tag "${pathway_vep_sorted_gnomad_vcf.baseName}_variant_table.txt"

  input:
    file pathway_vep_sorted_gnomad_vcf from combined_vcf_files 
    file("/scratch/ucgd/lustre/work/u0806040/data/sfari_gene_score_dict.txt")
 
  output:
    file("${pathway_vep_sorted_gnomad_vcf.baseName}_variant_table.txt") into variant_tables
    file("${pathway_vep_sorted_gnomad_vcf.baseName}_variant_metadata.txt")
 
  script:
  """
  /scratch/ucgd/lustre/work/u0806040/src/variant_table $pathway_vep_sorted_gnomad_vcf
 """
}

// calculate the bins based on quartiles of variant counts, after applying filters
process generate_bins {
  publishDir "$baseDir/clinco_input", mode:"copy"
  tag "${variant_table.baseName}_bins"
  cpus 2
  memory 120.GB
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
