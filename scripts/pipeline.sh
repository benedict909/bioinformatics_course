#!/bin/bash

workdir=~/ngs_course/dnaseq_pipeline/

# 1. Perform trimming with Trimmomatic on raw sequencing data
mkdir $workdir/data/trimmed_fastq
trimmomatic PE  \
  -threads 4 \
  -phred33 \
  $1 $2 \
  -baseout $workdir/data/trimmed_fastq/WES01_chr22m_trimmed_R \
  ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/NexteraPE-PE.fa:2:30:10 \
  TRAILING:25 MINLEN:50

# 2. FastQC Analysis
fastqc -t 4 $workdir/data/trimmed_fastq/WES01_chr22m_trimmed_R_1P \
  $workdir/data/trimmed_fastq/WES01_chr22m_trimmed_R_2P
mkdir $workdir/results/fastqc_trimmed_reads
mv $workdir/data/trimmed_fastq/*fastqc* $workdir/results/fastqc_trimmed_reads/

# 3.1 Alignment
alignment_dir=$workdir/data/aligned_data
mkdir $alignment_dir

bwa mem -t 4 -v 1 -R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.WES01\tSM:WES01\tPL:ILLUMINA\tLB:nextera-wes01-blood\tDT:2017-02-23\tPU:HWI-D00119' \
  -I 250,50  $workdir/data/reference/hg19.fa.gz $workdir/data/trimmed_fastq/WES01_chr22m_trimmed_R_1P \
  $workdir/data/trimmed_fastq/WES01_chr22m_trimmed_R_2P > $workdir/data/aligned_data/WES01_chr22m.sam

# 3.2 Convert, process and index SAM file
samtools view -h -b $alignment_dir/WES01_chr22m.sam > $alignment_dir/WES01_chr22m.bam
samtools sort $alignment_dir/WES01_chr22m.bam > $alignment_dir/WES01_chr22m_sorted.bam
samtools index $alignment_dir/WES01_chr22m_sorted.bam

# 3.3 Post Alignment QC and Filtering
picard MarkDuplicates I=$alignment_dir/WES01_chr22m_sorted.bam O=$alignment_dir/WES01_chr22m_sorted_marked.bam \
  M=$alignment_dir/marked_dup_metrics.txt
samtools index $alignment_dir/WES01_chr22m_sorted_marked.bam

# 3.4 Filter BAM based on mapping quality and bitwise flags using samtools
samtools view -F 1796  -q 20 -o $alignment_dir/WES01_chr22m_sorted_marked_filtered.bam \
  $alignment_dir/WES01_chr22m_sorted_marked.bam

samtools index $alignment_dir/WES01_chr22m_sorted_marked_filtered.bam

# 3.5 Alignment Statistics
mkdir $alignment_dir/alignment_stats
stats_dir=$alignment_dir/alignment_stats

samtools flagstats $alignment_dir/WES01_chr22m_sorted_marked_filtered.bam \
  > $stats_dir/flagstats_output.txt
samtools view -h $alignment_dir/WES01_chr22m_sorted_marked_filtered.bam > $alignment_dir/WES01_chr22m_sorted_marked_filtered.sam
samtools idxstats  $alignment_dir/WES01_chr22m_sorted_marked_filtered.bam \
  > $stats_dir/idxstats_output.txt

picard CollectInsertSizeMetrics -H $stats_dir/CollectInsertSizeMetrics_histogram.pdf -I $alignment_dir/WES01_chr22m_sorted_marked_filtered.bam -O $stats_dir/CollectInsertSizeMetrics_output.txt

coverageBed -d -a $workdir/data/chr22.genes.hg19.bed -b $alignment_dir/WES01_chr22m_sorted_marked_filtered.bam \
  > $stats_dir/coverageBed_output.txt


# 4.1 Variant calling with Freebayes 
freebayes --bam $alignment_dir/WES01_chr22m_sorted_marked_filtered.bam \
  --fasta-reference $workdir/data/reference/hg19.fa --vcf $workdir/results/WES01_chr22m.vcf
bgzip -f $workdir/results/WES01_chr22m.vcf
tabix -p vcf $workdir/results/WES01_chr22m.vcf.gz

# 4.2 Filtering the VCF
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
  $workdir/results/WES01_chr22m.vcf.gz > $workdir/results/WES01_chr22m_filtered.vcf

bedtools intersect -header -wa -a $workdir/results/WES01_chr22m_filtered.vcf -b $workdir/data/chr22.genes.hg19.bed \
  > $workdir/results/WES01_chr22m_filtered_chr22.vcf
bgzip -f $workdir/results/WES01_chr22m_filtered_chr22.vcf
tabix -p vcf $workdir/results/WES01_chr22m_filtered_chr22.vcf.gz 


# 5.1 Variant Annotation and Prioritisation
~/tools/annovar/convert2annovar.pl -format vcf4 $workdir/results/WES01_chr22m_filtered_chr22.vcf.gz  \
  > $workdir/results/WES01_chr22m_filtered_chr22.avinput

~/tools/annovar/table_annovar.pl $workdir/results/WES01_chr22m_filtered_chr22.avinput ~/tools/annovar/humandb/ -buildver hg19 \
  -out $workdir/results/WES01_chr22m_filtered_chr22 -remove \
  -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout