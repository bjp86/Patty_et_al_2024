#!/bin/bash
#
#SBATCH -N 1
#SBATCH --job-name=CR_analysis
#SBATCH --output=CR_out.txt
#SBATCH --time=02-00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bjp86@pitt.edu
#SBATCH --cpus-per-task=4
#SBATCH --mem=60g
#SBATCH -p scavenger

#############################

module load gcc/8.2.0
module load samtools/1.14

for f in *.fastq.gz; do gunzip $f; done 
for f in *.fastq; do awk '{if(NR%4==1){print $1} else{print substr($1, 1, 25)}}' $f > ${f/_001.fastq/_trim25.fastq}; done

module load fastqc

fastqc -t 6 *trim25.fastq 

module load gcc/8.2.0
module load samtools
module load bowtie2/2.4.5 

for f in *R1_trim25.fastq; do bowtie2 -p 6 -I 10 -X 1000 -N 1 --very-sensitive -x /ix1/shainer/mouse/mm10/mm10 -1 $f -2 ${f/R1/R2} -S ${f/_R1_trim25.fastq/.sam}; done

module load picard

for f in *.sam; do java -Xmx45g -jar /ix1/shainer/picard-tools-2.5.0/picard.jar SortSam INPUT=$f OUTPUT=${f/.sam/_picard.bam} VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp SORT_ORDER=coordinate; java -Xmx45g -jar /ix1/shainer/picard-tools-2.5.0/picard.jar MarkDuplicates INPUT=${f/.sam/_picard.bam} OUTPUT=${f/.sam/_picard2.bam} VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp METRICS_FILE=dup.txt REMOVE_DUPLICATES=true; done

module purge
module load gcc/8.2.0
module load samtools/1.14


# filter low-quality reads (MAPQ < 10) and convert sam to bam

for f in *_picard2.bam; do samtools view -h -@ 5 -O SAM -Sq 10 -o ${f/_picard2.bam/_filtered.sam} $f; done
 
for f in *_filtered.sam; do awk ' $9 <= 120 && $9 >= 1 || $9 >= -120 && $9 <= -1 ' $f > ${f/.sam/.1_120.sam}; cp /ix1/shainer/Dave/Master_Reference_Files/mm10/bowtie2.mm10.header ${f/_filtered.sam/_filtered.1_120.header}; cat ${f/.sam/.1_120.sam} >> ${f/_filtered.sam/_filtered.1_120.header}; rm ${f/_filtered.sam/_filtered.1_120.sam}; mv ${f/_filtered.sam/_filtered.1_120.header} ${f/_filtered.sam/_filtered.1_120.sam}; samtools view -@ 5 -S -t /ix1/shainer/mouse/mm10/mm10.chrom.sizes -b -o ${f/_filtered.sam/.1_120.bam} ${f/_filtered.sam/_filtered.1_120.sam}; done

for f in *.1_120.sam; do samtools sort -@ 5 -O BAM -o ${f/_filtered.1_120.sam/1_120_sorted.bam} $f; samtools index -@ 23 ${f/_filtered.1_120.sam/1_120_sorted.bam}; done

for f in *_filtered.sam; do awk ' $9 <= 200 && $9 >= 100 || $9 >= -200 && $9 <= -100 ' $f > ${f/.sam/.100_200.sam}; cp /ix1/shainer/Dave/Master_Reference_Files/mm10/bowtie2.mm10.header ${f/_filtered.sam/_filtered.100_200.header}; cat ${f/.sam/.100_200.sam} >> ${f/_filtered.sam/_filtered.100_200.header}; rm ${f/_filtered.sam/_filtered.100_200.sam}; mv ${f/_filtered.sam/_filtered.100_200.header} ${f/_filtered.sam/_filtered.100_200.sam}; samtools view -@ 5 -S -t /ix1/shainer/mouse/mm10/mm10.chrom.sizes -b -o ${f/_filtered.sam/.100_200.bam} ${f/_filtered.sam/_filtered.100_200.sam}; done

for f in *.100_200.sam; do samtools sort -@ 5 -O BAM -o ${f/_filtered.100_200.sam/100_200_sorted.bam} $f; samtools index -@ 23 ${f/_filtered.100_200.sam/100_200_sorted.bam}; done

for f in *_filtered.sam; do awk ' $9 <= 500 && $9 >= 130 || $9 >= -500 && $9 <= -130 ' $f > ${f/.sam/.130_500.sam}; cp /ix1/shainer/Dave/Master_Reference_Files/mm10/bowtie2.mm10.header ${f/_filtered.sam/_filtered.130_500.header}; cat ${f/.sam/.130_500.sam} >> ${f/_filtered.sam/_filtered.130_500.header}; rm ${f/_filtered.sam/_filtered.130_500.sam}; mv ${f/_filtered.sam/_filtered.130_500.header} ${f/_filtered.sam/_filtered.130_500.sam}; samtools view -@ 5 -S -t /ix1/shainer/mouse/mm10/mm10.chrom.sizes -b -o ${f/_filtered.sam/.130_500.bam} ${f/_filtered.sam/_filtered.130_500.sam}; done

for f in *.130_500.sam; do samtools sort -@ 5 -O BAM -o ${f/_filtered.130_500.sam/130_500_sorted.bam} $f; samtools index -@ 23 ${f/_filtered.130_500.sam/130_500_sorted.bam}; done

#Make bigwigs for viewing browser tracks
# module purge
module load gcc/8.2.0
module load deeptools

for f in *_sorted.bam; do bamCoverage -b $f -bs 5 --maxFragmentLength 1000 --smoothLength 20 --normalizeUsing RPGC --effectiveGenomeSize 2407883318 -e -p "max" -of bigwig -o ${f/.bam/.bw}; done
for f in *.bw; do computeMatrix reference-point -R mm10_refseq_select.bed -S $f -b 2000 -a 2000 -bs 20 -o ${f/.bw/.mat}; done
for f in *.mat; do plotHeatmap -m $f -o ${f/.bw/.mat}_TSSHeatmap.png --sortRegions descend --colorMap Reds --zMin 0.0 --zMax 5.0 --xAxisLabel Distance_From_Peak ; done

for f in *.bw; do computeMatrix scale-regions -R mm10_refseq_select.bed -S $f -m 10000 -b 1000 -a 1000 -bs 20 -o ${f/.bw/_SR.mat}; done
for f in *_SR.mat; do plotHeatmap -m $f -o ${f/.bw/.mat}_gene_body_Heatmap.png --sortRegions descend --colorMap Reds --zMin 0.0 --zMax 5.0 --xAxisLabel Gene_body ; done

