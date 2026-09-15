#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status

# --- 0. USER CONFIGURATION (ADJUST THESE!) ---

# Set the number of threads (cores)
THREADS=64

# --- Paths to genome and annotation files ---
GENOME_FASTA="02_genome/GRCm38.p4.genome.fa"
GENOME_GTF="02_genome/gencode.vM10.annotation.gtf"

# --- Library Strandness (CRITICAL for featureCounts) ---
# 0 = Unstranded
# 1 = Stranded
# 2 = Reversely Stranded (Most common for Illumina)
STRANDNESS=2

# --- END OF CONFIGURATION ---

echo "--- STARTING RNA-SEQ PIPELINE ---"
echo "Using $THREADS threads."

# --- 1. INITIAL QC (FASTQC) ---
echo "STEP 1: Running FastQC on raw data..."
fastqc -t $THREADS -o 03_qc_raw 01_data_raw/*.fastq.gz

# --- 2. TRIMMING (FASTP) ---
echo "STEP 2: Running fastp for adapter trimming..."
for R1 in 01_data_raw/*_R1_*.fastq.gz; do
    R2=$(echo $R1 | sed 's/_R1_/_R2_/')
    SAMPLE=$(basename $R1 | sed 's/_R1_.*//')
    
    echo "Processing sample: $SAMPLE"
    
    fastp \
        -i $R1 \
        -I $R2 \
        -o 04_trimmed/${SAMPLE}_R1.trimmed.fastq.gz \
        -O 04_trimmed/${SAMPLE}_R2.trimmed.fastq.gz \
        --html 04_trimmed/${SAMPLE}.fastp.html \
        --json 04_trimmed/${SAMPLE}.fastp.json \
        --thread $THREADS \
        --detect_adapter_for_pe
done

# --- 3. POST-TRIMMING QC (FASTQC) ---
echo "STEP 3: Running FastQC on trimmed data..."
fastqc -t $THREADS -o 05_qc_trimmed 04_trimmed/*.trimmed.fastq.gz

# --- 4. GENOME INDEXING (STAR) ---
echo "STEP 4: Generating STAR genome index (if it does not exist)..."
if [ ! -f "06_star_index/SA" ]; then
    STAR --runMode genomeGenerate \
         --runThreadN $THREADS \
         --genomeDir 06_star_index \
         --genomeFastaFiles $GENOME_FASTA \
         --sjdbGTFfile $GENOME_GTF \
         --sjdbOverhang 100
else
    echo "STAR index already exists. Skipping."
fi

# --- 5. ALIGNMENT (STAR) ---
echo "STEP 5: Aligning reads with STAR..."
mkdir -p 07_star_aligned/STAR_temp

for R1 in 04_trimmed/*_R1.trimmed.fastq.gz; do
    R2=$(echo $R1 | sed 's/_R1/_R2/')
    SAMPLE=$(basename $R1 | sed 's/_R1.trimmed.fastq.gz//')
    
    echo "Aligning sample: $SAMPLE"
    
    STAR --runMode alignReads \
         --runThreadN $THREADS \
         --genomeDir 06_star_index \
         --readFilesIn $R1 $R2 \
         --readFilesCommand zcat \
         --outFileNamePrefix 07_star_aligned/${SAMPLE}_ \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes Standard \
         --outTmpDir 07_star_aligned/STAR_temp/${SAMPLE}_
done

# --- 6. QUANTIFICATION (featureCounts) ---
echo "STEP 6: Generating count matrix with featureCounts..."
featureCounts \
    -T $THREADS \
    -p \
    -s $STRANDNESS \
    -a $GENOME_GTF \
    -o 08_featurecounts/counts.txt \
    07_star_aligned/*.bam

# --- 7. AGGREGATED QC (MultiQC) ---
echo "STEP 7: Generating final MultiQC report..."
multiqc . -o 09_multiqc -f

echo "--- PIPELINE COMPLETED! ---"
echo "Final count matrix: 08_featurecounts/counts.txt"
echo "Final QC report: 09_multiqc/multiqc_report.html"