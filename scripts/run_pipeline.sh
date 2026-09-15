#!/bin/bash
set -e # El script fallará si un comando falla

# --- 0. USER CONFIGURATION (ADJUST THESE!) ---

# Configura el número de hilos (cores)
THREADS=64

# --- Rutas a tus archivos de genoma y anotación ---
# (¡Estos ya están configurados con los nombres que proporcionaste!)
GENOME_FASTA="02_genome/GRCm38.p4.genome.fa"
GENOME_GTF="02_genome/gencode.vM10.annotation.gtf"

# --- Library Strandness (CRÍTICO para featureCounts) ---
# 0 = Sin hebra (Unstranded)
# 1 = Con hebra (Stranded)
# 2 = Con hebra reversa (Reversely Stranded) (Más común para Illumina)
# ¡Déjalo en 2! Si obtienes 0 conteos, cambia esto a 1 o 0.
STRANDNESS=2

# --- END OF CONFIGURATION ---

echo "--- INICIANDO PIPELINE DE RNA-SEQ ---"
echo "Usando $THREADS hilos."

# --- 1. QC INICIAL (FASTQC) ---
echo "PASO 1: Ejecutando FastQC en datos crudos..."
fastqc -t $THREADS -o 03_qc_raw 01_data_raw/*.fastq.gz

# --- 2. TRIMMING (FASTP) ---
echo "PASO 2: Ejecutando fastp para trimming..."
for R1 in 01_data_raw/*_R1_*.fastq.gz; do
    R2=$(echo $R1 | sed 's/_R1_/_R2_/')
    SAMPLE=$(basename $R1 | sed 's/_R1_.*//')
    
    echo "Procesando muestra: $SAMPLE"
    
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

# --- 3. QC POST-TRIMMING (FASTQC) ---
echo "PASO 3: Ejecutando FastQC en datos trimmeados..."
fastqc -t $THREADS -o 05_qc_trimmed 04_trimmed/*.trimmed.fastq.gz

# --- 4. ÍNDICE DEL GENOMA (STAR) ---
echo "PASO 4: Creando índice del genoma STAR (si no existe)..."
if [ ! -f "06_star_index/SA" ]; then
    STAR --runMode genomeGenerate \
         --runThreadN $THREADS \
         --genomeDir 06_star_index \
         --genomeFastaFiles $GENOME_FASTA \
         --sjdbGTFfile $GENOME_GTF \
         --sjdbOverhang 100 # Ideal: ReadLength-1. 100 es un valor seguro.
else
    echo "Índice de STAR ya existe. Saltando."
fi

# --- 5. ALINEAMIENTO (STAR) ---
echo "STEP 5: Aligning reads with STAR..."

# --- AÑADE ESTA LÍNEA ---
# Crea un directorio temporal base para la clasificación de STAR
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
         --outTmpDir 07_star_aligned/STAR_temp/${SAMPLE}_ # <-- AÑADE ESTA LÍNEA (incluye la '\' en la línea anterior)
done
# --- 6. CONTEO (featureCounts) ---
echo "PASO 6: Generando tabla de conteos con featureCounts..."
featureCounts \
    -T $THREADS \
    -p # -p especifica lecturas paired-end
    -s $STRANDNESS \
    -a $GENOME_GTF \
    -o 08_featurecounts/counts.txt \
# --- 7. QC AGREGADO (MultiQC) ---
echo "PASO 7: Generando reporte final MultiQC..."
multiqc . -o 09_multiqc -f

echo "--- ¡PIPELINE COMPLETADO! ---"
echo "Tabla de conteos final en: 08_featurecounts/counts.txt"
echo "Reporte de QC final en: 09_multiqc/multiqc_report.html"
