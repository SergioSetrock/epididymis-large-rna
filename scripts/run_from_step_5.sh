#!/bin/bash
set -e # El script fallará si un comando falla

# --- 0. USER CONFIGURATION (¡Necesario para los Pasos 5, 6, 7!) ---

# Configura el número de hilos (cores)
THREADS=64

# Rutas a tus archivos de genoma y anotación
GENOME_FASTA="02_genome/GRCm38.p4.genome.fa"
GENOME_GTF="02_genome/gencode.vM10.annotation.gtf"

# Library Strandness
STRANDNESS=2

# --- END OF CONFIGURATION ---

echo "--- RE-INICIANDO PIPELINE DESDE EL PASO 5 ---"
echo "Usando $THREADS hilos."

# --- 5. ALINEAMIENTO (STAR) ---
echo "PASO 5: Alineando lecturas con STAR..."

# Bucle para alinear. NO creamos el directorio temp, dejamos que STAR lo haga.
for R1 in 04_trimmed/*_R1.trimmed.fastq.gz; do
    R2=$(echo $R1 | sed 's/_R1/_R2/')
    SAMPLE=$(basename $R1 | sed 's/_R1.trimmed.fastq.gz//')
    
    echo "Alineando muestra: $SAMPLE"
    
    # Comprobamos si el BAM de esta muestra YA existe
    if [ -f "07_star_aligned/${SAMPLE}_Aligned.sortedByCoord.out.bam" ]; then
        echo "Archivo BAM para $SAMPLE ya existe. Saltando."
    else
        echo "Alineando $SAMPLE..."
        STAR --runMode alignReads \
             --runThreadN $THREADS \
             --genomeDir 06_star_index \
             --readFilesIn $R1 $R2 \
             --readFilesCommand zcat \
             --outFileNamePrefix 07_star_aligned/${SAMPLE}_ \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMattributes Standard \
             --outTmpDir 07_star_aligned/STAR_temp/
    fi
done

# --- 6. CONTEO (featureCounts) ---
echo "PASO 6: Generando tabla de conteos con featureCounts..."

# Este es el bloque de comando corregido.
# Asegúrate de que CADA LÍNEA, excepto la última, termine con una '\'
featureCounts \
    -T $THREADS \
    -p \
    -s $STRANDNESS \
    -a $GENOME_GTF \
    -o 08_featurecounts/counts.txt \
    07_star_aligned/*Aligned.sortedByCoord.out.bam

# --- 7. QC AGREGADO (MultiQC) ---
echo "PASO 7: Generando reporte final MultiQC..."
multiqc . -o 09_multiqc -f

echo "--- ¡PIPELINE (PASOS 5-7) COMPLETADO! ---"
echo "Tabla de conteos final en: 08_featurecounts/counts.txt"
echo "Reporte de QC final en: 09_multiqc/multiqc_report.html"
