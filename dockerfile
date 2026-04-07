# --------------------------
# Use stable LTS base
# --------------------------
FROM ubuntu:22.04

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV THREADS=4
ENV PATH="/usr/bin:${PATH}"

# --------------------------
# Install dependencies with root access
# --------------------------
RUN apt-get update && apt-get install -y \
    sudo wget curl unzip gzip \
    bwa hisat2 bowtie2 \
    samtools bcftools fastqc fastp multiqc \
    sra-toolkit \
    default-jre \
    python3 python3-pip \
    ca-certificates \
    pigz \
    openssl \
    libssl-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# --------------------------
# Update certificates + OpenSSL config
# --------------------------
RUN update-ca-certificates

# --------------------------
# Create directories (root owns everything)
# --------------------------
RUN mkdir -p /data /app/scripts

# Optional: create non-root user for later use
RUN useradd -m bioinfo && chown -R bioinfo:bioinfo /data /app

# --------------------------
# Set working directory
# --------------------------
WORKDIR /data

# --------------------------
# Download GRCh38 reference securely
# --------------------------
RUN wget --https-only -q https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz && \
    gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz && \
    mv Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38.fa

# --------------------------
# Indexing for BWA, HISAT2, Bowtie2
# --------------------------
RUN bwa index GRCh38.fa
RUN hisat2-build GRCh38.fa GRCh38_hisat2_index

# --------------------------
# Create pipeline script
# --------------------------
RUN printf '#!/bin/bash\n\
set -e\n\
THREADS=${THREADS:-4}\n\
SRRS=("SRR22044200" "SRR22044199")\n\
for SRR in "${SRRS[@]}"; do\n\
    echo "Processing $SRR with $THREADS threads..."\n\
    fasterq-dump $SRR -O /data --split-files --threads $THREADS\n\
    pigz -p $THREADS /data/${SRR}_1.fastq\n\
    pigz -p $THREADS /data/${SRR}_2.fastq\n\
    fastqc /data/${SRR}_1.fastq.gz /data/${SRR}_2.fastq.gz -o /data\n\
    fastp -i /data/${SRR}_1.fastq.gz -I /data/${SRR}_2.fastq.gz \\\n\
          -o /data/${SRR}_1.trimmed.fastq.gz \\\n\
          -O /data/${SRR}_2.trimmed.fastq.gz \\\n\
          -w $THREADS\n\
    hisat2 -p $THREADS -x GRCh38_hisat2_index \\\n\
        -1 /data/${SRR}_1.trimmed.fastq.gz \\\n\
        -2 /data/${SRR}_2.trimmed.fastq.gz \\\n\
        -S /data/${SRR}_hisat2.sam\n\
    bwa mem -t $THREADS GRCh38.fa \\\n\
        /data/${SRR}_1.trimmed.fastq.gz \\\n\
        /data/${SRR}_2.trimmed.fastq.gz > /data/${SRR}.sam\n\
    samtools view -@ $THREADS -Sb /data/${SRR}.sam | \\\n\
        samtools sort -@ $THREADS -o /data/${SRR}.bam\n\
    samtools index /data/${SRR}.bam\n\
    bcftools mpileup --threads $THREADS -f GRCh38.fa /data/${SRR}.bam | \\\n\
        bcftools call -mv -o /data/${SRR}.vcf\n\
    echo "$SRR completed"\n\
done\n\
multiqc /data -o /data/multiqc_report\n\
echo "Pipeline completed successfully!"\n' > /app/scripts/pipeline.sh

# --------------------------
# Give execution permission
# --------------------------
RUN chmod 750 /app/scripts/pipeline.sh

# --------------------------
# Default command
# --------------------------
CMD ["/app/scripts/pipeline.sh"]
