module load apps/cellranger/6.1.2

cellranger count \
    --id=rca1-1 \
    --sample=rca1-1 \
    --fastqs=/home/pingxr/Atherosis_0723/data_human/FASTQ/rca1-1 \
    --transcriptome=/DATA/public/refdata-gex-GRCh38-2020-A \
    --localmem=60 \
    --localcores=10
    
cellranger count \
--id=rca1-2 \
--sample=rca1-2 \
--fastqs=/home/pingxr/Atherosis_0723/data_human/FASTQ/rca1-2 \
--transcriptome=/DATA/public/refdata-gex-GRCh38-2020-A \
--localmem=60 \
--localcores=10

cat rca2-1_S1_L004_R1_001.fastq rca2-1_S1_L005_R1_001.fastq rca2-1_S1_L006_R1_001.fastq >rca2-1_S1_L006_R1_001.fastq
cat rca2-1_S1_L004_R2_001.fastq rca2-1_S1_L005_R2_001.fastq rca2-1_S1_L006_R2_001.fastq >rca2-1_S1_L006_R2_001.fastq
cellranger count \
    --id=rca2-1 \
    --sample=rca2-1 \
    --fastqs=/home/pingxr/Atherosis_0723/data_human/FASTQ/rca2-1 \
    --transcriptome=/DATA/public/refdata-gex-GRCh38-2020-A \
    --localmem=60 \
    --localcores=10

cat rca2-2_S2_L004_R1_001.fastq rca2-2_S2_L005_R1_001.fastq rca2-2_S2_L006_R1_001.fastq >rca2-2_S2_L006_R1_001.fastq
cat rca2-2_S2_L004_R2_001.fastq rca2-2_S2_L005_R2_001.fastq rca2-2_S2_L006_R2_001.fastq >rca2-2_S2_L006_R2_001.fastq
cellranger count \
    --id=rca2-2 \
    --sample=rca2-2 \
    --fastqs=/home/pingxr/Atherosis_0723/data_human/FASTQ/rca2-2 \
    --transcriptome=/DATA/public/refdata-gex-GRCh38-2020-A \
    --localmem=60 \
    --localcores=10

cat rca3-1_S6_L004_R1_001.fastq rca3-1_S6_L005_R1_001.fastq rca3-1_S6_L006_R1_001.fastq >rca3-1_S6_L006_R1_001.fastq
cat rca3-1_S6_L004_R2_001.fastq rca3-1_S6_L005_R2_001.fastq rca3-1_S6_L006_R2_001.fastq >rca3-1_S6_L006_R2_001.fastq
cellranger count \
    --id=rca3-1 \
    --sample=rca3-1 \
    --fastqs=/home/pingxr/Atherosis_0723/data_human/FASTQ/rca3-1 \
    --transcriptome=/DATA/public/refdata-gex-GRCh38-2020-A \
    --localmem=60 \
    --localcores=10

cat rca3-2_S7_L004_R1_001.fastq rca3-2_S7_L005_R1_001.fastq rca3-2_S7_L006_R1_001.fastq >rca3-2_S7_L006_R1_001.fastq
cat rca3-2_S7_L004_R2_001.fastq rca3-2_S7_L005_R2_001.fastq rca3-2_S7_L006_R2_001.fastq >rca3-2_S7_L006_R2_001.fastq
cellranger count \
    --id=rca3-2 \
    --sample=rca3-2 \
    --fastqs=/home/pingxr/Atherosis_0723/data_human/FASTQ/rca3-2 \
    --transcriptome=/DATA/public/refdata-gex-GRCh38-2020-A \
    --localmem=60 \
    --localcores=10

cat rca3-3_S8_L004_R1_001.fastq rca3-3_S8_L005_R1_001.fastq rca3-3_S8_L006_R1_001.fastq >rca3-3_S8_L006_R1_001.fastq
cat rca3-3_S8_L004_R2_001.fastq rca3-3_S8_L005_R2_001.fastq rca3-3_S8_L006_R2_001.fastq >rca3-3_S8_L006_R2_001.fastq
cellranger count \
    --id=rca3-3 \
    --sample=rca3-3 \
    --fastqs=/home/pingxr/Atherosis_0723/data_human/FASTQ/rca3-3 \
    --transcriptome=/DATA/public/refdata-gex-GRCh38-2020-A \
    --localmem=60 \
    --localcores=10

cellranger count \
    --id=rca4-1 \
    --sample=rca4-1 \
    --fastqs=/home/pingxr/Atherosis_0723/data_human/FASTQ/rca4-1 \
    --transcriptome=/DATA/public/refdata-gex-GRCh38-2020-A \
    --localmem=60 \
    --localcores=10
