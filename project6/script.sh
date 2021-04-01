#!/usr/bin/bash

echo Start on $(date --iso-8601='seconds')

[ ! -d rawdata ] && mkdir rawdata
cd rawdata

echo
echo 1 Downloading data

for filename in SRR941816 SRR941817 SRR941818 SRR941819; do
    [ ! -f "$filename.fastq.gz" ] && wget "http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR941/$filename/$filename.fastq.gz"
done

if [ ! -f GCF_000146045.2_R64_genomic.fna ]
then
    [ ! -f GCF_000146045.2_R64_genomic.fna.gz ] && wget http://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
    gunzip GCF_000146045.2_R64_genomic.fna.gz
fi

if [ ! -f GCF_000146045.2_R64_genomic.gff ]
then
    [ ! -f GCF_000146045.2_R64_genomic.gff.gz ] && wget http://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz
    gunzip GCF_000146045.2_R64_genomic.gff.gz
fi

echo
echo 2 Build index and align
[ ! -f genome_index.1.ht2 ] && python3 ~/software/hisat2-2.1.0/hisat2-build GCF_000146045.2_R64_genomic.fna genome_index

for filename in SRR941816 SRR941817 SRR941818 SRR941819; do
    echo "$filename"
    [ ! -f "$filename.bam" ] && ~/software/hisat2-2.1.0/hisat2 -p 4 -x genome_index -U "$filename.fastq.gz" | samtools sort > "$filename.bam"
done

echo
echo 3 Feature count
[ ! -f GCF_000146045.2_R64_genomic.gtf ] && gffread GCF_000146045.2_R64_genomic.gff -T -o GCF_000146045.2_R64_genomic.gtf
[ ! -f features.tsv ] && featureCounts -g gene_id -a GCF_000146045.2_R64_genomic.gtf -o features.tsv *.bam
[ ! -f simple_counts.txt ] && cat features.tsv | cut -f 1,7-10 > simple_counts.txt

echo
echo 4 Find differentially expressed genes with Deseq2
[ ! -f result.txt ] && cat simple_counts.txt | R -f ../deseq2.r
[ ! -f output.pdf ] && cat norm-matrix-deseq2.txt | R -f ../draw-heatmap.r
[ ! -f genes.txt ] && head result.txt -n 50 | cut -f 1 > genes.txt
[ ! -f genes_clean.txt ] && sed 's/gene-//' genes.txt > genes_clean.txt

echo
echo End on $(date --iso-8601='seconds')
