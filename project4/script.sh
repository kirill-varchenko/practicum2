#!/usr/bin/bash

echo Start on $(date --iso-8601='seconds')

[ ! -d rawdata ] && mkdir rawdata
cd rawdata

echo
echo 1 Downloading data

[ ! -f YOKOZUNA-1.scaffolds.fa ] && wget http://kumamushi.org/data/YOKOZUNA-1.scaffolds.fa
[ ! -f scaffold015.fa ] && wget http://public.dobzhanskycenter.ru/mrayko/BIMM185/scaffold015.fa
[ ! -f augustus.whole.gff ] && wget http://public.dobzhanskycenter.ru/mrayko/BIMM185/augustus.whole.gff
[ ! -f getAnnoFasta.pl ] && wget http://augustus.gobics.de/binaries/scripts/getAnnoFasta.pl

[ ! -f augustus.whole.aa ] && ./getAnnoFasta.pl augustus.whole.gff

[ ! -f peptides.fa ] && wget http://public.dobzhanskycenter.ru/mrayko/BIMM185/peptides.fa

stats_fasta.awk augustus.whole.aa

echo
echo 2 Protein localization

makeblastdb -in augustus.whole.aa -parse_seqids -blastdb_version 5 -title "R. varieornatus" -dbtype prot
blastp -db augustus.whole.aa -query peptides.fa -outfmt '6 qseqid sseqid' -out results.out
cat results.out | cut -f 2 | sort | uniq > results.txt
xargs samtools faidx augustus.whole.aa < results.txt > blasted.fa

cat "../WoLF PSORT Prediction.html" | cut --field=1,4 --delimiter=' ' | grep ':' | sed -e 's/://' -e 's/ /\t/' > psort.out
awk -F"\t" '$3 != 0 {print $2,$6;}' ../HMMER.txt > HHMER.out
../parse_TargetP.py ../TargetP-output.json > TargetP.out
blastp -db swissprot -query blasted.fa -negative_taxids 947166 -outfmt '6 qacc sacc evalue bitscore pident ssciname stitle qcovs' -out blast.tsv


