echo "bowtie -p 4 -n 0 -m 1 -k 1 --best ../hg19/bt1 ALL.fastq.trimmed > ALL.fastq.trimmed.n0m1k1b.bt"
bowtie -p 4 -n 0 -m 1 -k 1 --best ../hg19/bt1 ALL.fastq.trimmed > ALL.fastq.trimmed.n0m1k1b.bt 2> ALL.fastq.trimmed.n0m1k1b.info
echo "bowtie -p 4 -n 0 -m 20 -k 20 --best ../hg19/bt1 ALL.fastq.trimmed > ALL.fastq.trimmed.n0m20k20b.bt"
bowtie -p 4 -n 0 -m 20 -k 20 --best ../hg19/bt1 ALL.fastq.trimmed > ALL.fastq.trimmed.n0m20k20b.bt 2> ALL.fastq.trimmed.n0m20k20b.info
