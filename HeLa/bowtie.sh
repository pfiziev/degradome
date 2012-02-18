echo "bowtie -p 4 -n 0 -m 1 -k 1 --best /home/pf/UCLA/databases/bowtie/hg19 ALL.fastq.trimmed > ALL.fastq.trimmed.n0m1k1b.bt"
bowtie -p 4 -n 0 -m 1 -k 1 --best /home/pf/UCLA/databases/bowtie/hg19 ALL.fastq.trimmed > ALL.fastq.trimmed.n0m1k1b.bt 2> ALL.fastq.trimmed.n0m1k1b.info
