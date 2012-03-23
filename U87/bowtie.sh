echo "bowtie -p 4 -n 0 -m 1 -k 1 --best ../hg19/bt1 U87.fq.noAdapters > U87.fq.noAdapters.n0m1k1b.bt"
bowtie -p 4 -n 0 -m 1 -k 1 --best ../hg19/bt1 U87.fq.noAdapters > U87.fq.noAdapters.n0m1k1b.bt
echo "bowtie -p 4 -n 0 -m 20 -k 20 --best ../hg19/bt1 U87.fq.noAdapters > U87.fq.noAdapters.n0m20k20b.bt"
bowtie -p 4 -n 0 -m 20 -k 20 --best ../hg19/bt1 U87.fq.noAdapters > U87.fq.noAdapters.n0m20k20b.bt
