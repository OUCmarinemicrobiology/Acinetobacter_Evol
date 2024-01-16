# Acinetobacter_Evol
# Interspecific host tropism and functional differentiation in the genus Acinetobacter

This repository contains a collection of code and scripts used in the paper **Interspecific host tropism and functional differentiation in the genus Acinetobacter** by Chun-Xu Xue, Kai Zhou,  et al..

The links below the sub-headings lead to the scripts needed for the corresponding steps. Most of the scripts were developed for running on the Huawei FusionServer Pro 5885H V5 server. You may download and adapt the scripts to suit your requirements.

## 1. Software used in this workflow

- [Perl](https://www.perl.org/)
- [Python3](https://www.python.org/)
- [Trimmomatic](https://github.com/timflutre/trimmomatic)
- [SPAdes](https://github.com/ablab/spades)
- [Unicycle](https://github.com/rrwick/Unicycler)
- [Prokka](https://github.com/tseemann/prokka)
- [Prodigal](https://github.com/hyattpd/Prodigal)
- [HMMER](http://hmmer.org/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [RAxML](https://evomics.org/learning/phylogenetics/raxml/)
- [Pandas](https://pandas.pydata.org/)
- [Numpy](https://numpy.org/)
- [CheckM](https://ecogenomics.github.io/CheckM/)
- [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk)
- [fastANI](https://github.com/ParBLiSS/FastANI)
- [Abricate](https://github.com/tseemann/abricate)
- [Resfinder](https://github.com/cadms/resfinder)
- [Snippy](https://github.com/tseemann/snippy)
- [Gubbins](https://github.com/nickjcroucher/gubbins)
- [Phytools](https://cran.r-project.org/web/packages/phytools/index.html)
- [Seqkit](https://bioinf.shenwei.me/seqkit/)
- [Hyphy](https://github.com/veg/hyphy-analyses/tree/master/FitMG94)

>Take the A-baumannii-104 isolate as an example.

## 2. Dataset
All assembled Illumina sequence data have been deposited in GenBank under the BioProject accession number [PRJNA663756](https://ncbi.nlm.nih.gov/bioproject/663756) and [PRJNA1063895](https://ncbi.nlm.nih.gov/bioproject/063895).

## 3. read trimming
### Trimmomatic
```bash
java -jar trimmomatic-0.36.jar PE -threads 5 A-baumannii-104_raw_1.fq.gz A-baumannii-104_raw_2.fq.gz A-baumannii-104_clean_1.fq.gz A-baumannii-104__unpaired_1.fq.gz A-baumannii-104_clean_2.fq.gz A-baumannii-104_unpaired_2.fq.gz
```

## 4. Assembly
### SPAdes
```bash
spades.py -1 A-baumannii-104_clean_1.fq.gz -2 A-baumannii-104_clean_2.fq.gz --isolate --cov-cutoff auto -o A-baumannii-104.fasta
```
### Unicycle
```bash
unicycler -1 A-baumannii-104_1.clean_1.fq.gz -2 A-baumannii-104_2.clean_1.fq.gz -l A-baumannii-104.nanopore.fq.gz -o A-baumannii-104.unicycle.fasta
```

## 5. Genome dereplication
### dRep
```bash
nohup dRep dereplicate fasta_dir.dRep_comp_0.75_nc_30.out -sa 0.95 -nc 0.30 -p 24 -comp 70 -con 10 -g fasta_dir/*.fasta &
# fasta_dir, the input directory containing a set of genomic assembly sequences.
# fasta_dir.dRep_comp_0.75_nc_30.out, output directory.
```

## 6. Taxonomy assignment
### GTDB
```bash
nohup gtdbtk classify_wf --genome_dir fasta_dir/ --out_dir fasta_dir.GTDB.out --extension fasta &
# fasta_dir, the input directory containing a set of genomic assembly sequences.
# fasta_dir.GTDB.out, output directory
```

## 7. Amino acid identity (ANI) calculation
### fastANI
```bash
fastANI --ql quer_genome.list --rl ref_genome.list -o FastANI.out -t 40
```

## 8. Genome annotation
### Prokka
```bash
prokka A-baumannii-104.fasta --prefix A-baumannii-104 --outdir A-baumannii-104.prokka.out/KP16932 --compliant
```

## 9. Identification of ARGs,
### Abricate
```bash
mkdir ARG_dir
for f in `ls fasta_dir`; do abricate -db resfinder --nopath --minid 50 --mincov 70 --quiet fasta_dir/${f} > ARG_dir/${f%%.fasta}.tab; done
abricate --nopath --summary ARG_dir/*tab > ARG.tab
# fasta_dir, the input directory containing a set of genomic assembly sequences.
```

## 10. Core SNP Phylogenetic analysis
### Snippy, Gubbins, RAxML
```bash
# call SNPs for multiple isolates from the same reference KP16932.
snippy-multi input.tab --ref KP16932.fa  --cpu 24  > runme.sh
# input.tab, a tab separated input file as follows
# input.tab = ID assembly.fasta
# Isolate	/path/to/contigs.fasta
less runme.sh   # check the script makes sense
sh ./runme.sh   # leave it running over lunch

# remove all the "weird" characters and replace them with N
snippy-clean_full_aln core.full.aln > clean.full.aln 

### Gubbins
# detect recombination region
run_gubbins.py -f 50 -p gubbins clean.full.aln

# remove recombination region
snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln
# -c only output columns containing exclusively ACGT

### RAxML
# build core SNP tree
raxmlHPC -f a -x 12345 -p 12345 -# 100 -m GTRGAMMAX -s clean.core.aln -n tree
```

## 11. Ancestral state reconstruction of gene number
### Phytools, R
```R
# phylogenetic_tree.nwk and gene_count.csv can be found in `Gene_ancestral_state` dictionary.
setwd("path/work_dictionary")
library(phytools)

tree <- read.tree("phylogenetic_tree.nwk")

mge <- read.csv("gene_count.csv",row.names=1) # input gene number of each isolates.
mge<-as.matrix(mge)[,1] # selected first column as example

# estimate ancestral states and compute variances & 95% confidence intervals for each node:
fit<-fastAnc(tree,mge,vars=TRUE,CI=TRUE)
fit

# projection of the reconstruction onto the edges of the tree
obj<-contMap(tree,mge,plot=FALSE)
plot(obj,legend=0.7*max(nodeHeights(tree)),
     fsize=c(0.1,0.9), lwd=1, outline = F, leg.txt="gene_count",ftype="off")
# fsize, set font size
# outline, logical value indicated whether or not to outline the plotted color bar with a 1 pt
line.
# leg.txt, set title of legend.
# ftype="off", don't show leave's name.

# OR set colors manually
obj<-setMap(obj,c("red", "#fffc00", "green", "purple", "blue", "#d7ff00", "black"))
plot(obj,legend=0.7*max(nodeHeights(tree)),
     fsize=c(0.1,0.9), lwd=1, outline = F, leg.txt="gene_count",ftype="off")
```

## 11. Habitat switching
### Phytools, R
```R
# phylogenetic_tree.nwk and genome_habitat.csv can be found in `Habitat_switching` dictionary.
setwd("path/work_dictionary")
library(phytools)

tree <- read.tree("phylogenetic_tree.nwk")

mge <- read.csv("genome_habitat.csv",row.names=1) # input gene number of each isolates.
mge<-as.matrix(mge)[,1] # selected first column as example
mtree<-make.simmap(tree,x,model="ER",,nsim=100)
mtrees
par(mfrow=c(10,10))
null<-sapply(mtrees,plot,colors=cols,lwd=1,ftype="off")
pd<-summary(mtrees,plot=FALSE)
pd
cols <- read.csv("cols.csv",row.names =1) #set color
cols <- as.matrix(cols)[,1]
plot(mtrees[[1]],cols,type="fan",fsize=0.8,ftype="off")
```

## 12. Selection analysis
### Hyphy
#### Whole genome selection analysis
```bash
hyphy FitMG94.bf --alignment gene_families_alignment_concat.fasta --tree phylogenetic_tree.nwk --lrt Yes --type local
# gene_families_alignment_concat.fasta file contain core gene sequence alignments
```

#### KatG selection analysis
```bash
hyphy absrel --alignment katG_gene_alignment.fasta --tree phylogenetic_tree.nwk --output katG.hyphy-absrel.out > katG.hyphy-absrel.nohup.out
#katG_gene_alignment.fasta file contain katG gene sequence aligemnts
```




