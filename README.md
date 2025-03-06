# About
### This is the main script used for metagenomic data processing and statistical analysis in the following study:
Decade-Long Experimental Warming Accelerates Antibiotic Resistances in Grassland Soils. *Unpublished*.

For simplicity, we take sample S1 as an example:
# 1 Metagenomic data preprocessing
## 1.1 quality control
We used bbduk.sh in BBMap (https://sourceforge.net/projects/bbmap/) for QC.
```
bbduk.sh \
in1=S1_R1_001.fastq.gz in2=S1_R2_001.fastq.gz \
out=S1.clean.fastq \
ref=/path/to/adapters.fa,/path/to/phix_adapters.fa.gz \
ktrim=r k=23 mink=11 hdist=1 qtrim=r trimq=20 maq=20 -t=4 tbo
```

## 1.2 Assembly
We used MEGAHIT （https://github.com/voutcn/megahit）to assemble the contigs.
```
megahit --12 S1.clean.fastq \
-t 20 \
--k-list 31,51,71,91,111,131 \
--kmin-1pass \
--min-contig-len 500 \
--continue \
--out-prefix S1 \
-o S1
```

## 1.3 Binning
The BASALT pipeline (https://github.com/EMBL-PKU/BASALT) was employed for binning. To ensure compatibility with the input format of BASALT, the interleaved FASTQ file was initially converted into paired files by reformat.sh in BBMap.
```
reformat.sh in=S1.clean.fastq out1=S1.R1.fq out2=S1.R2.fq
BASALT -a S1.contigs.fa -s S1.R1.fq,S1.R2.fq -t 20 -m 230 --mode continue
```

# 2 ARG annotation
## 2.1 read-based ARG annotation
We used RGI bwt (https://github.com/arpcard/rgi/blob/master/docs/rgi_bwt.rst) to align metagenomic reads against CARD’s canonical reference and CARD’s Resistomes & Variants database.
```
rgi clean --local
rgi load --wildcard_annotation wildcard_database_v4.0.1.fasta \
--card_json /path/to/card.json \
--wildcard_index /path/to/index-for-model-sequences.txt \
--card_annotation card_database_v3.2.7.fasta \
--local

rgi bwt --read_one S1.R1.fq \
--read_two S1.R2.fq \
--output_file S1.card.bwt \
--local -n 20 --include_wildcard
```

## 2.2 contig-based ARG annotation
ARGs were also annotated from assembled contigs or MAGs.
ORFs were first predicted from the assembled contigs of each metagenome or each MAG by Prodigal:
```
prodigal -p meta -q -f gff -i S1.contigs.fa -o  S1.gff -a S1.faa -d S1.fna

#remove '*' at the end of some protein sequences
sed -i 's/\*$//' S1.faa
```
Then, ORFs were assigned as ARGs by applying the CARD RGI main software:
```
rgi clean --local
rgi load --card_json /path/to/card.json --local

rgi main --input_sequence S1.faa \
--output_file S1.card_orf \
--local --clean --low_quality -t protein --num_threads 8
```
The remaining unannotated genes were filtered and subsequently annotated with Resfams core database (https://www.dantaslab.org/resfams):
```
hmmsearch --cut_ga --cpu 4 --tblout S1.resfams /path/to/Resfams.hmm S1.faa
```
To calculate the abundance of args, coverage depth values were assigned to each ARG ORF based on the coverage depth of the corresponding contig, which was determined by CoverM (https://github.com/wwood/CoverM) using default settings:
```
coverm contig --interleaved S1.clean.fastq \
--reference S1.contigs.fa -t 10 \
--bam-file-cache-directory orf_bam \
-o S1.contigs.coverm.tsv
```

# 3 Host predictions
## 3.1 Contig taxonomic annotation by Kraken 2
To predict the taxonomy of ARG hosts, contigs containing ARG-encoding ORFs were classified using Kraken 2 (https://github.com/DerrickWood/kraken2):
```
kraken2 --db /path/to/k2_database/ S1.contigs.fa --threads 20 --use-names --output S1.kraken2
```

## 3.2 MAG taxonomy by GTDB-tk
```
gtdbtk classify_wf \
--genome_dir /path/to/MAGs/ \
--out_dir gtdbtk --cpus 40 -x fasta
```

# 4 ARG sequences clustering
Potential mobile ARGs were identified by detecting identical or nearly identical assembled sequences shared across different species. We first pooled the nucleotide sequences of ARG ORFs from contigs and MAGs together, and then clustered the sequences at 99% identity and 90% coverage using MMseqs2 (https://github.com/soedinglab/MMseqs2).
```
# fna = {single fasta file containing all ARG ORF nucleotide sequences}

mmseqs createdb ${fna} ${fna}.mmdb
mkdir nt_cluster_99_tmp nt_clusters_i99_c90
mmseqs cluster ${fna}.mmdb nt_clusters_i99_c90/nt_cluster_99 nt_cluster_99_tmp --threads 12 --min-seq-id 0.99 -c 0.90 --cov-mode 1
mmseqs createtsv ${fna}.mmdb ${fna}.mmdb nt_clusters_i99_c90/nt_cluster_99 nt_clusters_i99_c90/nt_cluster_99.tsv
mmseqs result2repseq ${fna}.mmdb nt_clusters_i99_c90/nt_cluster_99 nt_clusters_i99_c90/nt_cluster_99_rep
mmseqs result2flat ${fna}.mmdb ${fna}.mmdb nt_clusters_i99_c90/nt_cluster_99_rep nt_clusters_i99_c90/nt_cluster_99_rep.fasta
```
