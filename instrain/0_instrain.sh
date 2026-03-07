## 0. install the environment
conda create -n instrain
conda activate instrain
mamba install -c conda-forge -c bioconda -c defaults instrain
mamba install bioconda::samtools
mamba install bioconda::bowtie2
mamba install bioconda::prodigal
pip install drep --upgrade

#Script 1: 1_prepare_references.sh (One-time setup)
#This script handles all the one-off prep work: merging MAGs, gene prediction, building BWA2/Bowtie2 indices, and creating the scaffold_to_bin file.

#Script 2: 2_process_single_sample.sh (Per-sample processing)
#This script processes a single sample, going from alignment all the way to inStrain profiling. It's designed to be called by GNU Parallel.

## 1. Run the prep script (just once):
#Using dereplicated MAGs as reference genomes (95% ANI)
cd /home/linwei/nws0920/19_instrain
chmod +x 1_prepare_references.sh
./1_prepare_references.sh

## 2. Prepare the sample list and run in parallel:
conda activate instrain
chmod +x 2_process_single_sample.sh


cat sample_list.txt | parallel -j 15 --bar ./2_process_single_sample.sh {}
