#!/bin/bash

# Prepares all reference files needed for the inStrain analysis.
# This script should be run only ONCE before processing individual samples.
# It now uses the official parse_stb.py script for scaffold-to-bin mapping.

set -e

# --- 1. CONFIGURATION ---
readonly CONDA_ENV_NAME="instrain"
# Directory containing your high-quality MAGs in .fa or .fna format
readonly MAG_DIR="/home/wulab/project/new_warming/derep_bins_95/dereplicated_genomes"
# Main output directory where all results will be stored
readonly OUTPUT_DIR="/home/linwei/nws0920/19_instrain/selection_analysis_v2"
# Number of threads for multithreaded tools
readonly THREADS=60

# --- Activate Conda Environment ---
# Ensure conda is initialized in your shell environment.
# If this command fails, you might need to run 'conda init bash' once and restart your shell.
source $(conda info --base)/etc/profile.d/conda.sh
conda activate "$CONDA_ENV_NAME"

# --- 2. SETUP DIRECTORIES ---
echo "--- Setting up directories ---"
readonly GENE_DIR="$OUTPUT_DIR/genes"
mkdir -p "$OUTPUT_DIR/bowtie2_index" "$GENE_DIR"

# --- 3. CONSOLIDATE MAGS ---
readonly ALL_MAGS_FA="$OUTPUT_DIR/all_mags_reference.fasta"
echo "--- Consolidating MAGs into a single reference file ---"
if [ ! -f "$ALL_MAGS_FA" ]; then
    # Using 'find' is more robust for different file structures
    find "$MAG_DIR" -type f \( -name "*.fasta" -o -name "*.fna" \) -print0 | xargs -0 cat > "$ALL_MAGS_FA"
    echo "MAGs consolidated into $ALL_MAGS_FA"
else
    echo "$ALL_MAGS_FA already exists. Skipping."
fi

# --- 4. PREDICT GENES WITH PRODIGAL (PARALLEL) ---
readonly ALL_GENES_FNA="$OUTPUT_DIR/all_mags.genes.fna"
echo "--- Predicting genes for all MAGs using Prodigal ---"
if [ ! -f "$ALL_GENES_FNA" ]; then
    echo "Running Prodigal in parallel..."
    # Using find and parallel together for robustness
    find "$MAG_DIR" -type f \( -name "*.fasta" -o -name "*.fna" \) | \
    parallel -j "$THREADS" --bar "prodigal -i {} -o $GENE_DIR/{/.}.gff -a $GENE_DIR/{/.}.faa -d $GENE_DIR/{/.}.fna -m -p single > /dev/null"


    echo "Concatenating gene files..."
    find "$GENE_DIR" -type f -name "*.fna" -print0 | xargs -0 cat > "$ALL_GENES_FNA"
    echo "Gene prediction complete. All genes at $ALL_GENES_FNA"
else
    echo "Concatenated gene file $ALL_GENES_FNA already exists. Skipping prediction."
fi

# --- 5. BUILD BOWTIE2 INDEX ---
readonly BT2_INDEX="$OUTPUT_DIR/bowtie2_index/all_mags_bt2"
echo "--- Building Bowtie2 index ---"
if [ ! -f "${BT2_INDEX}.1.bt2" ]; then
    bowtie2-build "$ALL_MAGS_FA" "$BT2_INDEX" --large-index --threads "$THREADS"
    echo "Bowtie2 index built at $BT2_INDEX"
else
    echo "Bowtie2 index already exists. Skipping."
fi

# --- 6. CREATE SCAFFOLD-TO-BIN MAP  ---
readonly STB_FILE="$OUTPUT_DIR/scaffold_to_bin.tsv"
echo "--- Creating scaffold_to_bin.tsv map using parse_stb.py ---"
if [ ! -f "$STB_FILE" ]; then
    # This command uses the official inStrain helper script.
    # It takes all MAG files in the specified directory as input.
    # The --reverse flag makes the output format: bin_id -> scaffold_name
    # The -f flag specifies the input genome files. The wildcard '*' handles all MAGs.
    parse_stb.py --reverse -f "$MAG_DIR"/* -o "$STB_FILE"
    echo "scaffold_to_bin.tsv created at $STB_FILE"
else
    echo "scaffold_to_bin.tsv already exists. Skipping."
fi

echo -e "\n--- Preparation complete! ---"
echo "You can now run the per-sample processing script using GNU Parallel."
echo "Example command:"
echo "cat sample_list.txt | parallel -j 6 --bar ./process_single_sample.sh {}"
echo "--------------------------"

conda deactivate
