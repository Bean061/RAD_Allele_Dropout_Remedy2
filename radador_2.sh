#!/bin/bash

set -e  # Exit if any command fails

##############################
# User-defined paths
##############################
# make sure you have 50.loci file from ipyrad.
# make sure you installed geneminer.py

input="/work/users/w/z/wzhou10/venus_flytrap/50_outfiles/50.loci"          # Your loci input from ipyrad
output="./output_dir"              # Output directory
split_flag="2"                     # paired-end is 2; single-end is 1
reads1="../../Aldrovanda_vesiculosa/trimmed-rename_12.fastq"    # QC-ed read file_R1
reads2="../../Aldrovanda_vesiculosa/trimmed-rename_22.fastq"    # QC-ed read file_R2
geneminer="/work/users/w/z/wzhou10/tools/GeneMiner/geneminer.py"    # Full path of your geneminer.py
sample_name="Al_ves"               # Change your outgroup name with less than 9 characters
threads=20                         # Threads of your CPUs

##############################
# Step 1: Split loci to fasta
##############################
# run python script from the folder ../s1_Loci_2_fasta.py
# no external module is needed for this step.

mkdir -p "$output"
cd "$output"

if [ ! -d "R1R2_output" ] ||  [ -z "$(ls -A "R1R2_output")" ] ; then
    echo "Step1: Running s1_Loci_2_fasta.py..."
    python ../s1_Loci_2_fasta.py -i "$input"
    touch Step1_loci_2_fasta_complete.flag
else
    echo "Step 1 already completed, skipping."
fi

##############################
# Step 2: Optional split by NNN+
##############################
# run python script from the folder ../s2_split_R1_R2.py
# install biopython 1.79

if [ ! -d "R1_fasta" ] || [ -z "$(ls -A "R1_fasta")" ] || [ ! -d "R2_fasta" ] || [ -z "$(ls -A "R2_fasta")" ]; then
    echo "Step2: Running s2_split_R1_R2.py..."
    python ../s2_split_R1_R2.py -i R1R2_output -o R1_fasta -p R2_fasta
    touch Step2_split_complete.flag
else
    echo "Step 2 outputs exist, skipping."
fi

##############################
# Step 3: Run GeneMiner
##############################
# install GeneMiner https://github.com/happywithxpl/GeneMiner
# confirm you changed the geneminer path in the input command line.

if [ ! -d geneMiner_output_R1 ] || [ -z "$(ls -A "geneMiner_output_R1")" ]; then
    echo "Step3: Running GeneMiner for R1..."
    python "$geneminer" -1 "$reads1" -2 "$reads2" -rtfa ./R1_fasta -o geneMiner_output_R1 -t $threads -k1 21 -k2 21
    touch Step3_geneminer_R1_complete.flag
else
    echo "GeneMiner R1 output exists, skipping."
fi

if [ ! -d geneMiner_output_R2 ] || [ -z "$(ls -A "geneMiner_output_R2")" ]; then
    echo "Step3: Running GeneMiner for R2..."
    python "$geneminer" -1 "$reads1" -2 "$reads2" -rtfa ./R2_fasta -o geneMiner_output_R2 -t $threads -k1 21 -k2 21
    touch Step3_geneminer_R2_complete.flag
else
    echo "GeneMiner R2 output exists, skipping."
fi

##############################
# Step 4: Combine and format loci with more than 5 read covered
##############################
#

if [ ! -d split_genes_R1 ] || [ -z "$(ls -A "split_genes_R1")" ] || [ ! -d split_genes_R2 ] || [ -z "$(ls -A "split_genes_R2")" ]; then
    mkdir -p split_genes_R1 split_genes_R2
    echo "Processing R1 loci..."
    touch passed_genes_R1.txt
    for file in geneMiner_output_R1/GM_results/*.fasta; do
        name=$(basename "${file}" .fasta)
        filtered_file="geneMiner_output_R1/filtered_out/${name}.fasta"

        # Check if filtered file exists and has more than 5 reads
        if [ -f "$filtered_file" ]; then
            read_count=$(grep -c '^>' "$filtered_file")
            if [ "$read_count" -ge 5 ]; then
                # Add locus to alignment
                awk '/^>/ {print ">'"$sample_name"'"} !/^>/ {print}' "$file" > split_genes_R1/${name}.fasta
                cat "./R1_fasta/${name}.fasta" >> split_genes_R1/${name}.fasta
		echo "$name" >> passed_genes_R1.txt
            fi
        fi
    done

    touch Step4_combined_outgroup_R1_complete.flag

    echo "Processing R2 loci..."
    touch passed_genes_R2.txt
    for file in geneMiner_output_R2/GM_results/*.fasta; do
        name=$(basename "${file}" .fasta)
	filtered_file="geneMiner_output_R2/filtered_out/${name}.fasta"

        # Check if filtered file exists and has more than 5 reads
        if [ -f "$filtered_file" ]; then
            read_count=$(grep -c '^>' "$filtered_file")
            if [ "$read_count" -ge 5 ]; then
                # Add locus to alignment
                awk '/^>/ {print ">'"$sample_name"'"} !/^>/ {print}' "$file" > split_genes_R2/${name}.fasta
		cat "./R2_fasta/${name}.fasta" >> split_genes_R2/${name}.fasta
	        echo "$name" >> passed_genes_R2.txt
	    fi
        fi
    done

    touch Step4_combined_outgroup_R2_complete.flag
else
    echo "combined outgroup sequences output exists, skipping."
fi

##############################
# Step 5: Align with MAFFT
##############################
# install MAFFT to get the MSA
#

if [ ! -d aligned_R1 ] || [ -z "$(ls -A "aligned_R1")" ] || [ ! -d aligned_R2 ] || [ -z "$(ls -A "aligned_R2")" ]; then
    mkdir -p aligned_R1 aligned_R2
    echo "Aligning R1 loci with MAFFT..."
    for file in split_genes_R1/*.fasta; do
        name=$(basename "${file}" .fasta)
        if [ ! -f aligned_R1/${name}.aln ]; then
            mafft --adjustdirection --localpair "$file" > aligned_R1/${name}.aln
            sed 's/_R_//g' "aligned_R1/${name}.aln" > "aligned_R1/${name}.tmp"
            mv "aligned_R1/${name}.tmp" "aligned_R1/${name}.aln"
        fi
    done
    
    touch Step5_MAFFT_R1_complete.flag
    
    echo "Aligning R2 loci with MAFFT..."
    for file in split_genes_R2/*.fasta; do
        name=$(basename "${file}" .fasta)
        if [ ! -f aligned_R2/${name}.aln ]; then
            mafft --adjustdirection --localpair --thread $threads "$file" > aligned_R2/${name}.aln
            sed 's/_R_//g' "aligned_R2/${name}.aln" > "aligned_R2/${name}.tmp"
            mv "aligned_R2/${name}.tmp" "aligned_R2/${name}.aln"
        fi
    done
    
    touch Step5_MAFFT_R2_complete.flag
else
    echo "MAFFT output exists, skipping."
fi

#########################################
# Step 6: Alignment Trimming with trimAl
#########################################

# install trimAl to trim the alignment.
# Trim R1 alignments

if [ ! -d trimmed_aligned_R1 ] || [ -z "$(ls -A "trimmed_aligned_R1")" ]; then
    mkdir -p trimmed_aligned_R1
    for file in aligned_R1/*.aln; do
        name=$(basename "$file" .aln)
        trimal -in "$file" -out "trimmed_aligned_R1/${name}.fasta" -automated1
    done
    touch Step6_trimmed_R1_complete.flag
else
    echo "trimmed_aligned_R1 exists, skipping trimming for R1."
fi

# Trim R2 alignments
if [ ! -d trimmed_aligned_R2 ] || [ -z "$(ls -A "trimmed_aligned_R2")" ]; then
    mkdir -p trimmed_aligned_R2
    for file in aligned_R2/*.aln; do
        name=$(basename "$file" .aln)
        trimal -in "$file" -out "trimmed_aligned_R2/${name}.fasta" -automated1
    done
    touch Step6_trimmed_R2_complete.flag
else
    echo "trimmed_aligned_R2 exists, skipping trimming for R2."
fi


echo "Pipeline completed successfully."
