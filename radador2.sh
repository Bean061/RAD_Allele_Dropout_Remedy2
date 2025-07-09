#!/bin/bash

set -euo pipefail

# Required variables (you must define them before running this script):
# - reads: folder with *.fastq files
# - sample_name_prefix: prefix for sample names
# - output: main output directory path
# - geneminer: path to GeneMiner script
# - threads: number of threads to use
# - locifile: loci file from ipyrad

reads="/work/users/w/z/wzhou10/venus_flytrap/genome/Aldrovanda_vesiculosa/"
sample_name_prefix="outgroup"
output="./results2/"
geneminer="/work/users/w/z/wzhou10/tools/GeneMiner/geneminer.py"
threads=4
locifile="/work/users/w/z/wzhou10/venus_flytrap/50_outfiles/50.loci"

# Start index for naming
i=1
mkdir -p "$output"
cd $output

##############################
# Step 1: Split loci to fasta
##############################

if [ ! -d "R1R2_output" ] || [ -z "$(ls -A "R1R2_output")" ]; then
    echo "Step1: Running s1_Loci_2_fasta.py..."
    python ../s1_Loci_2_fasta.py -i $locifile
    touch Step1_loci_2_fasta_complete.flag
else
    echo "Step 1 already completed, skipping."
fi

##############################
# Step 2: Optional split by NNN+
##############################

if [ ! -d "R1_fasta" ] || [ -z "$(ls -A "R1_fasta")" ] || [ ! -d "R2_fasta" ] || [ -z "$(ls -A "R2_fasta")" ]; then
    echo "Step2: Running s2_split_R1_R2.py..."
    python ../s2_split_R1_R2.py -i R1R2_output -o R1_fasta -p R2_fasta
    touch Step2_split_complete.flag
else
    echo "Step 2 outputs exist, skipping."
fi

# Loop over all read1 fastq files
for fqfile in "$reads"/*1.fastq; do
    fastq_name=$(basename "$fqfile" 1.fastq)
    reads1="${fastq_name}1.fastq"
    reads2="${fastq_name}2.fastq"

    sample_name="${sample_name_prefix}_${i}"
    output_1="${sample_name}"

    echo "🔄 Processing $sample_name"
    mkdir -p "$output_1"
    #pushd "$output_1" > /dev/null
    cd $output_1

    ##############################
    # Step 3: Run GeneMiner
    ##############################

    if [ ! -d geneMiner_output_R1 ] || [ -z "$(ls -A "geneMiner_output_R1")" ]; then
        echo "Step3: Running GeneMiner for R1..."
        python "$geneminer" -1 "${reads}/${reads1}" -2 "${reads}/${reads2}" -rtfa ../R1_fasta -o geneMiner_output_R1 -t "$threads" -k1 21 -k2 21
        touch Step3_geneminer_R1_complete.flag
    else
        echo "GeneMiner R1 output exists, skipping."
    fi

    if [ ! -d geneMiner_output_R2 ] || [ -z "$(ls -A "geneMiner_output_R2")" ]; then
        echo "Step3: Running GeneMiner for R2..."
        python "$geneminer" -1 "${reads}/$reads1" -2 "${reads}/$reads2" -rtfa ../R2_fasta -o geneMiner_output_R2 -t "$threads" -k1 21 -k2 21
        touch Step3_geneminer_R2_complete.flag
    else
        echo "GeneMiner R2 output exists, skipping."
    fi

    #popd > /dev/null
    cd ..
    i=$((i + 1))
    pwd
done

echo "✅ All samples processed."

##############################
# Step 4: Combine and format loci with more than 5 read covered
##############################
#

if [ ! -d ./split_genes_R1 ] || [ -z "$(ls -A "./split_genes_R1")" ] || [ ! -d ./split_genes_R2 ] || [ -z "$(ls -A "./split_genes_R2")" ]; then
    mkdir -p split_genes_R1 split_genes_R2
    echo "Processing R1 loci..."
    touch passed_genes_R1.txt
    for sp_file in ./"${sample_name_prefix}"*/; do
        sp_name=$(basename "${sp_file}")
        for gene_file in "${sp_name}"/geneMiner_output_R1/GM_results/*.fasta; do
            gene_name=$(basename "${gene_file}" .fasta)
            filtered_file="${sp_name}/geneMiner_output_R1/filtered_out/${gene_name}.fasta"

            # Check if filtered file exists and has more than 5 reads
            if [ -f "$filtered_file" ]; then
                read_count=$(grep -c '^>' "$filtered_file")
                if [ "$read_count" -ge 5 ]; then
                    # Add locus to alignment
                    awk '/^>/ {print ">'"$sp_name"'"} !/^>/ {print}' "$gene_file" >> ./split_genes_R1/${gene_name}.fasta
                    #cat "./R1_fasta/${name}.fasta" >> split_genes_R1/${name}.fasta
                    echo "$gene_name" >> "${sp_name}"_passed_genes_R1.txt
                fi
            fi
        done
    done

    touch Step4_combined_outgroup_R1_complete.flag

    echo "Processing R2 loci..."
    touch passed_genes_R2.txt
    for sp_file in ./"${sample_name_prefix}"*/; do
        sp_name=$(basename "${sp_file}")
        for gene_file in "${sp_name}"/geneMiner_output_R2/GM_results/*.fasta; do
            gene_name=$(basename "${gene_file}" .fasta)
            filtered_file="${sp_name}/geneMiner_output_R2/filtered_out/${gene_name}.fasta"

            # Check if filtered file exists and has more than 5 reads
            if [ -f "$filtered_file" ]; then
                read_count=$(grep -c '^>' "$filtered_file")
                if [ "$read_count" -ge 5 ]; then
                    # Add locus to alignment
                    awk '/^>/ {print ">'"$sp_name"'"} !/^>/ {print}' "$gene_file" >> ./split_genes_R2/${gene_name}.fasta
                    #cat "./R1_fasta/${name}.fasta" >> split_genes_R1/${name}.fasta
                    echo "$gene_name" >> "${sp_name}"_passed_genes_R2.txt
                fi
            fi
        done
    done


    for file in ./R1_fasta/*.fasta; do
        gene_name=$(basename "${file}" .fasta)
        cat $file >> ./split_genes_R1/${gene_name}.fasta
    done
    echo "R1 loci combined."

    for file in ./R2_fasta/*.fasta; do
        gene_name=$(basename "${file}" .fasta)
        cat $file >> ./split_genes_R2/${gene_name}.fasta
    done

    echo "R2 loci combined."

else
    echo "Combined R1 and R2 files exist, skipping."

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

