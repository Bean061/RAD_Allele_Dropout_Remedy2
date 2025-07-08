# RAD-seq Allele Dropout Remedy (RADADOR) version 2
A Pipeline Using Whole Genome Sequencing Data (~10x) To Obtain Shared RAD-seq Loci From Newly Sequenced Samples.
by Wenbin Zhou

It can help obtain shared loci from newly sequenced samples and help resolve the phylogeny of new samples using WGS data.

# Prerequisites:
Run [ipyrad](https://ipyrad.readthedocs.io/en/latest/) for RAD-seq data to get the .loci file

# Dependencies
python 3
Biopython 1.79 or later
GeneMiner
MAFFT
TtrimAl

# installation
[argparse](https://pypi.org/project/argparse/). Easy installation from [conda](https://anaconda.org/conda-forge/argparse)

[MAFFT](https://mafft.cbrc.jp/alignment/software/). Easy installation from [conda](https://anaconda.org/bioconda/mafft)

# Environment
Examples can be run on Mac or Linux.


## Details and Usage

  This script is used to recapture loci in RAD-seq outgroup from RNA-seq or Genome data. It requires .loci file from ipyrad, .phy file from ipyrad, .fasta file from RNA-seq (Trinity) or Genomic data. You also need to define one outgroup name and output partition file name. '-i', '-o', '-iph', '-itr', '-og' are required arguments. The default output file will be generated at current working dirctory.
  
  Make sure your working dirctory contain all the RADADOR python scripts. You can copy all files to RAD_Allele_Dropout_Remedy folder.
  
  ``` 
  cd RAD_Allele_Dropout_Remedy2/
  vim radador_2.sh

  ```
  Modify all parameters in User-defined paths
  ```
  input="/work/users/w/z/wzhou10/venus_flytrap/50_outfiles/50.loci"          # Your loci input from ipyrad
  output="./output_dir"              # Output directory
  split_flag="2"                     # paired-end is 2; single-end is 1
  reads1="../../Aldrovanda_vesiculosa/trimmed-rename_12.fastq"    # QC-ed read file_R1
  reads2="../../Aldrovanda_vesiculosa/trimmed-rename_22.fastq"    # QC-ed read file_R2
  geneminer="/work/users/w/z/wzhou10/tools/GeneMiner/geneminer.py"    # Full path of your geneminer.py
  sample_name="Al_ves"               # Change your outgroup name with less than 9 characters
  threads=20                         # Threads of your CPUs
  ```

  Then run command
  
  ```bash
  sh radador_2.sh
  ```

## Output
* R1R2_output
    fasta files for every locus
    
* R1_fasta and R2_fasta
    splited fasta files for every locus
  
* geneMiner_output_R1 and geneMiner_output_R2
    geneMiner output for newly sequenced individual from Whole Genome Sequencing data (~10x)
  
* split_genes_R1 and split_genes_R2
    splited files combined with newly sequenced sample. 
 
* aligned_R1 and aligned_R2
    Maftt align every locus

* trimmed_aligned_R1 and trimmed_aligned_R2
    trim all loci using trimAl

Finally you can concatenate trimmed_aligned_R1 and trimmed_aligned_R2 as the final matrix for phylogenetic analyses.


## Citation

* Zhou, W. in prep.
