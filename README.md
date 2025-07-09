# RAD-seq Allele Dropout Remedy (RADADOR) version 2
A Pipeline Using Whole Genome Sequencing Data (~10x) To Obtain Shared RAD-seq Loci From Newly Sequenced Samples.
by Wenbin Zhou

It can help obtain shared loci from newly sequenced samples and help resolve the phylogeny of new samples using WGS data.

# Prerequisites:
Run [ipyrad](https://ipyrad.readthedocs.io/en/latest/) for RAD-seq data to get the .loci file

# Dependencies
python 3
Biopython 1.79 or later. via conda
[MAFFT](https://mafft.cbrc.jp/alignment/software/) via [conda](https://anaconda.org/bioconda/mafft)
[TtrimAl](https://trimal.readthedocs.io/en/latest/) via [conda](https://anaconda.org/bioconda/trimal).
[argparse](https://pypi.org/project/argparse/). Easy installation from [conda](https://anaconda.org/conda-forge/argparse)

```
conda create --name radador python=3.7
source activate radador
conda install biopython=1.79
conda install bioconda::mafft
conda install bioconda::trimal
conda install conda-forge::argparse

```

[GeneMiner](https://github.com/happywithxpl/GeneMiner).
```
wget -c https://github.com/happywithxpl/GeneMiner/releases/download/v1.0.0/GeneMiner_v1.0.0_linux.tar.gz
tar GeneMiner_v1.0.1_linux.tar.gz
cd GeneMiner_v1.0.1_linux
chmod 755 geneminer.py

```

# Environment
Examples can be run on Mac or Linux.


## Details and Usage

  This script is used to recapture loci in RAD-seq outgroup from RNA-seq or Genome data. It requires .loci file from ipyrad, .phy file from ipyrad, .fasta file from RNA-seq (Trinity) or Genomic data. You also need to define one outgroup name and output partition file name. '-i', '-o', '-iph', '-itr', '-og' are required arguments. The default output file will be generated at current working dirctory.
  
  Make sure your working dirctory contain all the RADADOR python scripts. You can copy all files to RAD_Allele_Dropout_Remedy folder.
  
  ``` 
  cd RAD_Allele_Dropout_Remedy2/
  vim radador_2.sh

  ```
  Modify all parameters below
  ```
# Required variables (you must define them before running this script):
# - reads: folder with *.fastq files, end up with 1.fastq and 2.fastq
# - sample_name_prefix: prefix for sample names
# - output: main output directory path
# - geneminer: path to GeneMiner script
# - threads: number of threads to use

reads="/FULL/PATH/OF/venus_flytrap/genome/Aldrovanda_vesiculosa/"
sample_name_prefix="outgroup"
output="./results2/"
geneminer="/FULL/PATH/OF/GeneMiner/geneminer.py"
threads=4
locifile="/FULL/PATH/OF/LOCIFILE/50.loci"
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
  
* all individual folder with geneMiner_output_R1 and geneMiner_output_R2
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
