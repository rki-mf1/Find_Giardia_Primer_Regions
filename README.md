# Find Giardia Primer Regions
The Snakemake pipeline aims to find potential giardia primer regions 
## General Workflow
- Pairwise alignment
    - one of the input sample is set as reference sequence
         - pairwise alignment of each input file against the reference using the all vs all comparison of nucleotide sequences tool nucmer from MUMer (hmmer v3.3)
        - output in delta format
    - the alignments are filtered using the MUMer delta-filter function
        - this script is a wrapper around nucmer that filters 1-to1 alignments
    - extract the filtered alignment information and summarize them in a *.txt
    - combine the alignment information of all samples into an alignment information matrix
    - after sorting the matrix, filter by length and whether the alignment positions can be found in all input samples

- creating multiple sequence alignments (MSA) of the filtered genomic regions

    - extracting the filtered sequence regions per sample in fasta format
    - computing MSAs for each sequence region using MAFFT (v7.471)
    - using a sliding window filter on the MSAs to determine if the input samples can be distinguished by DNA variations within the window 
        - can wie finde differentiating regions in the MSAs
     
- read mapping to determine the allele sequence heterogeneity (ASH) per position and sample

    - map the reads against the initial fasta files per sample using Bowtie2 (v2.4.1)
    - extract specific regions from the bam files using samtools (v1.3)
    - count reads and bases per position using pysamstats (v1.1.2)

- combine the data

    - combine the statistics in a final table
    - rank the MSAs by their statistics and report number of position of interest (POI) per alignment
 
  ## Set up
  The pipeline needs Snakemake and Biopython. We recommend to create a Conda environment:
  ```bash
  cd designated/path
  git clone https://github.com/rki-mf1/Find_Giardia_Primer_Regions.git
  cd Find_Giardia_Primer_Regions
  conda create -f pipeline_conda_env.yaml -n GiardiaEnv
  conda activate GiardiaEnv
  ```

  ## Run the Pipeline
  The Pipeline needs a config file (see config_template.yaml) and a samplesheet (see sample_template.txt - a list of the samples and sample names) as input.
   ```bash
  snakemake -s Snakefile --configfile /path/to/config_template.yaml --use-conda 
  ```
