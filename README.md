# Genome Recombination Identification (Under development)
This program aims to discover the relation between a recombined genome sequence versus multiple reference sequences. Currently the program is under development and may containing unexpected errors.


## Dependencies
The program requires following dependencies. If you use anaconda/miniconda, you can run the following to initialize the environment.

```bash
# add channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create conda environment
conda create --name gene-recomb-env

# activate conda environment
conda activate gene-recomb-env

conda install -c bioconda -c conda-forge -c anaconda python=3 minimap2 scikit-learn kneed numpy samtools bcftools
```
If you face any issues on installing any dependencies via conda, try to install them directly and input their execute path into the program manually.

Otherwise, please manually install the following dependencies.
* scikit-learn
* kneed
* numpy
* minimap2
* bcftools
* samtools

## Usage
```
usage: gene_recomb.py [-h] -ref REF_FILE -read READ_FILE -out OUTDIR [-pm MINIMAP2] [-ps SAMTOOLS] [-pb BCFTOOLS]

optional arguments:
  -h, --help            show this help message and exit
  -ref REF_FILE, --ref_file REF_FILE
                        reference sequence, .fasta format
  -read READ_FILE, --read_file READ_FILE
                        read sequence, .fasta format
  -out OUTDIR, --out_directory OUTDIR
                        output directory
  -pm MINIMAP2, --path_minimap2 MINIMAP2
                        path to minimap2
  -ps SAMTOOLS, --path_samtools SAMTOOLS
                        path to samtools
  -pb BCFTOOLS, --path_bcftools BCFTOOLS
                        path to bcftools
```

## Output
* `<outdir>/bin_{base_aln}_{0-9|noise}.fasta` contains the clustered reads for each bin based on the sequencing divergence. The header line includes the read id, map quality and identity against consensus sequence, and a list a referenece assignments separated by `,` for each nucleotide, if several reference names are found with in one assignment separated by `:`, it represents the nucleotide may comes from multiple references.
* `<outdir>/consensus_bin_{base_aln}_{0-9|noise}.fasta` contains the consensus sequence generated from such group using bcftools.
* `<outdir>/work.log` the log file includes some configuration such as optimal k-mer size, number of clusters, DBSCAN paramter settings, and etc.

## Algorithm Pipeline
1. Fix the strandedness of the reads to consistent with the references based on minimap2 alignment (side-result: base alignment for each read).
2. Find an optimal k-mer size based on k-size vs unmapped rate trade (median unmapped rate of k-mers in a read among all the reads), using the cutoff error rate (default with 0.05).
3. Construct a hash table data structure: chunk the reference sequences into k-mers and merge identical k-mers with labels and keep its position from sequence.
4. For each read, do a sliding window approach using the optimal k value from (2) and decide its relation to references via identical matches for every nucleotide.
5. Use base alignment from step 1 to primarily cluster the reads.
6. For each base cluster, generate a pairwise distance matrix using pairwise divergence distance, which is defined as the average of the exclusive-or difference for each mapped reference by two reads. Run DBSCALL clustering to form a secondary cluster results.
7. For each secondary cluster, pick a read with highest identity to its base alignment reference, use minimap2 to generate a multiple alignment within the cluster and use bcftools to generate a consensus sequence, run minimap2 against such consensus sequence and record the alignment score.
8. Output the results by cluster.