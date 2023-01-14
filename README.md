# Genome Recombination Identification (Under development)
This program aims to discover the relation between a recombined genome sequence versus multiple reference sequences. Currently the program is under development and may containing unexpected errors.


## Usage
The program requires [minimap2](https://github.com/lh3/minimap2) for preprocessing, either modify the `minimap_exec` variable stored in `src/utils/utils.py` to specified the location of `minimap2` or ensure it can be directly accessed using `minimap2` command without absolute path.

```
python3 src/gene_recomb.py <ref_file.fasta> <read_file.fastq> <outdir>
```

## Output
* `<outdir>/bin{0_9}.fasta` contains the clustered reads for each bin based on the sequencing divergence. The header line includes the read id and a list a referenece assignments separated by `,` for each nucleotide, if several reference names are found with in one assignment, it represents the nucleotide may comes from multiple references.
* `<outdir>/work.log` the log file includes some configuration such as optimal k-mer size, number of clusters, and etc.

## Algorithm Pipeline
1. Fix the strandedness of the reads to consistent with the references based on minimap2 alignment.
2. Find an optimal k-mer size based on k-size vs unmapped rate trade (median unmapped rate of k-mers in a read among all the reads), using the cutoff error rate (default with 0.05).
3. Construct a hash table data structure: chunk the reference sequences into k-mers and merge identical k-mers with labels and keep its position from sequence.
4. For each read, do a sliding window approach using the optimal k value from (2) and decide its relation to references via identical matches for every nucleotide.
5. Iteratively binning the reads into clusters based on the divergence rate (default as 0.30) between two reads, two reads be classified into same cluster if they have divergence score within the cutoff bound, the divergence rate is calculated as the average of the exclusive-or difference for each mapped reference by two reads.
