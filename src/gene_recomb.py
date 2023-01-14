import os
import sys

from utils.file_parser import *
from utils.utils import *


# python gene_recomb.py reference.fasta 20220413_XX.Sample01.BC48.Qmin15.min2500max3500.fq

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("{0} <ref_file.fasta> <read_file.fastq> <outdir>".format(sys.argv[0]))
        sys.exit(0)

    ref_file = sys.argv[1]
    read_file = sys.argv[2]
    outdir = sys.argv[3]

    if os.path.exists(outdir):
        print("output directory exists, please remove it first or use another directory name to avoid overwrite")
        sys.exit(0)
    
    if outdir[-1] != '/':
        outdir += '/'
    os.system("mkdir {0}".format(outdir))

    # parsing files
    read_dict = get_read_file(read_file)
    ref_dict = get_ref_file(ref_file)

    # Step 1. fix reads strandedness
    base_alignment, paf_reader = get_base_alignment(read_dict, read_file, ref_file, outdir)

    # Step 2. compute the optimal k value based on unmapped rate
    final_k, reftree, read_dist, read_pos = get_best_k(ref_dict, read_dict)
    print("optimal k selection: ", final_k)

    # Step 3. project the region-reference selection to each of the reference
    ref_projection = {}
    for read_id in read_dist.keys():
        dist, poses = read_dist[read_id], read_pos[read_id]
        ref_assignments = {}
        for ref_id, ref_seq in ref_dict.items():
            ref_assignments[ref_id] = [0 for _ in range(len(ref_seq))]
        
        # project the result to 1d arr-es
        for (res, pos) in zip(dist, poses):
            if res != None:
                for (re, p) in zip(res, pos):
                    ref_assignments[re][p] = 1

        ref_projection[read_id] = ref_assignments

    # Step 4. bin the reads based on mapping divergence
    bins = read_binning(ref_projection, ref_dict.keys())

    # Store the result
    output_result_by_cluster(read_dict, read_dist, bins, outdir)

    config_file = outdir + "work.log"
    os.system("echo "" > " + config_file)
    with open(config_file, "w") as cfd:
        cfd.write("k-mer size: {0}\nNumber of bins: {1}\nNumber of reads: {2}\nNumber of references: {3}\n".format(
            final_k, len(bins), len(read_dict), len(ref_dict)
        ))
        cfd.write("Bin id\tNumber of reads\n")
        for bid, bin in enumerate(bins):
            cfd.write("{0}\t{1}\n".format(bid, len(bin)))
        cfd.close()

    sys.exit(0)