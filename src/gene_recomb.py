import os, sys
import argparse

from utils.file_parser import *
from utils.utils import *

from utils.clustering import *
from utils.consensus import get_consensus, aln_consensus
import utils.external as external

# example: python src/gene_recomb.py reference.fasta 20220413_XX.Sample01.BC48.Qmin15.min2500max3500.fq

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog=sys.argv[0])
    parser.add_argument("-ref", "--ref_file", dest="ref_file", required=True, default=None, 
        type=str, help="reference sequence, .fasta format")
    
    parser.add_argument("-read", "--read_file", dest="read_file", required=True, default=None, 
        type=str, help="read sequence, .fasta format")   

    parser.add_argument("-out", "--out_directory", dest="outdir", required=True, default="outdir/", 
        type=str, help="output directory")
    
    parser.add_argument("-pm", "--path_minimap2", dest="minimap2", required=False, default="minimap2", 
        type=str, help="path to minimap2")
    
    parser.add_argument("-ps", "--path_samtools", dest="samtools", required=False, default="samtools", 
        type=str, help="path to samtools")
    
    parser.add_argument("-pb", "--path_bcftools", dest="bcftools", required=False, default="bcftools", 
        type=str, help="path to bcftools")
    
    args = parser.parse_args()

    if os.path.exists(args.outdir):
        print("output directory exists, please remove it first or use another directory name to avoid overwrite")
        sys.exit(0)
    
    if args.outdir[-1] != '/':
        args.outdir += '/'
    os.system("mkdir {0}".format(args.outdir))

    external.minimap_exec = args.minimap2
    external.samtools_exec = args.samtools
    # external.cd_hit_exec = "cd-hit"
    external.bcftools_exec = args.bcftools
    external.outdir = args.outdir

    # parsing files
    read_dict = get_read_file(args.read_file)
    ref_dict = get_ref_file(args.ref_file)

    # Step 1. fix reads strandedness
    base_alignments, paf_reader = get_base_alignment(read_dict, args.read_file, args.ref_file)

    # Step 2. compute the optimal k value based on unmapped rate
    final_k, reftree, read_dist, read_pos = get_best_k(ref_dict, read_dict)
    print("optimal k selection: ", final_k)

    # Step 3. project the region-reference selection to each of the reference
    # ref_projection = {}
    ref_ids = list(ref_dict.keys())
    read_feature_dict = {} # used for clustering
    for read_id in read_dist.keys():
        ref_assignments_list, ref_assignments_int = read_projection(read_dist, read_pos, ref_dict, read_id, ref_ids)
        # ref_projection[read_id] = ref_assignments_list
        read_feature_dict[read_id] = ref_assignments_int

    # base binning
    print("Base-alignment clustering")
    base_bins = {}
    for read_id, aln_obj in base_alignments.items():
        if aln_obj.rname not in base_bins:
            base_bins[aln_obj.rname] = []
        base_bins[aln_obj.rname].append(read_id)

    print("Number of base bins: ", len(base_bins))
    print("DBSCALL clustering on each base bin..")
    sec_bins = {}
    dbscan_config = {}
    for base_id, base_bin in base_bins.items():
        print("Constructing pairwise-distance matrix for base bin: ", base_id)
        labels, n_clusters_, n_noise_,  min_pts, eps = dbscan_clustering(get_distance_matrix(read_feature_dict, base_bin))
        sec_dict = {}
        for i, sec_id in enumerate(labels):
            if sec_id not in sec_dict:
                sec_dict[sec_id] = []
            sec_dict[sec_id].append(base_bin[i])
        for sec_id, sec_bin in sec_dict.items():
            sec_bins["bin_{0}_{1}".format(base_id, (str(sec_id) if sec_id != -1 else 'noise'))] = sec_bin
        dbscan_config[base_id] = (min_pts, eps)

    print("Number of secondary bins: ", len(sec_bins))

    # construct consensus sequence for each cluster
    print("Constructing consensus sequences..")
    aln_dict = {}
    for sec_id, sec_bin in sec_bins.items():
        con_fname = get_consensus(sec_id, read_dict, sec_bin, base_alignments)
        aln_consensus(aln_dict, read_dict, sec_bin, con_fname)

    # Store the result
    output_result_by_cluster(read_dict, aln_dict, read_dist, sec_bins)

    config_file = external.outdir + "work.log"
    os.system("echo "" > " + config_file)
    with open(config_file, "w") as cfd:
        cfd.write("k-mer size: {0}\nNumber of bins: {1}\nNumber of reads: {2}\nNumber of references: {3}\n".format(
            final_k, len(sec_bins), len(read_dict), len(ref_dict)
        ))
        cfd.write(">DBSCAN clustering: base_id\tmin_pts\teps\n")
        for base_id, (min_pts, eps) in dbscan_config.items():
            cfd.write("{0}\t{1}\t{2}\n".format(base_id, min_pts, eps))
        cfd.write("Bin_id\tNumber_of_reads\tRead_ids\n")
        for bid, sec_bin in sec_bins.items():
            cfd.write("{0}\t{1}\t{2}\n".format(bid, len(sec_bin), ', '.join(sec_bin)))
        cfd.close()

    os.system("rm {0}/temp*".format(external.outdir))
    sys.exit(0)