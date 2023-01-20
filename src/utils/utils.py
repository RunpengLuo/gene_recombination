import os
from utils.paf_reader import *
from utils.transformer import *
from utils.ktree import *
import utils.external as external

def run_minimap(ref_file: str, query_file: str, out_file="aln.paf", xflag="map-ont", use_secondary=False, flags="-c"):
    print("Running minimap2 alignment..")
    os.system("{0} -x {1} --secondary={2} {3} {4} {5} > {6}".format(
        external.minimap_exec, xflag, "yes" if use_secondary else "no", flags, ref_file, query_file, out_file))
    print("Done")
    return

def get_base_alignment(read_dict: dict, read_file: str, ref_file: str):
    paf_file = external.outdir + "aln.paf"
    run_minimap(ref_file, read_file, paf_file, "map-ont", False, "-c")

    nfix = 0
    paf_reader = Paf_reader(paf_file)
    base_alignment = dict.fromkeys(read_dict.keys())
    for aln_obj in paf_reader.get_lines():
        # fix strandedness
        if not aln_obj.is_forward:
            read_dict[aln_obj.qname] = reverse_alpha(read_dict[aln_obj.qname])
            nfix += 1
        
        base_alignment[aln_obj.qname] = aln_obj

    print("number of fixed strandedness read: ", nfix, len(read_dict))
    return base_alignment, paf_reader

def read_projection(read_dist: dict, read_pos: dict, ref_dict: dict, read_id: str, ref_ids: list):
    # obtain reference-level alignments for a single read
    dist, poses = read_dist[read_id], read_pos[read_id]
    ref_assignments = {}
    for ref_id, ref_seq in ref_dict.items():
        ref_assignments[ref_id] = [0 for _ in range(len(ref_seq))]
    
    # project the result to 1d arr-es
    for (res, pos) in zip(dist, poses):
        if res != None:
            for (re, p) in zip(res, pos):
                ref_assignments[re][p] = 1

    ref_assignments_int = []
    for ref_id in ref_ids:
        int_str = [str(a) for a in ref_assignments[ref_id]]
        ref_assignments_int.append(''.join(int_str))

    return ref_assignments, int(''.join(ref_assignments_int), 2)

def merge_list(unmerged_arr: list):
    counter = 1
    prev = unmerged_arr[0]
    merged_arr = []
    for elem in unmerged_arr[1:]:
        if elem == prev:
            counter += 1
            prev = elem
        else:
            merged_arr.append((prev, counter))
            counter = 1
            prev = elem
    merged_arr.append((prev, counter))
    return merged_arr

def get_best_k(ref_dict: dict, read_dict: dict, err_rate=0.05, min_k=3, max_k=23):
   # evaluate current k choice
    k_unmapped = []
    k_single = []
    k_multi = []

    k_choices = [i for i in range(min_k, max_k, 1)]

    reftree = None
    read_ref_dist = None
    read_ref_pos = None

    final_k = (min_k + max_k) // 2

    for k in k_choices:
        med_unmapped = []
        med_single = []
        med_multi = []

        # construct reference table
        reftree = KTree(k)
        for ref_id, ref_seq in ref_dict.items():
            for ind in range(0, len(ref_seq) - k + 1):
                kmer = ref_seq[ind: ind + k]
                reftree.insert_by_val(kmer, ref_id, ind)

        read_ref_dist = {}
        read_ref_pos = {}
        # query the reads
        for read_id, read_seq in read_dict.items():
            dist = []
            pos = []
            nkmer = len(read_seq) - k + 1
            nunmap = 0
            nsingle = 0
            nmulti = 0
            for ind in range(0, nkmer):
                kmer = read_seq[ind: ind + k]
                if kmer not in reftree.key_dict:
                    dist.append(None)
                    nunmap += 1
                else:
                    dist.append(reftree.key_dict[kmer].dist)
                    pos.append(reftree.key_dict[kmer].pos)
                    if len(reftree.key_dict[kmer].dist) > 1:
                        nmulti += 1
                    else:
                        nsingle += 1
            read_ref_dist[read_id] = dist
            read_ref_pos[read_id] = pos
            med_unmapped.append(float(nunmap)/nkmer)
            med_single.append(float(nsingle)/nkmer)
            med_multi.append(float(nmulti)/nkmer)

        k_unmapped.append(median(med_unmapped))
        k_single.append(median(med_single))
        k_multi.append(median(med_multi))
        if median(med_unmapped) >= err_rate:
            final_k = k
            break

    # print("Total number of reads: ", len(read_dict))
    # print("K choices:", k_choices)
    # print("Unmapped trade: ", k_unmapped)
    # print("Single trade: ", k_single)
    # print("Double trade: ", k_multi)   
    return final_k, reftree, read_ref_dist, read_ref_pos