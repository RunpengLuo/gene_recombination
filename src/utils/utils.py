import os
from utils.paf_reader import *
from utils.transformer import *
from utils.ktree import *

minimap_exec = "minimap2"

def run_minimap(ref_file: str, query_file: str, out_file="aln.paf", xflag="map-ont", use_secondary=False, flags="-c"):
    print("Running minimap2 alignment..")
    os.system("{0} -x {1} --secondary={2} {3} {4} {5} > {6}".format(
        minimap_exec, xflag, "yes" if use_secondary else "no", flags, ref_file, query_file, out_file))
    print("Done")
    return

def get_base_alignment(read_dict: dict, read_file: str, ref_file: str, outdir):
    paf_file = outdir + "aln.paf"
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

def get_divergence_score(proj1dict: dict, proj2dict: dict, keys: list):
    total_len = 0
    divergence = 0
    for key in keys:
        div = 0
        ref_a, ref_b = proj1dict[key], proj2dict[key]
        total_len += len(ref_a)
        for (a,b) in zip(ref_a, ref_b):
            div += (a+b) % 2
        divergence += float(div)/len(ref_a)

    return divergence/len(keys)        

def read_binning(ref_projection: dict, ref_keys: list, div_rate=0.30):
    """Construct the read bins by iteratively adding reads to bins, create a new bin if none of existing bins has divergence rate
    less than <div_rate>, order invariant since computing xor for divergence check

    Args:
        ref_projection (dict): 1d projection on references for each reads
        div_rate (float, optional): divergence ratio. Defaults to 0.10.
    """
    read_count = 0

    bins = []
    bin_len = 0
    for rid, ref_proj in ref_projection.items():
        read_count += 1
        assign_bin = -1
        min_rate = sys.maxsize
        for bid in range(bin_len):
            added_rid = bins[bid][random.randint(0, len(bins[bid])-1)]
            rate = get_divergence_score(ref_projection[added_rid], ref_proj, ref_keys)
            if rate < min_rate:
                min_rate = rate
                assign_bin = bid
        
        if min_rate <= div_rate:
            bins[assign_bin].append(rid)
        else:
            # not added by any of the bins, create new bin 
            bins.append([rid])
            bin_len += 1

    print("Total clusters: ", bin_len)
    print("Total processed reads: ", read_count)

    return bins

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