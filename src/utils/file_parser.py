import os

def get_ref_file(ref_file: str):
    ref_dict = {}
    with open(ref_file, "r") as ref_fd:
        i = 0
        info = ""
        ref_id = None
        ref_seq = None
        for line in ref_fd.readlines():
            if line == "\n":
                break
            info += line
            i = (i+1) % 2
            if i == 0:
                ref_id = str((info.split("\n")[0]).split()[0][1:])
                ref_seq = str(info.split("\n")[1])
                ref_dict[ref_id] = ref_seq
                info =  ""
        ref_fd.close()
    return ref_dict


def get_read_file(read_file: str):
    read_dict = {}
    with open(read_file, "r") as read_fd:
        i = 0
        read = ""
        read_id = None
        read_seq = None
        for line in read_fd.readlines():
            if line == "\n":
                break
            read += line
            i = (i+1) % 4
            if i == 0:
                # previous read is ready
                # get first no-space read id, ignore symbol @
                read_id = str((read.split("\n")[0]).split()[0][1:])
                read_seq = str(read.split("\n")[1])
                read_dict[read_id] = read_seq
                # idx2read.append(read.split("\n")[1])
                read = ""
        read_fd.close()
    return read_dict

def output_result_by_cluster(read_dict: dict, read_dist: dict, bins: list, outdir):
    for bid, bin in enumerate(bins):
        outfile = outdir + "bin{0}.fasta".format(bid)
        os.system("echo "" > {0}".format(outfile))
        with open(outfile, "w") as fd:
            for read_id in bin:
                dist = read_dist[read_id]
                out_st = ""
                for d in dist:
                    if d == None or len(d) == 0:
                        out_st += "*,"
                    else:
                        for elem in set(d):
                            out_st += str(elem) + ":"
                        out_st = out_st[:-1] + ","
                fd.write(">{0}\t{1}\n{2}\n".format(read_id, out_st, read_dict[read_id]))
            fd.close()
    print("result stored in: "+outdir+"bin\{#\}.fasta")
    return