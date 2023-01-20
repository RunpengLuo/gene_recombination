import os, sys
from utils.utils import run_minimap
import utils.external as external
from utils.paf_reader import Paf_reader

from Bio import AlignIO
from Bio.Align import AlignInfo


def get_consensus_v2(bid: str, glb_read_dict: dict, read_ids: list):
    gfname = external.outdir + "temp_reads.fasta"
    os.system("echo "" > {0}".format(gfname))

    alignment = AlignIO.read(sys.argv[1], 'fasta')
    summary_align = AlignInfo.SummaryInfo(alignment)
    summary_align.dumb_consensus(float(sys.argv[2]))
    return

def get_consensus(bid: str, glb_read_dict: dict, read_ids: list, base_alignments: dict):
    rfname = external.outdir + "temp_ref.fasta"
    gfname = external.outdir + "temp_reads.fasta"
    sam_fname = external.outdir + "temp_aln.sam"
    bam_fname = external.outdir + "temp_aln.sam"
    vcf_fname = external.outdir + "temp_calls.vcf.gz"

    # init files
    con_fname = external.outdir + "consensus_{0}.fasta".format(bid)
    for fname in [rfname, gfname, sam_fname, bam_fname, vcf_fname, con_fname]:
        os.system("echo "" > {0}".format(fname))

    # pick the reference seq as consensus base
    con_id = ""
    con_read = ""
    best_aln = sys.maxsize

    gfd = open(gfname, "w")
    for read_id in read_ids:
        if base_alignments[read_id].tags['NM'] < best_aln:
            best_aln = base_alignments[read_id].tags['NM']
            con_id = read_id
            con_read = glb_read_dict[read_id]
        gfd.write(">{0}\n{1}\n".format(read_id, glb_read_dict[read_id]))
    gfd.close()

    rfd = open(rfname, "w")
    rfd.write(">{0}\n{1}\n".format(con_id, con_read))
    rfd.close()

    # get alignment informations
    run_minimap(rfname, gfname, sam_fname, flags="-a -D")
    os.system("{0} view -bS {1} | {0} sort - -o {2}".format(external.samtools_exec, sam_fname, bam_fname))
   
    # variant calling using bcftools
    os.system("{0} mpileup -X ont -Ou -f {1} {2} | {0} call -Ou -mv | {0} norm -f {1} -Oz -o {3}".format(
        external.bcftools_exec, rfname, bam_fname, vcf_fname)
    )
    
    os.system("{0} index {1}".format(external.bcftools_exec, vcf_fname))
    os.system("{0} consensus -f {1} {2} > {3}".format(external.bcftools_exec, rfname, vcf_fname, con_fname))

    os.system("rm {0}/temp*".format(external.outdir))
    return con_fname

def aln_consensus(glb_aln_dict: dict, glb_read_dict: dict, read_ids: list, con_fname: str):
    gfname = external.outdir + "temp_reads.fasta"
    paf_fname = external.outdir + "temp_aln.paf"
    
    os.system("echo "" > {0}".format(gfname))
    gfd = open(gfname, "w")
    for read_id in read_ids:
        gfd.write(">{0}\n{1}\n".format(read_id, glb_read_dict[read_id]))
    gfd.close()

    run_minimap(con_fname, gfname, paf_fname)
    paf_reader = Paf_reader(paf_fname)
    for aln_obj in paf_reader.get_lines():
        glb_aln_dict[aln_obj.qname] = (aln_obj.mapq, round(100*aln_obj.nmatch/float(aln_obj.nregion),2))

    return