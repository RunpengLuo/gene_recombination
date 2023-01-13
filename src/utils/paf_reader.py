import os
import sys
import re

class Paf_reader:
        
    def __init__(self, paf_file: str) -> None:
        if not os.path.exists(paf_file):
            print("{0} not exists".format(paf_file))
            sys.exit(1)
        paf_fd = open(paf_file, "r")
        self.lines = []
        for line in paf_fd.readlines():
            if line == "\n":
                break
            paf_line = Paf_line(line)
            self.lines.append(paf_line)
        paf_fd.close()

        self.cigar_parser = re.compile("[0-9]+[MIDNSHP=X]")

    def get_lines(self):
        return self.lines

    def _parse_cigar(self, cigar_str, read_str, ctg_str, qstart, cstart):
        """ cigar parsing is borrowed from Flye assembler: flye/utils/sam_parser.py/ line 214

        Args:
            cigar_str (_type_): cigar alignment string
            read_str (_type_): query string
            ctg_str (_type_): reference string
            qstart (_type_): query starting position (0-based)
            cstart (_type_): reference starting position (0-based)

        Raises:
            AlignmentException: _description_

        Returns:
            _type_: _description_
        """
        #ctg_str = self.ref_fasta[ctg_name]
        trg_seq = []
        qry_seq = []
        trg_start = cstart
        trg_pos = cstart
        qry_start = qstart
        qry_pos = qstart

        left_hard = True
        left_soft = True
        hard_clipped_left = 0
        hard_clipped_right = 0
        soft_clipped_left = 0
        soft_clipped_right = 0
        for token in self.cigar_parser.findall(cigar_str):
            size, op = int(token[:-1]), token[-1:]
            if op == "H":
                if left_hard:
                    qry_start += size
                    hard_clipped_left += size
                else:
                    hard_clipped_right += size
            elif op == "S":
                qry_pos += size
                if left_soft:
                    soft_clipped_left += size
                else:
                    soft_clipped_right += size
            elif op in "MX=":
                qry_seq.append(read_str[qry_pos : qry_pos + size].upper())
                trg_seq.append(ctg_str[trg_pos : trg_pos + size].upper())
                qry_pos += size
                trg_pos += size
            elif op == "I":
                qry_seq.append(read_str[qry_pos : qry_pos + size].upper())
                trg_seq.append("-" * size)
                qry_pos += size
            elif op == "D":
                qry_seq.append("-" * size)
                trg_seq.append(ctg_str[trg_pos : trg_pos + size].upper())
                trg_pos += size
            else:
                raise AlignmentException("Unsupported CIGAR operation: " + str(op))
            left_hard = False
            if op != "H":
                left_soft = False

        trg_seq = "".join(trg_seq)
        qry_seq = "".join(qry_seq)
        matches = 0
        for i in range(len(trg_seq)):
            if trg_seq[i] == qry_seq[i]:
                matches += 1
        err_rate = 1 - matches / len(trg_seq)

        trg_end = trg_pos
        qry_end = qry_pos + hard_clipped_left
        qry_len = qry_end + hard_clipped_right
        qry_start += soft_clipped_left
        qry_end -= soft_clipped_right

        return (trg_start, trg_end, len(ctg_str), trg_seq,
                qry_start, qry_end, qry_len, qry_seq, err_rate)

class Paf_line:
    def __init__(self, line: str) -> None:
        items = line[:-1].split()
        self.qname = items[0]
        self.qlen = int(items[1])
        self.qstart = int(items[2])
        self.qend = int(items[3])

        self.is_forward = True if items[4] == '+' else False

        self.rname = items[5]
        self.rlen = int(items[6])
        self.rstart = int(items[7])
        self.rend = int(items[8])

        self.nmatch = int(items[9])
        self.nregion = int(items[10])
        self.mapq = int(items[11])

        self.tags = dict()
        for tag in items[12:]:
            tid, tt, tval = tag.split(":")
            if tt == "i":
                self.tags[tid] = int(tval)
            else:
                self.tags[tid] = str(tval)
        return

class AlignmentException(Exception):
    pass