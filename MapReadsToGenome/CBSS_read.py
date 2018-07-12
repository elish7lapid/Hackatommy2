import pandas as pd
import numpy as np


name_to_index ={"QNAME":0, "FLAG":1, "RNAME":2, "POS":3, "MAPQ":4, "CIGAR":5, "RNEXT":6, "PNEXT":7, "TLEN":8, "SEQ":9,
                                  "QUAL":10,
                                  "*":11, " DIFF":12, "METH":13, "**":14}


THRESHOLD_PRECENT = 0.25
class CBSS_read:

   def __init__(self, chr, pos, seq, meth):
       self.chr = chr
       self.pos = int(pos)
       self.seq = seq
       self.methylation = meth
       self.countMethylation()


   def countMethylation(self):
       num_z = self.methylation.count("z")
       num_Z = self.methylation.count("Z")
       if (num_z == 0 and num_Z == 0):
           self.methy_precentage = 0
       else:
           self.methy_precentage =num_Z /(num_z + num_Z  ) # new field)
       self.methy_count = num_Z
       self.unmethy_count =num_z   # new field


def create_reads_vector(filename):
    file = open(filename)
    # table_file = pd.read_table(filename,header = None, sep = "\t", index_col=False)
    # total_rows = table_file.shape[0]


    meth_SAM_file = open(r"D:\LEA\BIOINFORMATICS\Year_3\Sem_B\TOMMY2\HAckatommy\methylated_SAM", "w")
    unmeth_SAM_file = open(r"D:\LEA\BIOINFORMATICS\Year_3\Sem_B\TOMMY2\HAckatommy\unmethylated_SAM", "w")
    reads_list = []
    lines = file.readlines()
    for line in lines:
        parsed = line.split()
        chr = parsed[name_to_index["RNAME"]]
        pos = int(parsed[name_to_index["POS"]])
        seq = parsed[name_to_index["SEQ"]]
        meth = parsed[name_to_index["METH"]]
        cbss_Read = CBSS_read(chr, pos, seq, meth)
        cbss_Read.countMethylation()

        end = pos+len(seq)
        if cbss_Read.methy_precentage <THRESHOLD_PRECENT:
            unmeth_SAM_file.write(line)
        else:
            meth_SAM_file.write(line)




        reads_list.append(cbss_Read)

    reads_vector = np.array(reads_list)
    return reads_vector

create_reads_vector(r"D:\LEA\BIOINFORMATICS\Year_3\Sem_B\TOMMY2\HAckatommy\cat_H3K36me3_2D.sam_ch9.unique")