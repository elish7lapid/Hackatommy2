import pandas as pd
import numpy as np
import sys



name_to_index ={"QNAME":0, "FLAG":1, "RNAME":2, "POS":3, "MAPQ":4, "CIGAR":5, "RNEXT":6, "PNEXT":7, "TLEN":8, "SEQ":9,
                                  "QUAL":10,
                                  "*":11, " DIFF":12, "METH":13, "**":14}
class CBSS_read:

   def __init__(self, chr, pos, seq, meth):
       self.chr = chr
       self.pos = pos
       self.seq = seq
       self.methylation = meth



   def countMethylation(self):
       self.methy_precentage = self.methylation.count("Z") /self.methylation.count("z)")   # new field
       self.methy_count = self.methylation.count("Z")
       self.unmethy_count = self.methylation.count("z)")   # new field




def main():
    if len(sys.argv )!= 2:
        print("USAGE : file name containing the reads")


    filename = sys.argv[1]


    file = open(filename)


    # table_file = pd.read_table(filename,header = None, sep = "\t", index_col=False)
    # total_rows = table_file.shape[0]

    reads_list = []
    lines = file.readlines()
    for line in lines:
        parsed = line.split()
        chr = parsed[name_to_index["RNAME"]]
        pos = parsed[name_to_index["POS"]]
        seq = parsed[name_to_index["SEQ"]]
        meth = parsed[name_to_index["METH"]]
        cbss_Read = CBSS_read(chr, pos, seq, meth)
        reads_list.append(cbss_Read)

    reads_vector = np.array(reads_list)
    return reads_vector


if __name__ == '__main__':
    main()