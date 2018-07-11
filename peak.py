import pysam
import numpy as np
import sys

CHROMO_IDEX = 0
START_INDEX = 1
END_INDEX=2
LENGTH_INDEX = 3
SMMIT_POS_INDEX = 4

HEIGHT_INDEX = 8

class Peak:

   def _init_(self, chr, start, end, length, summit_pos,summit_height, index):
       self.chr = chr
       self.start = int(start)
       self.end = int(end)
       self.length = length
       self.summit_pos = summit_pos
       self.summit_height = summit_height
       self.reads = []
       # self.index = index
       self.name = chr + "_" + str(index)
       self.methylation_percentage = 0


def create_peak_vector(filename):
    file = open(filename)
    index = 0
    peaks_list = []
    lines = file.readlines()
    for i in range(24, len(lines)):
        line = lines[i]
        parsed = line.split()
        chr = parsed[CHROMO_IDEX]
        start = parsed[START_INDEX]
        end = parsed[END_INDEX]
        length = parsed[LENGTH_INDEX]
        summit_pos = parsed[SMMIT_POS_INDEX]
        summit_height = parsed[HEIGHT_INDEX]

        cbss_Read = Peak(chr, start, end, length, summit_pos,summit_height, index)
        index += 1
        peaks_list.append(cbss_Read)

    peaks_vector = np.array(peaks_list)
    return peaks_vector





def main():
    if len(sys.argv )!= 2:
        print("USAGE : file name containing the peaks (xls)")
    filename = sys.argv[1]

    peal_list = create_peak_vector(filename)
    samfile = pysam.AlignmentFile("ex1.bam", "rb")
    for peak in peal_list:
        for read in samfile.fetch('chr1', peak.start, peak.end):
            peak.reads.append[read]

    samfile.close()


if __name__ == '__main__':
    main()