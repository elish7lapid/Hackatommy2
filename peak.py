import pysam
import numpy as np
import sys
from MapReadsToGenome import CBSS_read

CHROMO_IDEX = 0
START_INDEX = 1
END_INDEX=2
LENGTH_INDEX = 3
SMMIT_POS_INDEX = 4

HEIGHT_INDEX = 8

H_27_4D = "/mnt/5D37B39465DF0588/TEMP/hack/cat_H3K36me3_4D.bam.unique.bam"
WCE_4D =


class Peak:

   def __init__(self, chr, start, end, length, summit_pos,summit_height, index):
       self.chr = chr
       self.start = int(start)
       self.end = int(end)
       self.length = length
       self.summit_pos = summit_pos
       self.summit_height = summit_height
       self.reads_from_chip = []
       self.reads_from_wce = []
       # self.index = index
       self.name = chr + "_" + str(index)
       self.methylation_percentage = 0
       self.methylation_percentage_wce = 0

   def methylation_percent_in_peaks(self):
        """
        :param peaks_list: list of all peaks objects
        :return: null
        """
        counter_methylated = 0
        counter_non_methylated = 0
        for read in self.reads_from_chip:
            counter_methylated += read.methy_count
            counter_non_methylated += read.unmethy_count
        self.methyaltion_percent = counter_methylated / (
        counter_methylated + counter_non_methylated)

        for read in self.reads_from_wce:
            counter_methylated += read.methy_count
            counter_non_methylated += read.unmethy_count
        self.methyaltion_percent_wce = counter_methylated / (
        counter_methylated + counter_non_methylated)

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
    if len(sys.argv )!= 3:
        print("USAGE : file name containing the peaks locations (xls)")

    # File with the peaks locations from MACS in a BED format
    # /mnt/5D37B39465DF0588/TEMP/hack/NA_peaks_short.xls
    filename = sys.argv[1]

    peak_list = create_peak_vector(filename)[:2]
    samfile_chipseq = pysam.AlignmentFile(H_27_4D, "rb")
    samfile_WCE = pysam.AlignmentFile(WCE_4D, "rb")

    for peak in peak_list:
        print( peak.start)
        print(peak.end)
        for read in samfile_chipseq.fetch(peak.chr, peak.start, peak.end):
            read_arr = str(read).split()
            pos = read_arr[3]
            seq = read_arr[9]
            meth = read_arr[-5]
            # print(pos, seq, meth)
            read_obj = CBSS_read.CBSS_read(peak.chr, pos, seq, meth)
            # print(read_obj.methylation, read_obj.methy_count)
            peak.reads_from_chip.append(read_obj)

        for read in samfile_WCE.fetch(peak.chr, peak.start, peak.end):
            read_arr = str(read).split()
            pos = read_arr[3]
            seq = read_arr[9]
            meth = read_arr[-5]
            # print(pos, seq, meth)
            read_obj2 = CBSS_read.CBSS_read(peak.chr, pos, seq, meth)
            # print(read_obj.methylation, read_obj.methy_count)
            peak.reads_from_wce.append(read_obj2)

        peak.methylation_percent_in_peaks()
        print("percentage: ", peak.methyaltion_percent)

    samfile_chipseq.close()
    samfile_WCE.close()

if __name__ == '__main__':
    main()