import pysam
import numpy as np
import sys
from MapReadsToGenome import Read
from scipy.stats import binom

CHROMO_IDEX = 0
START_INDEX = 1
END_INDEX=2
LENGTH_INDEX = 3
SMMIT_POS_INDEX = 4
HEIGHT_INDEX = 8


N_INDEX = 0
P_INDEX = 1

peak_number = 0

wce_N_p = []
chip_Z = []
pval_of_peak = []

# H_27_4D = "/cs/phd/alonap/www/hackathon2018/ChIP-BS-seq/cat_H3K27ac_4D.bam.unique.bam"
# # H_27_4D = "/mnt/5D37B39465DF0588/TEMP/hack/cat_H3K36me3_4D.bam.unique.bam"
# WCE_4D = "/cs/phd/alonap/www/hackathon2018/ChIP-BS-seq/cat_WCE_4D.bam.unique.bam"
# #arg: python3 Hackatommy2/peak.py "peaks H3K27ac/NA_peaks.xls

# cs/phd/alonap/www/hackathon2018/ChIP-BS-seq/

# for K4ME3
# H_27_4D = "/cs/phd/alonap/www/hackathon2018/ChIP-BS-seq/cat_H3K36me3_4D.bam.unique.bam"
# WCE_4D = "/cs/phd/alonap/www/hackathon2018/ChIP-BS-seq/cat_WCE_4D.bam.unique.bam"
# # arg: python3 Hackatommy2/peak.py "peaks H3K36me3/NA_peaks.xls

# for stem cell
# H_27_4D = "/cs/phd/alonap/www/hackathon2018/ChIP-BS-seq/cat_H3K27ac_ES.bam.unique.bam"
# WCE_4D = "/cs/phd/alonap/www/hackathon2018/ChIP-BS-seq/cat_WCE_ES.bam.unique.bam"
#arg: python3 Hackatommy2/peak.py "ES/peaks H3K27ac/NA_peaks.xls" todo

H_27_4D = "/cs/phd/alonap/www/hackathon2018/ChIP-BS-seq/cat_H3K36me3_ES.bam.unique.bam"
WCE_4D = "/cs/phd/alonap/www/hackathon2018/ChIP-BS-seq/cat_WCE_ES.bam.unique.bam"
#arg: python3 Hackatommy2/peak.py "ES/peaks H3K36me3/NA_peaks.xls" todo


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
       self.index = index
       self.methylation_percentage = 0
       self.methylation_percentage_wce = 0


   def methylation_percent_in_peaks(self):
        """
        :param peaks_list: list of all peaks objects
        :return: null
        """
        # with open("peak_meth_percent_wce", "at") as writerWCE, open("peak_meth_percent_H3K27ac", "at") as writerChIP:
        counter_methylated = 0
        counter_non_methylated = 0
        for read in self.reads_from_chip:
            counter_methylated += read.methy_count
            counter_non_methylated += read.unmethy_count
        if (counter_methylated + counter_non_methylated) != 0:
            self.methylation_percentage = counter_methylated / (counter_methylated + counter_non_methylated)
        # writerChIP.write(self.name + " " + str(self.methylation_percentage) + "\n")
        chip_Z.append(counter_methylated)

        for read in self.reads_from_wce:
            counter_methylated += read.methy_count
            counter_non_methylated += read.unmethy_count
        if (counter_methylated + counter_non_methylated) != 0:
            self.methylation_percentage_wce = counter_methylated / (counter_methylated + counter_non_methylated)
        # writerWCE.write(self.name + " " + str(self.methylation_percentage_wce) + "\n")
        wce_N_p.append((counter_methylated + counter_non_methylated, self.methylation_percentage_wce))

   def write_to_file(self):
        # with open("peak_meth_percent_wce_ES", "at") as writerWCE, open("peak_meth_percent_H3K27ac_ES", "at") as writerChIP:
        with open("peak_meth_percent_wce_ES_36", "at") as writerWCE, open("peak_meth_percent_H3K36me3_ES", "at") as writerChIP:
        # with open("peak_meth_percent_wce_36", "at") as writerWCE, open("peak_meth_percent_H3K36me3", "at") as writerChIP:
            writerChIP.write(self.name + " " + str(self.methylation_percentage) + " " + str(pval_of_peak[self.index]) + "\n")
            writerWCE.write(self.name + " " + str(self.methylation_percentage_wce) + "\n")

def create_peak_vector(filename):
    global peak_number
    file = open(filename)
    index = 0
    peaks_list = []
    lines = file.readlines()
    for i in range(24, len(lines)):
    # for i in range(24, 27): #todo
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
    peak_number = len(peaks_list)
    return peaks_vector

def calc_pval():
    global wce_N_p, chip_Z, peak_number, pval_of_peak
    print("number: ", peak_number)

    N = [x[0] for x in wce_N_p]
    p = [x[1] for x in wce_N_p]
    # print("N: ", N)
    # print("p: ", p)
    wce_N_p = np.zeros((peak_number, 2))
    wce_N_p[:, N_INDEX] = np.array(N).reshape((peak_number))
    wce_N_p[:, P_INDEX] = p
    chip_Z = np.array(chip_Z)
    print("chip_Z: ", chip_Z)
    pval_of_peak = np.zeros((peak_number))

    # filename_chipk27 = "peak_meth_percent_H3K27ac"
    # filename_wce = "peak_meth_percent_wce"
    #
    # chip = np.loadtxt(filename_chipk27, dtype=str)
    # chip = chip[:, 1].astype(np.float)
    # wce = np.loadtxt(filename_wce, dtype=str)
    # p_vec_wce = wce[:, 1].astype(np.float)

    # the methylation percentage = the probability for methylation in each peak
    cdf_prob = binom.cdf(chip_Z, wce_N_p[:, N_INDEX], wce_N_p[:, P_INDEX])
    # print("cdf: ", cdf_prob)
    pval_of_peak = cdf_prob
    pval_of_peak[chip_Z >= (wce_N_p[:, N_INDEX]*wce_N_p[:, P_INDEX])] = 1 - cdf_prob[chip_Z >= (wce_N_p[:, N_INDEX]*wce_N_p[:, P_INDEX])]
    print("pval_of_peak:" , pval_of_peak)

def main():
    if len(sys.argv )!= 2:
        print("USAGE : file name containing the peaks locations (xls)")

    # File with the peaks locations from MACS in a BED format
    # /mnt/5D37B39465DF0588/TEMP/hack/NA_peaks_short.xls
    filename = sys.argv[1]

    peak_list = create_peak_vector(filename)

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
            read_obj = Read.Read(peak.chr, pos, seq, meth)
            # print(read_obj.methylation, read_obj.methy_count)
            peak.reads_from_chip.append(read_obj)

        for read in samfile_WCE.fetch(peak.chr, peak.start, peak.end):
            read_arr = str(read).split()
            pos = read_arr[3]
            seq = read_arr[9]
            meth = read_arr[-5]
            # print(pos, seq, meth)
            read_obj2 = Read.Read(peak.chr, pos, seq, meth)
            # print(read_obj.methylation, read_obj.methy_count)
            peak.reads_from_wce.append(read_obj2)

        peak.methylation_percent_in_peaks()
        # print("percentage: ", peak.methyaltion_percent)

    samfile_chipseq.close()
    samfile_WCE.close()
    calc_pval()

    for peak in peak_list:
        peak.write_to_file()


if __name__ == '__main__':
    main()
