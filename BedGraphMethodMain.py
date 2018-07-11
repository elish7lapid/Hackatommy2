import CBSS_read
import ChromCoverage
import sys

CHROM_LENGTHS_MM9 = {"chr1": 197195432, "chr2": 181748087,
                     "chr3": 159599783, "chr4": 155630120,
                     "chr5": 152537259, "chr6": 149517037,
                     "chr7": 152524553, "chr8": 131738871,
                     "chr9": 124076172, "chr10": 129993255,
                     "chr11": 121843856, "chr12": 121257530,
                     "chr13": 120284312, "chr14": 125194864,
                     "chr15": 103494974, "chr16": 98319150,
                     "chr17": 95272651, "chr18": 90772031,
                     "chr19": 61342430, "chrX": 166650296,
                     "chrY": 15902555}

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("USAGE : file name containing the reads")
	filename = sys.argv[1]
	
	# Assuming all reads come from same chromosome:
	reads = CBSS_read.create_reads_vector(filename)
	positions = [read.pos for read in reads]
	ChromCoverage.preform_coverage_analysis(reads[0].chr,
	                                        CHROM_LENGTHS_MM9[reads[0].chr],
	                                        positions)