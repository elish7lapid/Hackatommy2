
import numpy as np

BIN_SIZE = 10
READ_SIZE = 100
BED_GRAPH_HEADER = "chrom\tstart\tend\tvalue\n"

def create_chrom_coverage_vec(positions_vec, chrom_length):
	"""
	For a given chromosome, this method returns for each binned
	position in the chromosome, how many reads got mapped to it. The
	method assumes the reads are of fixed sixe: READ_SIZE.
	To change bin and read sizes, change the global constants
	READ_SIZE and BIN_SIZE.
	:param positions_vec: A vector of start position of reads over
	the chromosome
	:param chrom_length: The length of the chromosome
	:return: a vector (of length chrom_length//BIN_SIZE) counting how
	many reads got mapped to each binned-position in the chromosome.
	"""
	scaled_read_size = READ_SIZE // BIN_SIZE
	scaled_positions = np.array(positions_vec) // BIN_SIZE
	coverage_vec = np.zeros((chrom_length//BIN_SIZE))
	
	read_positions = np.tile(np.arange(scaled_read_size),
	                         [len(scaled_positions)]) + \
	                 np.repeat(scaled_positions, scaled_read_size)
	unique, counts = np.unique(read_positions, return_counts=True)
	coverage_vec[unique] += counts
	return coverage_vec
	
	
def output_coverage_as_bedGraph(chrom_name, coverage_vec):
	bed_graph_file = open("chrom_coverage_bedGraph.txt", "w")
	bed_graph_file.write(BED_GRAPH_HEADER)
	for i in range(len(coverage_vec)):
		if coverage_vec[i] > 0:
			bed_graph_file.write("{}\t{}\t{}\t{}\n".format(chrom_name,
			                                               i*BIN_SIZE,
			                                               i*BIN_SIZE +
			                                               BIN_SIZE,
			                                               coverage_vec[i]))
	bed_graph_file.close()
	return


def preform_coverage_analysis(chrom_name, chrom_length,
                              positions_vec):
	coverage_vec = create_chrom_coverage_vec(positions_vec,
	                                         chrom_length)
	output_coverage_as_bedGraph(chrom_name, coverage_vec)
	