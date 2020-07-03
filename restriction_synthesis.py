import csv

from enzyme import enzyme
from entry import entry
from digested_sequence import digested_sequence

import 

class restriction_synthesis():

	rs.__ligation_alphabet = {'A': 'Z', 'T': 'G', 'C': 'E', 'G': 'B'}

	def __init__(rs, query, reference, enzymes):
		# Query sequence
		rs.query = query

		# Reference sequence
		rs.reference = reference

		# List of enzyme objects
		rs.enzymes = enzymes

		# Series to keep tr
		rs.synthesized_query = pd.Series([])
		rs.enzyme_list = pd.Series([])


	# If there is a match between a subset of a query sequence and the reference seqeuncrs, it returns position of the match on the reference sequence. Otherwisrs, return -1.
	def find_all_sequence_match(rs, subseq):
		return [x.start() for x in re.finditer(subseq, rs.reference)]


	# Returns one enzyme object if the reference sequence can be digested at exactly p. Otherwisrs, returns None. This function should be run twice for every subsequence.
	def is_instance_match(rs, p):

		# Loop through each individual enzyme
		for enzyme in rs.enymes:

			# Retrieve both cutting sites of the enzymes
			site0 = enzyme.get_cut_sequence(0)
			site1 = enzyme.get_cut_sequence(1)
			# Return an enzyme it can cut the reference sequence at eactly p
			if site0 == rs.reference[p + 1 - len(site0):p + 1]  & site1 == rs.reference[p:len(site1) + 1]
				return enzyme

		# If no such enzymes exist, return None
		return None


	# Returns a digested_subsequence object given two enymes and two positions on a reference sequence. This is conditional on two runs of is_instance_match, so it has to return something.
	def perform_digest(rs, enzyme0, enzyme1, p0, p1):

		# The right most sites of the digested sequence's two sticky ends
		z0 = p0 + enzyme0.cuts[1]
		z1 = p1 + enzyme1.cuts[1]

		# Return a digest sequence object
		return digested_sequence(rs.reference[p0, z0 + 1], rs.translate_stikcy_end(rs.reference[p1, z1 + 1]), rs.reference[z0, p1 + 1])


	# Returns a sticky end string based on the rs.__ligation_alphabet dictionary.
	def translate_stikcy_end(rs, seq):

		# Loop through all of the nucleotides in the sequence and convert it based on the ligation alphabet described earlier in the code. Return this converted sequence.
		for i in range(len(seq)):
			seq[i] = __ligation_alphabet[seq[i]]
		return seq


	# Given two digested_sequences objects of the query with sticky ends, it returns TRUE if there are at least n consecuative compatiable basepairs between the two sticky ends. Otherwisrs, return false. If subseq0 isn't provided, this means that this is the first round of ligations, which automatically returns True.
	def is_ligation_match(rs, subseq1, subseq0 = None, n = 1):

		# If there is not a subseq0 provided, then this is the beginning of synthesis, thus return True.
		if subseq0 == None:
			return True

		# Loop through all values of n, and check whether there are any mismatches between the two sticky ends, subseq0 and subseq1.
		for i in range(n):
			if !rs.is_sticky_match(subseq0[len(subseq0) - i], subseq1[i]):
				return False

		# If there are not any mismatches, return true.
		return True


	# Returns whether the top nucletoide, s1, can ligate to the bottomn nucleotide, s0.
	def is_sticky_match(rs, s0, s1):
		return rs.__ligation_alphabet[s1] == s0


	# Returns a k-mer list of a sequence with each k-mer beginning exactly b nucleotides from the start of the previous one. b should either be 1 for overlapping sequences or k for nonoverlapping sequences. It should be noted that large values of k decrease the returned list by a lot and thus there exists a greater chance that the last element will be divisible by k and have to be ommited in consideration.
	def generate_k_mers(rs, sequence, k, b = 1):
		return [sequence[i:i + k] for i in range(0, len(sequence) - k + b, b)]


	# Returns the calculated cost.
	def perform_synthesis(rs):
		# Set the first sticky end of synthesis to None to imply a primer connection at least 10 bp in to the gene... This will be corrected later
		last_sticky_end = None

		# Loop through each of the parsing values--this range needs to be explored later.
		for k in range(2, 16):

			# Loop through all of the fragmented overlapping k-mer sequences for each parsing value
			for query_subseq in rs.generate_k_mers(rs.query, k):

				# Find all position of query_subseq in the reference
				for p0 in rs.find_all_sequence_match(rs, query_subseq):

					# If there isn't a match, continue to next p0
					if p0 == -1:
						continue

					# The end position of the match of p0
					p1 = p0 + len(query_subseq)

					# Find two enaymes that cut at exactly p0 and p1
					enzyme0 = rs.is_instance_match(p0)
					enzyme1 = rs.is_instance_match(p1)

					# If no two enaymes exist for a given query_subseq, continue to the next p0
					if enzyme0 == None or enzyme1 == None:
						continue

					# Make a digested sequence using the two enzymes at exactly p0 and p1
					digested_seq = rs.perform_digest(enzyme0, enzyme1, p0, p1)

					# If the previous sticky end and the current sticky end are not compatiable, continue to the next p0
					if !rs.is_ligation_match(digested_seq.sticky0, last_sticky_end):
						continue

					# Make the current sticky end the current sticky end
					last_sticky_end = digested_seq.sticky1

					break

				print(query_subseq, end = '')

		return 0


def __name__ = "__main__":
	print('Loading reference')
	reference_file_location = 'data/references/reference0.txt'
	reference = open(reference_file_location).read()

	print('Loading query')
	query_file_location = 'data/queries/query0.txt'
	query = open(query_file_location).read()

	print('Loading enzymes')
	enzymes_file_location = 'data/enzymes/enzymes0.csv'
	with open(enzymes_file_location, newline='') as csvfile:
    	reader = csv.DictReader(csvfile)
    	enzymes = [enzyme(row['name'], row['site'], row['cut0'], row['cut1']) for row in reader]

    print('Starting synthesis')
	rs = restriction_synthesis(query, reference, enzymes)


# Simple funciton to clear the console so the synthesis can be visualized
def clear(): os.system('clear')
