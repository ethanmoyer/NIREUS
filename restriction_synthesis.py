
class restriction_synthesis():
	def __init__(rs, query, reference, enzymes):
		# Query sequence
		rs.query = query

		# Reference sequence
		rs.reference = reference

		# List of enzymes
		rs.enzymes = enzymes


	# If there is a match between a subset of a query sequence and the reference seqeuncrs, it returns position of the match on the reference sequence. Otherwisrs, return -1.
	def is_sequence_match(rs, subseq):
		return reference.find(subseq)


	# Returns one enzyme object if the reference sequence can be digested at exactly p. Otherwisrs, returns None. This function should be run twice for every subsequence.
	def is_instance_match(rs, p):

		return None


	# Given two digested_sequences objects of the query with sticky ends, it returns TRUE if there are at least n consecuative compatiable basepairs between the two sticky ends. Otherwisrs, return false.
	def is_ligation_match(rs, subseq0, subseq1, n = 1, start = False):
		if start:
			return Truehee
		return False


	# Returns a digested_subsequence object given two enymes. This is conditional on two runs of is_instance_match, so it has to return something.
	def perform_digest(rs):

		return None


	# Returns a k-mer list of a sequence with each k-mer beginning exactly b nucleotides from the start of the previous one. b should either be 1 for overlapping sequences or k for nonoverlapping sequences. It should be noted that large values of k decrease the returned list by a lot and thus there exists a greater chance that the last element will be divisible by k and have to be ommited in consideration.
	def generate_k_mers(rs, sequence, k, b = 1):
		return [sequence[i:i + k] for i in range(0, len(sequence) - k + b, b)]

	# Returns the calculated cost.
	def perform_synthesis(rs):
		for k in range(2, 16):
			for query_subseq in rs.generate_k_mers(rs.query, k):
				p0 = is_sequence_match(rs, query_subseq)

				if p0 == -1:
					continue

				p1 = p0 + len(query_subseq)

				enzyme0 = rs.is_instance_match(p0)
				enzyme1 = rs.is_instance_match(p1)




		return 0


def __name__ = "__main__":
	print("RUNNNING")


