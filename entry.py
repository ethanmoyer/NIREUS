class entry:
	def __init__(e, sequence, length, p_a, p_t, p_c, p_c, p_g, r1, r2, c, gp_a, gp_t, gp_c, gp_g, gr1, gr2, en, output):

		# 'SEQ' feature
		e.sequence = sequence

		# Length of 'SEQ'
		e.length = length

		# Nucleotide proportions of the entry
		e.p_a = p_a
		e.p_t = p_t
		e.p_c = p_c
		e.p_g = p_g

		# Complexity ratings using the complete and proportional method of the entry
		e.r1 = r1
		e.r2 = r2

		# Number of times that the sequence appears in a given reference
		e.c = c

		# Nucleotide proportions of the query
		e.gp_a = gp_a
		e.gp_t = gp_t
		e.gp_c = gp_c
		e.gp_g = gp_g

		# Complexity ratings using the complete and proportional method of the query
		e.gr1 = gr1
		e.gr2 = gr2

	def get_separated_seq(e):
		return [x for x in e.sequence]

