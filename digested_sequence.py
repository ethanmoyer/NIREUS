class digested_sequuence:
	def __init__(ds, sticky0, sticky1, sequence):

		# Both sticky ends (left and right, respectively) encoded based on sticky stranded alphabet with respect to top sequence.
		ds.sticky0 = sticky0
		ds.sticky1 = sticky1

		# Top sequence between two sticky ends
		ds.sequence = sequence

