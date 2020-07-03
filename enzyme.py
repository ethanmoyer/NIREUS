class enzyme:
	def __init__(e, name, restriction_site, cut0, cut1, price = None, units = None):
		# Name of enzyme
		e.name = name

		# Restriction site of enzyme
		e.restriction_site = restriction_site

		# The points at which the enzyme cuts the restriction sequence
		e.cut0 = int(cut0)
		e.cut1 = int(cut1)

		# The left and right cut sites of the enzyme
		e.cut_site0 = e.get_cut_sequence(0)
		e.cut_site1 = e.get_cut_sequence(1)

		# Price and units were obtained using a restriction enzyme database from NEB
		# Price listed online
		e.price = price

		# Units listed online
		e.units = units

		e.restriction_alphabet = {'A': 'A', 'T': 'T', 'C': 'C', 'G': 'G', 'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'], 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'T', 'C', 'G']}


	# Returns the price per unit of the enzyme
	def get_price_per_unit(e):
		if e.price != None and e.units != None:
			return e.price / e.units

	
	# Returns the specified left (side == 0) or right (side == 1) cut site of the enzyme's restriction site with respect to the parent sequence.	
	def get_cut_sequence(e, side):
		if side == 0:
			return e.restriction_site[:e.cut0]
		elif side == 1:
			return e.restriction_site[e.cut0:e.cut1 + 1]


	# Returns whether the given sequence is equivalent to either the left (side == 0) or right (side == 1) cut site
	def equal_to_cut_site(e, seq, side):
		cut = ''

		# Store specified cut site
		if side == 0:
			cut = e.cut_site0
		elif side == 1:
			cut = e.cut_site1


		# Loop through all the characters in the given sequence and cut sequence and return False if there is a single discrepency based on the restriction_alphabet dictionary
		for i in range(len(seq)):
			if seq[i] not in e.restriction_alphabet[cut[i]]:
				return False

		# Otherwise return true
		return True
