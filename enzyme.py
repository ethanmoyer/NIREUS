
# name - Name of enzyme
# bp - Logical value for whether this object is storing a RestrictionType enzyme 
# cut0/cut1 - Store the points at which the enzyme cuts the restriction sequence
# restriction_site - Restriction site of enzyme
# cut_site0/cut_site1 - Store the left and right cut sites of the enzyme
# overlapping_seq - Store body of enzyme
# price - Price listed online
# units - units listed online

from Bio.Restriction import RestrictionBatch

class biopy_enzyme(enzyme):
	def __init__(be, bp_enzyme, price = None, units = None):
		enzyme.__init__(be, price, units)
		be.restriction_site = bp_enzyme.site
		be.name = RestrictionBatch([bp_enzyme]).as_string()[0]
		be.cut0 = bp_enzyme.charac[0]
		be.cut1 = bp_enzyme.charac[1] + be.site_len
		be.bp = True
		be.bp_enzyme = bp_enzyme


class csv_enzyme(enzyme):
	def __init__(ce, name, restriction_site, cut0, cut1, price = None, units = None):
		enzyme.__init__(ce, name, restriction_site, price, units)
		ce.restriction_site = restriction_site
		ce.name = name
		ce.cut0 = int(cut0)
		ce.cut1 = int(cut1)
		ce.bp = False

		ce.cut_site0 = ce.get_cut_sequence(0)
		ce.cut_site1 = ce.get_cut_sequence(1)

		ce.overlapping_seq = ce.get_overlapping_seq()


	# Returns the specified left (side == 0) or right (side == 1) cut site of the enzyme's restriction site with respect to the parent sequence.	
	def get_cut_sequence(ce, side):
		if side == 0:
			print()
			return e.restriction_site[:e.cut0 if e.cut0 <= e.cut1 else e.cut1]
		elif side == 1:
			if e.cut0 >= e.cut1:
				return e.restriction_site[e.cut0:len(e.restriction_site)]
			else:
				return e.restriction_site[e.cut1:len(e.restriction_site)]


	# Returns the overlapping body of the sequence. The body is defined as the stretch of the restriction sequence that creates those sticky ends. i.e. AATT in EcoRI's GAATTC 1 5 site.
	def get_overlapping_seq(ce):
		if e.cut0 <= e.cut1:
			return e.restriction_site[e.cut0:e.cut1]
		else:
			return e.restriction_site[e.cut1:e.cut0]


class enzyme:
	def __init__(e, restriction_site, price, units):
		e.restriction_site = restriction_site
		e.site_len = site_len
		e.price = price
		e.units = units
		e.bp = bp

		e.trans = {'A': 'A', 'T': 'T', 'C': 'C', 'G': 'G', 'R': ['A', 'G'], 
				   'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'], 
		    	   'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 
				   'D': ['A', 'G', 'T'], 'H': ['A', 'C', 'T'], 
				   'V': ['A', 'C', 'G'], 'N': ['A', 'T', 'C', 'G']}


	# Returns the price per unit of the enzyme
	def get_price_per_unit(e):
		if e.price != None and e.units != None:
			return e.price / e.units

	# Returns whether the given sequence is equivalent to either the left (side == 0) or right (side == 1) cut site
	def equal_to_cut_site(e, seq, side):
		cut = ''

		# Store specified cut site
		if side == 0:
			cut = e.cut_site0
		elif side == 1:
			cut = e.cut_site1


		# Loop through all the characters in the given sequence and cut sequence and return False if there is a single discrepency based on the trans dictionary
		for i in range(len(seq)):

			if seq[i] not in e.trans[cut[i]]:
				return False

		# Otherwise return true
		return True
