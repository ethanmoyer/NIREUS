class enzyme:
	def __init__(e, name, restriction_site, cuts, price, units):
		# Name of enzyme
		e.name = name

		# Restriction site of enzyme
		e.restriction_site = restriction_site

		# At least a pair of cuts the enzyme performs in a digest (top and bottom)
		e.cuts = cuts

		# Price and units were obtained using a restriction enzyme database from NEB
		# Price listed online
		e.price = price

		# Units listed online
		e.units = units

	# Returns the price per unit of the enzyme
	def get_price_per_unit(e):
		return e.price / e.units

	
		

