import csv

from enzyme import enzyme
from entry import entry
from digested_sequence import digested_sequence

import pandas as pd
import re

# Simple funciton to clear the console so the synthesis can be visualized
def clear(): os.system('clear')

# Returns a translated site given certain parameters: 'N#' is translated to 'N' * # and the site can be subscripted with both cuts.
def translate_site(site, cut0, cut1):

	# Find the number(s) following an N
	number_of_Ns = re.findall('N(\d+)', site)

	# If it doesn't exist, simply return the site.
	if number_of_Ns == []:

		# Extend the site toward the right side if need be to fit in both of the restriction sites
		if cut0 > len(site) or cut1 > len(site):
			site += 'N' * (cut0 - len(site)) if cut0 > cut1 else 'N' * (cut1 - len(site))

		return site
	
	# Substitute the 'N#' with 'N' * #
	site = re.sub('N\d+', 'N' * int(number_of_Ns[0]), site)

	# Extend the site toward the right side if need be to fit in both of the restriction sites
	if cut0 > len(site) or cut1 > len(site):
		site += 'N' * (cut0 - len(site)) if cut0 > cut1 else 'N' * (cut1 - len(site))

	return site


class restriction_synthesis():

	def __init__(rs, query, reference, enzymes):
		# Query sequence
		rs.query = query

		# Reference sequence
		rs.reference = reference

		# List of enzyme objects
		rs.enzymes = enzymes

		# Array to keep track of synthesis
		rs.synthesized_query = []
		rs.enzyme_list = []

		# This alphabet is used to traslate between sticky on bottom strand and sticky end on top strand.
		rs.ligation_alphabet = {'A': 'Q', 'T': 'W', 'C': 'E', 'G': 'R'}


	# If there is a match between a subset of a query sequence and the reference seqeuncrs, it returns position of the match on the reference sequence. Otherwisrs, return -1.
	def find_all_sequence_match(rs, subseq):
		return [x.start() for x in re.finditer(subseq, rs.reference)]


	# Returns one enzyme object if the reference sequence can be digested at exactly p. Otherwisrs, returns None. This function should be run twice for every subsequence.
	def is_instance_match(rs, p, ignore_enzymes = []):

		# Loop through each individual enzyme
		for enzyme in rs.enzymes:

			if enzyme.name in ignore_enzymes:
				continue

			# Retrieve both cutting sites of the enzymes
			site0 = enzyme.cut_site0
			site1 = enzyme.cut_site1

			# Return an enzyme it can cut the reference sequence at eactly p. This means that left restriction site, the right restriction site, and the body of the restriction site all match their relative segments in the reference relative to p.
			if enzyme.equal_to_cut_site(rs.reference[p - len(site0):p], 0) and rs.reference[p:p + len(enzyme.overlapping_seq)] == enzyme.overlapping_seq and enzyme.equal_to_cut_site(rs.reference[p + len(enzyme.overlapping_seq): p + len(enzyme.overlapping_seq) + len(site1)], 1):

				return enzyme

		# If no such enzymes exist, return None
		return None


	# Returns a digested_subsequence object given two enzymes and two positions on a reference sequence. This is conditional on two runs of is_instance_match, so it has to return something.
	def perform_digest(rs, enzyme0, enzyme1, p0, p1):

		# The right most sites of the digested sequence's two sticky ends
		z0 = p0 + len(enzyme0.overlapping_seq)
		z1 = p1 + len(enzyme1.overlapping_seq)

		# Return a digest sequence object
		return digested_sequence(rs.reference[p0:z0], rs.translate_stikcy_end(rs.reference[p1:z1]), rs.reference[z0:p1])


	# Returns a sticky end string based on the rs.ligation_alphabet dictionary.
	def translate_stikcy_end(rs, seq):

		# Loop through all of the nucleotides in the sequence and convert it based on the ligation alphabet described earlier in the code. Return this converted sequence.

		seq_ = []

		for i in range(len(seq)):
			seq_.append(rs.ligation_alphabet[seq[i]])
		return ''.join(seq_)


	# Given two digested_sequences objects of the query with sticky ends, it returns TRUE if there are at least n consecuative compatiable basepairs between the two sticky ends. Otherwisrs, return false. If subseq0 isn't provided, this means that this is the first round of ligations, which automatically returns True.
	def is_ligation_match(rs, subseq1, subseq0 = None, n = 1):
		# If there is not a subseq0 provided, then this is the beginning of synthesis, thus return True. Also, if both ligation ends are blunt ends, return True.
		if subseq0 == None or (subseq0 == "" and subseq1 == ""):
			return True

		# Return false if either one of them are blunt (conditioned that they both aren't) and if the new sticky overhand is shorter than the minimum value, n.
		if (subseq0 == "" or subseq1 == "") or (len(subseq0) < n):
			return False

		print(subseq0)
		print(subseq1)
		
		# Loop through all values of n, and check whether there are any mismatches between the two sticky ends, subseq0 and subseq1.
		for i in range(n):
			if not rs.is_sticky_match(subseq0[len(subseq0) - i - 1], subseq1[i]):
				return False

		# If there are not any mismatches, return true.
		return True


	# Returns whether the top nucletoide, s1, can ligate to the bottomn nucleotide, s0.
	def is_sticky_match(rs, s0, s1):
		return rs.ligation_alphabet[s1] == s0


	# Returns the calculated cost.
	def perform_synthesis(rs):

		# Set the first sticky end of synthesis to None to imply a primer connection at least 10 bp in to the gene... This will be corrected later
		last_sticky_end = None

		# Variable to dictate whether synthesis is running
		synthesis = True

		# Length of synthesized query
		synthesized_query_len = 0

		k_sug = None

		k = 16

		# Continue while synthesis is true.
		while (synthesis):
			# Loop through each of the parsing values--this range needs to be explored later.
			while k >= 2:

				ignore_enzymes = []

				# Loop through all of the fragmented overlapping k-mer sequences for each parsing value
				query_subseq = rs.query[synthesized_query_len:synthesized_query_len + k + 1]

				print('Parsing for:', query_subseq, sep='')

				# Find all position of query_subseq in the reference
				for p0 in rs.find_all_sequence_match(query_subseq):

					# If there isn't a match, continue to next p0
					if p0 == -1:
						continue

					# The end position of the match of p0
					p1 = p0 + len(query_subseq)

					# Find two enaymes that cut at exactly p0 and p1
					enzyme0 = rs.is_instance_match(p0, ignore_enzymes = ignore_enzymes)
					enzyme1 = rs.is_instance_match(p1)

					# If no two enaymes exist for a given query_subseq, continue to the next p0
					if enzyme0 == None or enzyme1 == None:
						continue

					#print(rs.reference[p0 - 5: p0 + len(query_subseq) + 5])

					#print(query_subseq)
					#print(p0)
					print('--------------------------------------')
					#print("Enzyme0")
					rs.display_enzyme(enzyme0)
					#print("Enzyme1")
					rs.display_enzyme(enzyme1)

					# Make a digested sequence using the two enzymes at exactly p0 and p1
					digested_seq = rs.perform_digest(enzyme0, enzyme1, p0, p1)

					# Display the results of the flanked digest
					rs.display_digested_seq(digested_seq)
					print('--------------------------------------')
					# If the previous sticky end and the current sticky end are not compatiable, continue to the next p0
					if not rs.is_ligation_match(digested_seq.sticky0, last_sticky_end):
						ignore_enzymes.append(enzyme0.name)
						continue

					# Make the current sticky end the current sticky end
					last_sticky_end = digested_seq.sticky1

					# Used to track the synthesis step by step in the display_synthesis funciton. 
					rs.synthesized_query.append(query_subseq)
					rs.enzyme_list.append([enzyme0, enzyme1])

					# Store length of rs.synthesized_query
					synthesized_query_len = len(''.join(rs.synthesized_query))

					ignore_enzymes = []

					# Displays the positon that query was found in the reference and the query sequence itself
					print('Percentage: {percent} --> {p0} \t {query_subseq}'.format(percent = synthesized_query_len / len(query), p0 = p0, query_subseq = query_subseq))
					print()

					# Synthesis is completed if the synthesized query is at least the size of the query
					if synthesized_query_len >= len(rs.query):
						print('Synthesis completed')
						return 0

					k = 16

					break
				else:
					k -= 1

			rs.synthesized_query.append('X')
			rs.enzyme_list.append([None, None])

		return 0


	# This function displays relevant information about the given enzyme.
	def display_enzyme(rs, enzyme_):
		print('----', enzyme_.name, '----', sep='')
		print('Restriciton site:\t', enzyme_.restriction_site)
		print('Top cut position:\t', enzyme_.cut0)
		print('Bottom cut position:\t',enzyme_.cut1)
		print('Body restriction site:\t', enzyme_.overlapping_seq)
		print('Left sequence:\t\t', enzyme_.cut_site0)
		print('Right sequence:\t\t', enzyme_.cut_site1)
		print()


	def display_digested_seq(rs, digested_sequence_):
		print('----', digested_sequence_.sticky0, digested_sequence_.sequence, digested_sequence_.sticky1, '----', sep='')
		print('Left stikcy overhang\t:', digested_sequence_.sticky0)
		print('Body sequence:\t', digested_sequence_.sequence)
		print('Right stikcy overhang:\t', digested_sequence_.sticky1)
		print()


	# Displays to console the results of the synthesis.
	def display_synthesis(rs):

		# Three arrays to display n lines of 10 subsequences created by flanking with two enzymes
		synthesized_query_display = []
		enzyme_display0 = []
		enzyme_display1 = []

		# Three temporary strings for each ith line of subsequences
		temp_synthesized_query_display = ''
		temp_enzyme_display0 = ''
		temp_enzyme_display1 = ''

		# Loop through all of the objects in in each iteration of digest
		for i in range(len(rs.synthesized_query)):

			# Append data to temporary variable with a tab as spacing
			temp_synthesized_query_display += rs.synthesized_query[i] + '\t'
			temp_enzyme_display0 = rs.enzyme_list[i][0] + '\t'
			temp_enzyme_display1 = rs.enzyme_list[i][1] + '\t'

			# Add string of 10 objects one at a time to the three arrays
			if (i + 1) % 10 == 0:
				temp_synthesized_query_display += '\n'
				temp_enzyme_display0 += '\n'
				temp_enzyme_display1 += '\n'

				synthesized_query_display.append(temp_synthesized_query_display)
				enzyme_display0.append(temp_enzyme_display0)
				enzyme_display1.append(temp_enzyme_display1)

				temp_synthesized_query_display = ''
				temp_enzyme_display0 = ''
				temp_enzyme_display1 = ''
			
		# Print each of the ith groups of data
		for i in range(len(synthesized_query_display)):
			print(synthesized_query_display[i])
			print(temp_enzyme_display0[i])
			print(temp_enzyme_display1[i])
			print()
			

if __name__ == '__main__':
	# Load the reference sequence
	print('Loading reference')
	reference_file_location = 'data/references/reference0.txt'
	reference = open(reference_file_location).read()
	re.sub('\n', '', reference)

	# Load the query sqeuence
	print('Loading query')
	query_file_location = 'data/queries/query0.txt'
	query = open(query_file_location).read()
	re.sub('\n', '', query)

	# Load the enzymes
	print('Loading enzymes')
	enzymes_file_location = 'data/enzymes/enzymes0.csv'
	with open(enzymes_file_location, newline='') as csvfile:
		reader = csv.DictReader(csvfile)
		enzymes = [enzyme(row['name'], translate_site(row['site'], int(row['cut0']), int(row['cut1'])), row['cut0'], row['cut1']) for row in reader]

	# Start the synthesis 
	print('Starting synthesis')

	# Initialize the restriction synthesis object
	rs = restriction_synthesis(query, reference, enzymes)

	# Perform the synthesis
	rs.perform_synthesis()
