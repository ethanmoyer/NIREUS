# pkill -f python
import csv

from enzyme import enzyme, csv_enzyme, biopy_enzyme
from entry import entry
from digested_sequence import digested_sequence
from restriction_map import restriction_map

import pandas as pd
import re

from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

# Simple function to clear the console so the synthesis can be visualized
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

	def __init__(rs, query, reference, enzymes, bp = False):
		# Query sequence
		rs.query = query

		# Reference sequence
		rs.reference = reference

		# List of enzyme objects
		rs.enzymes = enzymes

		# Array to keep track of synthesis
		rs.synthesized_query = []
		rs.enzyme_list = []
		rs.positions = []

		# This alphabet is used to translate between sticky on bottom strand and sticky end on top strand.
		rs.ligation_alphabet = {'A': 'Q', 'T': 'W', 'C': 'E', 'G': 'R'}


	# If there is a match between a subset of a query sequence and the reference sequence it returns position of the match on the reference sequence. Otherwise, return -1.
	def find_all_sequence_match(rs, subseq):
		return [x.start() for x in re.finditer(subseq, rs.reference)]


	# Returns one enzyme object if the reference sequence can be digested at exactly p. Otherwise, returns None. This function should be run twice for every subsequence. This is specifically designed for the BioPython package. 
	# Return list of enzymes
	def is_instance_match(rs, p, ignore_enzymes = []):
		enzyme_list = []
		# Loop through each individual enzyme
		for enzyme in rs.enzymes:

			# Loop through each individual enzyme
			if enzyme.name in ignore_enzymes:
				continue

			# Retrieve both cutting sites of the enzymes
			site_seq0 = enzyme.cut_site0
			site_seq1 = enzyme.cut_site1

			# Return an enzyme it can cut the reference sequence at exactly p. This means that left restriction site, the right restriction site, and the body of the restriction site all match their relative segments in the reference relative to p.
			if enzyme.equal_to_cut_site(rs.reference[p - len(site_seq0):p], 0) and rs.reference[p:p + len(enzyme.overlapping_seq)] == enzyme.overlapping_seq and enzyme.equal_to_cut_site(rs.reference[p + len(enzyme.overlapping_seq): p + len(enzyme.overlapping_seq) + len(site_seq1)], 1):

				enzyme_list.append(enzyme)

		return enzyme_list

	# Returns one enzyme object if the reference sequence can be digested at exactly p. Otherwise, returns None. This function should be run twice for every subsequence. This is specifically designed for the BioPython package. Return enzyme that is compatible with the most amount of enzymes.
	def is_instance_match_biopython(rs, p, ignore_enzymes = []):
		enzyme_list = []
		# Loop through each individual enzyme
		for enzyme in rs.enzymes:

			if enzyme.name in ignore_enzymes:
				continue

			# Retrieve both cutting site lengths of the enzymes
			site0 = enzyme.cut0
			site1 = enzyme.cut1

			seq = Seq(rs.reference[p - site0:p + enzyme.site_len - site0], amb)
			digestion = enzyme.bp_enzyme.catalyze(seq)

			if len(digestion) != 2:
				continue

			if len(digestion[0]) == len(enzyme.restriction_site[:site0]) and len(digestion[1]) == len(enzyme.restriction_site[site0:]):
				enzyme_list.append(enzyme)

		return enzyme_list


	# Given two digested_sequences objects of the query with sticky ends, it returns TRUE if there are at least n consecutive compatible base pairs between the two sticky ends. Otherwise, return false. If subseq0 isn't provided, this means that this is the first round of ligations, which automatically returns True.
	def is_ligation_match(rs, subseq1, subseq0 = None):
		# If there is not a subseq0 provided, then this is the beginning of synthesis, thus return True. Also, if both ligation ends are blunt ends, return True.
		if subseq0 == None or (subseq0 == "" and subseq1 == ""):
			return True

		# Return false if either one of them are blunt (conditioned that they both aren't) and if the new sticky overhand is shorter than the minimum value, n.
		if (subseq0 == "" or subseq1 == "") or (len(subseq0) < n):
			return False
		
		# Loop through all values of n, and check whether there are any mismatches between the two sticky ends, subseq0 and subseq1.
		for i in range(len(subseq1)):
			if not rs.is_sticky_match(subseq0[len(subseq0) - i - 1], subseq1[i]):
				return False

		# If there are not any mismatches, return true.
		return True


	# Given two enzymes, this function returns whether the two sites are compatible. This might result in a more conservative ligation match between the two enzymes. Where enzyme0 is last enzyme1 and enzyme1 is current enzyme0
	def is_ligation_match_biopython(rs, enzyme1, enzyme0):
		return enzyme0.bp_enzyme % enzyme1.bp_enzyme


	# Returns a digested_subsequence object given two enzymes and two positions on a reference sequence. This is conditional on two runs of is_instance_match, so it has to return something.
	def perform_digest(rs, enzyme0, enzyme1, p0, p1):

		# The right most sites of the digested sequence's two sticky ends
		z0 = p0 + len(enzyme0.overlapping_seq)
		z1 = p1 + len(enzyme1.overlapping_seq)

		# Return a digest sequence object
		return digested_sequence(rs.reference[p0:z0], rs.translate_stikcy_end(rs.reference[p1:z1]), rs.reference[z0:p1])


	# Returns a sticky end string based on the rs.ligation_alphabet dictionary.
	def translate_stikcy_end(rs, seq):

		seq_ = []
		# Loop through all of the nucleotides in the sequence and convert it based on the ligation alphabet described earlier in the code. Return this converted sequence.
		for i in range(len(seq)):
			seq_.append(rs.ligation_alphabet[seq[i]])
		return ''.join(seq_)


	# Returns whether the top nucleotide, s1, can ligate to the bottom nucleotide, s0.
	def is_sticky_match(rs, s0, s1):
		return rs.ligation_alphabet[s1] == s0


	def most_compatible_enzyme(rs, enzyme_list):
		enzyme_ = None
		max_count = 0
		for enzyme_i in enzyme_list:
			count = 0
			for enzyme_j in enzymes:
				if rs.is_ligation_match_biopython(enzyme_i, enzyme_j) and enzyme_i != enzyme_j:
					count += 1
			if count > max_count:
				enzyme_ = enzyme_i

		return enzyme_

	def perform_synthesis_bp(rs):
		# Set first enzyme1 to None to imply a primer connection in to the gene... This will be corrected later.
		last_enzyme1 = None

		# Variable to dictate whether synthesis is running
		synthesis = True

		# Length of synthesized query
		synthesized_query_len = 0

		# Continue while synthesis is true.
		while (synthesis):

			# Do we ever come back to this spot?
			k = 16
			# Loop through each of the parsing values--this range needs to be explored later.
			while k >= 4:

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

					# Return list of enzymes instead of only one... 
					# Find two enzymes that cut at exactly p0 and p1
					enzymes0 = rs.is_instance_match_biopython(p0, ignore_enzymes = ignore_enzymes)
					enzymes1 = rs.is_instance_match_biopython(p1)

					# If no two enzymes exist for a given query_subseq, continue to the next p0
					if enzymes0 == [] or enzymes1 == []:
						#print('Continue to next pair if present.')
						continue

					enzyme0 = rs.most_compatible_enzyme(enzymes0)
					enzyme1 = rs.most_compatible_enzyme(enzymes1)

					# If the previous sticky end and the current sticky end are not compatible, continue to the next p0.
					if last_enzyme1 is not None and not rs.is_ligation_match_biopython(last_enzyme1, enzyme0):
						ignore_enzymes.append(enzyme0.name)
						continue

					print('--------------------------------------')
					rs.display_enzyme(enzyme0)
					rs.display_enzyme(enzyme1)

					# Save the current sticky end
					last_enzyme1 = enzyme1 

					# Used to track the synthesis step by step in the display_synthesis function. 
					rs.enzyme_list.append([enzyme0, enzyme1])
					rs.synthesized_query.append(query_subseq)
					rs.positions.append([p0, p0 + len(query_subseq)])

					# Store length of rs.synthesized_query
					synthesized_query_len = len(''.join(rs.synthesized_query))

					ignore_enzymes = []

					# Displays the position that query was found in the reference and the query sequence itself
					print('Percentage: {percent} --> {p0} \t {query_subseq}'.format(percent = synthesized_query_len / len(query) * 100, p0 = p0, query_subseq = query_subseq))
					print()

					# Synthesis is completed if the synthesized query is at least the size of the query
					if synthesized_query_len >= len(rs.query):
						print('Synthesis completed')
						return 0

					# Since subsequence was found, reset parsing value to the max
					k = 16

					# Break out of loop back to beginning of while
					break
				else:
					# If query_subseq is not in the reference, subtract k by 1
					k -= 1

			# Reset the current sticky end
			last_enzyme1 = None

			# Used to track the synthesis step by step in the display_synthesis function. 
			rs.synthesized_query.append('X')
			rs.enzyme_list.append([enzyme(), enzyme()])
			rs.positions.append([None, None])

			# Store length of rs.synthesized_query
			synthesized_query_len = len(''.join(rs.synthesized_query))

	# Returns the calculated cost.
	def perform_synthesis(rs):

		# Set the first sticky end of synthesis to None to imply a primer connection at least 10 bp in to the gene... This will be corrected later
		last_sticky_end = None

		# Variable to dictate whether synthesis is running
		synthesis = True

		# Length of synthesized query
		synthesized_query_len = 0

		# Continue while synthesis is true.
		while (synthesis):

			# Do we ever come back to this spot?
			k = 16
			# Loop through each of the parsing values--this range needs to be explored later.
			while k >= 4:

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

					# Return list of enzymes instead of only one... 
					# Find two enzymes that cut at exactly p0 and p1
					enzyme0 = rs.is_instance_match(p0, ignore_enzymes = ignore_enzymes)
					enzyme1 = rs.is_instance_match(p1)

					# If no two enzymes exist for a given query_subseq, continue to the next p0
					if enzyme0 == [] or enzyme1 == []:
						#print('Continue to next pair if present.')
						continue

					# Make a digested sequence using the two enzymes at exactly p0 and p1
					digested_seq = rs.perform_digest(enzyme0, enzyme1, p0, p1)

					# If the previous sticky end and the current sticky end are not compatible, continue to the next p0
					if not rs.is_ligation_match(digested_seq.sticky0, last_sticky_end):
						ignore_enzymes.append(enzyme0.name)
						continue

					print('--------------------------------------')
					rs.display_enzyme(enzyme0)
					rs.display_enzyme(enzyme1)

					# Display the results of the flanked digest
					rs.display_digested_seq(digested_seq)
					print('--------------------------------------')

					# Save the current stick end
					last_sticky_end = digested_seq.sticky1

					# Used to track the synthesis step by step in the display_synthesis function. 
					rs.synthesized_query.append(query_subseq)
					rs.enzyme_list.append([enzyme0, enzyme1])
					rs.positions.append([p0, p0 + len(query_subseq)])

					# Store length of rs.synthesized_query
					synthesized_query_len = len(''.join(rs.synthesized_query))

					ignore_enzymes = []

					# Displays the position that query was found in the reference and the query sequence itself
					print('Percentage: {percent} --> {p0} \t {query_subseq}'.format(percent = synthesized_query_len / len(query) * 100, p0 = p0, query_subseq = query_subseq))
					print()

					# Synthesis is completed if the synthesized query is at least the size of the query
					if synthesized_query_len >= len(rs.query):
						print('Synthesis completed')
						return 0

					# Since subsequence was found, reset parsing value to the max
					k = 16

					# Break out of loop back to beginning of while
					break
				else:
					# If query_subseq is not in the reference, subtract k by 1
					k -= 1

			# Reset the current sticky end
			last_sticky_end = None

			# Used to track the synthesis step by step in the display_synthesis function. 
			rs.synthesized_query.append('X')
			rs.enzyme_list.append([enzyme(), enzyme()])
			rs.positions.append([None, None])

			# Store length of rs.synthesized_query
			synthesized_query_len = len(''.join(rs.synthesized_query))


	# This function displays relevant information about the given enzyme.
	def display_enzyme(rs, enzyme_):
		print('----', enzyme_.name, '----', sep='')
		print('Restriction site:\t', enzyme_.restriction_site)
		print('Top cut position:\t', enzyme_.cut0)
		print('Bottom cut position:\t',enzyme_.cut1)
		if not enzyme_.bp:
			print('Body restriction site:\t', enzyme_.overlapping_seq)
			print('Left sequence:\t\t', enzyme_.cut_site0)
			print('Right sequence:\t\t', enzyme_.cut_site1)
		print()


	def display_digested_seq(rs, digested_sequence_):
		print('----', digested_sequence_.sticky0, digested_sequence_.sequence, digested_sequence_.sticky1, '----', sep='')
		print('Left sticky overhang\t:', digested_sequence_.sticky0)
		print('Body sequence:\t', digested_sequence_.sequence)
		print('Right sticky overhang:\t', digested_sequence_.sticky1)
		print()


	# Displays to console the results of the synthesis.
	def display_synthesis(rs):
		temp_synthesized_query_display = ''
		temp_enzyme_display0 = ''
		temp_enzyme_display1 = ''

		# Print each of the ith groups of data
		for i in range(len(rs.synthesized_query)):

			# Append data to temporary variable with a tab as spacing
			temp_synthesized_query_display += rs.synthesized_query[i] + '\t'
			temp_enzyme_display0 += rs.enzyme_list[i][0].name + '\t'
			temp_enzyme_display1 += rs.enzyme_list[i][1].name + '\t'
			if len(rs.synthesized_query[i]) < 8:
				temp_synthesized_query_display += '\t'
			if len(rs.enzyme_list[i][0].name) < 8:
				temp_enzyme_display0 += '\t'
			if len(rs.enzyme_list[i][1].name) < 8:
				temp_enzyme_display1 += '\t'

			# Add string of 6 objects one at a time to the three arrays
			if (i + 1) % 6 == 0 or i == len(rs.synthesized_query) - 1:
				print(temp_synthesized_query_display)
				print(temp_enzyme_display0)
				print(temp_enzyme_display1)
				print()
				temp_synthesized_query_display = ''
				temp_enzyme_display0 = ''
				temp_enzyme_display1 = ''

		accuracy = (1 - rs.synthesized_query.count('X') / len(rs.query)) * 100
		print('Synthesis accuracy: ', str(round(accuracy, 2)))

		rm = restriction_map(len(reference))

		for i in range(len(rs.enzyme_list)):
			if rs.positions[i][0] is not None or rs.positions[i][1] is not None:
				rm.add_enzyme(rs.enzyme_list[i][0].name + ' ' + rs.enzyme_list[i][1].name, rs.positions[i][0])

		rm.show_map()

if __name__ == '__main__':

	# Add conformation for all three different steps

	# Load the reference sequence
	print('Loading reference sequence')
	reference_file_location = 'data/references/reference0.txt'
	print(f'Default reference file location: {reference_file_location}')
	ref_ans = str(input('Would you like to use the default reference file? [yes/no]: '))
	if False and ref_ans.lower() == 'no':
		reference_file_location = str(input('Please indicate which reference file you would like to use: '))

	reference = open(reference_file_location).read()
	reference = re.sub('\n', '', reference)
	print('Reference length: ', len(reference), sep = '')
	print()
	# Load the query sequence
	print('Loading query sequence')
	query_file_location = 'data/queries/query0.txt'
	print(f'Default query file location: {query_file_location}')
	query_ans = str(input('Would you like to use the default query file? [yes/no]: '))
	if False and query_ans.lower() == 'no' :
		query_file_location = str(input('Please indicate which query file you would like to use: '))

	query = open(query_file_location).read()
	query = re.sub('\n', '', query)
	query = query
	print('Query length: ', len(query), sep = '')

	print(f'\nReference length to query length ratio:',str(round(len(reference) / len(query), 2)))
	print()
	# Load the enzymes
	print('Loading enzyme list')
	data_ans = str(input('Would you like to use an enzyme file or a database from BioPython? [file/database]: '))
	if False and data_ans.lower() == 'file':	
		enzymes_file_location = 'data/enzymes/enzymes0.csv'
		print(f'Default enzyme file location: {enzymes_file_location}')
		enz_ans = str(input('Would you like to use the default enzyme file? [yes/no]: '))
		if enz_ans.lower() == 'no' :
			enzymes_file_location = str(input('Please indicate which enzyme file you would like to use: '))

		with open(enzymes_file_location, newline='') as csvfile:
			reader = csv.DictReader(csvfile)
			enzymes = [csv_enzyme(row['name'], translate_site(row['site'], int(row['cut0']), int(row['cut1'])), row['cut0'], row['cut1']) for row in reader]
		# Start the synthesis 
		print('Starting synthesis')

		# Initialize the restriction synthesis object
		rs = restriction_synthesis(query, reference, enzymes)

		# Perform synthesis
		rs.perform_synthesis()

	elif True or data_ans.lower() == 'database':
		amb = IUPACAmbiguousDNA()
		RestrictionBatch.show_codes()
		print('ALL = All of them')
		comp_ans = str(input('Please select which company/companies that will be used for synthesis. i.e. B,C: '))
		if comp_ans is 'N':
			enzymes = CommOnly
		else:
			companies = comp_ans.split(',')
			enzymes = [biopy_enzyme(e) for e in RestrictionBatch(first=[], suppliers=companies)]

		# Start the synthesis 
		print('Starting synthesis')

		# Initialize the restriction synthesis object
		rs = restriction_synthesis(query, reference, enzymes)

		# Perform synthesis
		rs.perform_synthesis_bp()


	# Display synthesis
	rs.display_synthesis()

