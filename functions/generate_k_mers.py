# Returns a k-mer list of a sequence with each k-mer beginning exactly b nucleotides from the start of the previous one. b should either be 1 for overlapping sequences or k for nonoverlapping sequences. It should be noted that large values of k decrease the returned list by a lot and thus there exists a greater chance that the last element will be divisible by k and have to be ommited in consideration.
def generate_k_mers(rs, sequence, k, b = 1):
	return [sequence[i:i + k] for i in range(0, len(sequence) - k + b, b)]

