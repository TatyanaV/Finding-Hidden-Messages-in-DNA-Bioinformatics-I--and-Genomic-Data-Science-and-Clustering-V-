#!/usr/bin/env python
'''
CODE CHALLENGE: Implement MotifEnumeration (reproduced below).
     Input: Integers k and d, followed by a collection of strings Dna.
     Output: All (k, d)-motifs in Dna.

    MotifEnumeration(Dna, k, d)
        Patterns ? an empty set
        for each k-mer Pattern in Dna
            for each k-mer Pattern’ differing from Pattern by at most d
              mismatches
                if Pattern' appears in each string from Dna with at most d
                mismatches
                    add Pattern' to Patterns
        remove duplicates from Patterns
        return Patterns


Return to main text

Sample Input:

3 1
ATTTGGC
TGCCTTA
CGGTATC
GAAAATT

Sample Output:

ATA ATT GTT TTT
'''
def MismatchList(kmer, d):
	'''Returns a list of all k-mers that mismatch a given k-mer by at most d characters.'''
	kmer_mismatches = [kmer]
	for i in range(1,d+1):
		# Each combination gives the indicies we want to mismatch.
		kmer_mismatches += CreateMismatches([[kmer, list(combo)] for combo in combinations(range(len(kmer)),i)])
	return kmer_mismatches

def CreateMismatches(swap_list):
	'''Generates k-mer mismatches by replacing the characters at given indicies with mismatching characters.'''
	nucleotides = 'ACGT'
	mismatch_list = []
	# Swap the i-th character of string with the character ch.
	swap = lambda string, ch, i: string[:index]+ch+string[index+1:]

	# If we have more than one index left to mismatch, repeat the process.
	if len(swap_list[0][1]) > 1:
		for kmer, indicies in swap_list:
			index = indicies[0]
			for nuc in filter(lambda n: n != kmer[index], nucleotides):
				mismatch_list.append([swap(kmer, nuc, index), indicies[1:]])

		return CreateMismatches(mismatch_list)

	# Otherwise, on the final mismatch return the list of k-mers.
	else:
		for kmer, [index] in swap_list:
			for nuc in filter(lambda n: n != kmer[index], nucleotides):
				mismatch_list.append(swap(kmer, nuc, index))

		return mismatch_list

with open('C:/Users/Tanya/PycharmProjects/untitled/file1.txt') as input_data:
	k, d = map(int, input_data.readline().split())
	dna_list = [line.strip() for line in input_data.readlines()]

def MotifEnumeration(k, d, dna_list):
    # Generate sets of (k,d)-motifs for each dna sequence in the list.
    motif_sets = [{kmer for i in range(len(dna)-k+1) for kmer in MismatchList(dna[i:i+k], d)} for dna in dna_list]

    # Intersect all sets to get the common elements.  The answers are displayed as sorted, so we'll sort too.
    motifs = sorted(list(reduce(lambda a,b: a&b, motif_sets)))


print (*MotifEnumeration(k, d, dna_list),sep=' ')