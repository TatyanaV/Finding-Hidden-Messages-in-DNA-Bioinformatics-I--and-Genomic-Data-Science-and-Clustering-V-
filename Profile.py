'''
Profile-most Probable k-mer Problem: Find a Profile-most probable k-mer in a string.
     Input: A string Text, an integer k, and a 4 × k matrix Profile.
     Output: A Profile-most probable k-mer in Text.

CODE CHALLENGE: Solve the Profile-most Probable k-mer Problem.

Return to main text

Sample Input:

ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT
5
0.2 0.2 0.3 0.2 0.3
0.4 0.3 0.1 0.5 0.1
0.3 0.3 0.5 0.2 0.4
0.1 0.2 0.1 0.1 0.2

Sample Output:

CCGAG
'''
from functools import reduce
text = 'GACCTGACTACTTAACGTCTACCCGCGCTTATTCAGCGTGGAGTGATTTCTGTTGAACCTTAAAGTGGTTGTGTCTGGATAGTGTAAGGGCGACGGACCCCAGTATTCACGGACTACCCAAGCGGATCTCCGATATCCAGTCCCGATGAAGGGTATAGTCTACCCTCGCGCTTCAAAACTGGGTGAGTTGGAATTTCAAGCAAATCCGGCGGTTGATTGAAGCCACCCTATCTATTTCTCCGCTCTTAACTTGCCGCTGAGCGGTAACCCCATCAGCAAACAGACTTCTGTATAGGGGCGATGCGACTTATAATATCTCCAGCGACAAGAAGTCATCGCCTCACGTCTTCAGTACGGCCGAGATGAGTGTCAGCTTTGAGTTCTGGGAGCTCCACACACAATTCGATGGCTAGTTCTGAAGTCCGCGCCCCTAGCATCGCGATCGGGCAGGCTGCAGCACTCTGGAGCTCTATACGTACAGCATCGTAAGGAGGACTGATCACCTGCTCCTAAGGATGTTACGAAACTCCATGCTTCTTGTGGGTATGGAGCCCGTGATGGGGGCGCTCGCTTCAGGGTCAAAAAGCGCGGTACTAGCTCACGTTGGACGTACCGCGCGCTGACATGCCACGGGAAGCTAATCCAAATCGAGTCACGATATGGCCGATTCACGGTTTCCAATCATTGAACACGCCGGACGTGCTACAGAATGGCTAATGCTTTAACTGGCAATTTAAGGCAGGGCAGGGGATCAAATGGACCATGAGAGTGGTTCTAACTAGACGGACGGTCAAACGGCCGTCACTAGGCGGATCTCGTCACTCACATCCGGGAAATGACGGGTGTTGGTACGATTTTCTTCAGAGGAAGAGACTGATTATCGATCACCGTTTACCCCCCGTCCGTCCCCATCATTTATATGATAGCGAGTCGTATTGTAACGTGGTGACCAGGTGGTTCACATATCGCCCTTGCTCGGCTTGGCCGCCTGCCCCAGACT'
k = 12
#A C G T
profile = """
0.253	0.193	0.169	0.386
0.181	0.229	0.265	0.325
0.265	0.181	0.241	0.313
0.289	0.217	0.241	0.253
0.145	0.289	0.301	0.265
0.217	0.265	0.253	0.265
0.253	0.313	0.229	0.205
0.349	0.277	0.205	0.169
0.229	0.241	0.373	0.157
0.205	0.241	0.169	0.386
0.277	0.241	0.229	0.253
0.229	0.325	0.241	0.205
"""

profile = [[float(f) for f in s.split()] for s in profile.split('\n') if s]
positions = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def get_prob(kmer, profile):
    probability = reduce(lambda x, y: x * y,
                  [profile[j][positions[c]]
                   for j, c in enumerate(list(kmer))])
    return probability

best_kmer, best_prob = text[:k], 0
for i in range(len(text) - k + 1):
    kmer = text[i: i+k]
    prob = get_prob(kmer, profile)
    #print(prob)
    if prob > best_prob:
        best_kmer, best_prob = kmer, prob
#    print '%s: %s' % (kmer, prob)
print (best_prob)
print (best_kmer)