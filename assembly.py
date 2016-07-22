# -*- coding: utf-8 -*-
def longest_common_subsequence(v, w):
    '''Returns the longest longest common subsequence of strings v and w.'''
    # Initialize the array S and iterate through all character of v and w.
    S = [[0]*(len(w)+1) for _ in range(len(v)+1)]
    for i in range(len(v)):
        for j in range(len(w)):
            if v[i] == w[j]:
                S[i+1][j+1] = S[i][j]+1
            else:
                S[i+1][j+1] = max(S[i+1][j],S[i][j+1])

    # Recover a maximum substring.
    longest_sseq = ''
    i,j = len(v), len(w)
    while i*j != 0:
        if S[i][j] == S[i-1][j]:
            i -= 1
        elif S[i][j] == S[i][j-1]:
            j -= 1
        else:
            longest_sseq = v[i-1] + longest_sseq
            i -= 1
            j -= 1

    return longest_sseq
v=" AGACTG"
w="GTACGA"
print (*longest_common_subsequence(v, w),sep=' ')

def frequentWords( s, k ):
    counts = {}
    for i in range(len(s)-k+1):
        if s[i:i+k] not in counts:
            counts[s[i:i+k]] = 0
        counts[s[i:i+k]] += 1
    m = max(counts.values())
    out = []
    for kmer in counts:
        if counts[kmer] == m:
            out.append(kmer)
    return out

s = 'Aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa'
k = 14
print (*frequentWords(s,k),sep=' ')
print ("JUST BLANK LINE")

# region CHARGING STATION: The Frequency Array
import itertools
#counts how many times certain pattern is observed in the sequence
def pattern_count(text, pattern):
    count = 0
    for i in range(0, len(text) - len(pattern) + 1):
        if text[i : len(pattern) + i] == pattern:
            count += 1
    return count
A = "GACCATCAAAACTGATAAACTACTTAAAAATCAGT"
B = "AAA"
print (pattern_count(A,B))


text= "CGCCTAAATAGCCTCGCGGAGCCTTATGTCATACTCGTCCT"
k= 3

def reversecomplement(text):
    base = ['A', 'C', 'G', 'T']
   # print (list(reversed(base)))
    #print((base.index('C')))
    #print(list(reversed(base))[base.index('C')])
    result = [list(reversed(base))[base.index(ch)] for ch in reversed(text)]
    return ''.join(result)

text2='GCTAGCT'
print("REVERSE COMPLEMENT")
print(reversecomplement(text2))

def frequentwords(text, k):
    frequentpatterns, count = [], []
    for i in range(len(text) - k + 1):
        count.append(pattern_count(text, text[i:i+k]))
    maxCount = max(count)
    for i in range(len(count)):
        if count[i] == maxCount:
            frequentpatterns.append(text[i:i+k])
    return set(frequentpatterns)
print(frequentwords(text, k))




def assembly(text,k):
    return sorted([text[i:i+k] for i in range(len(text)-k+1)])

#print(assembly('CAATCCAAC',5))
print(*assembly('CAATCCAAC',5), sep=' ')

def frequentWordsWithMismatches( s, k, d ):
    counts = {}
    for i in range(len(s)-k+1):
        for neighbor in neighbors(s[i:i+k],d):
            if neighbor not in counts:
                counts[neighbor] = 0
            counts[neighbor] += 1
    m = max(counts.values())
    return [kmer for kmer in counts if counts[kmer] == m]

def neighbors( s, d ):
    if d == 0:
        return [s]
    if len(s) == 1:
        return ['A','C','G','T']
    out = []
    for neighbor in neighbors(s[1:],d):
        if hamming(s[1:],neighbor) < d:
            out.extend(['A'+neighbor,'C'+neighbor,'G'+neighbor,'T'+neighbor])
        else:
            out.append(s[0] + neighbor)
    return out

def hamming( s, t ):
    return sum([s[i] != t[i] for i in range(len(s))])

s = 'aaaaa'
k = 5
d = 3
print (*frequentWordsWithMismatches(s,k,d),sep=' ')

#https://github.com/niemasd/Algorithm-Problem-Solutions/blob/79e248e9b1d824b95b7c000589ac9163861fcada/ROSALIND%20Solutions/Bioinformatics%20Textbook%20Track/1G%20-%20Frequent%20Words%20with%20Mismatches%20Problem/FrequentWordsWithMismatchesProblem.py

def frequentWordsWithMismatchesAndReverseComplements( s, k, d ):
    counts = {}
    for i in range(len(s)-k+1):
        for sub in [s[i:i+k],reverseComplement(s[i:i+k])]:
            for neighbor in neighbors(sub,d):
                if neighbor not in counts:
                    counts[neighbor] = 0
                counts[neighbor] += 1
    m = max(counts.values())
    return [kmer for kmer in counts if counts[kmer] == m]

def neighbors( s, d ):
    if d == 0:
        return [s]
    if len(s) == 1:
        return ['A','C','G','T']
    out = []
    for neighbor in neighbors(s[1:],d):
        if hamming(s[1:],neighbor) < d:
            out.extend(['A'+neighbor,'C'+neighbor,'G'+neighbor,'T'+neighbor])
        else:
            out.append(s[0] + neighbor)
    return out
print ("HAMMING")
z = "CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA"
y = "CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG"
def hamming( s, t ):
    return sum([s[i] != t[i] for i in range(len(s))])
print(hamming(z,y))

def reverseComplement( s ):
    return ''.join([complement(s[i]) for i in range(len(s)-1,-1,-1)])

def complement( s ):
    return {'A':'T','T':'A','C':'G','G':'C'}[s]

s = 'ATGATTTTTATGAAATTATGATGATGCAATTCAATGTTTATTATTTTTAAATGTTTATGATTATGATTTTTAAATTAACAATGATTATGAATTTCAATGTTTATTATTAACAATGATTTTTTTTCACACAATGATTTTTTTTCAATGATTCAATGTTTATTATTATTATTCAATGATGCAATTCATTTATGTTTATGTTTTTTATTCAAAATTATTTTTATGAAATTAAAAATTATGTTTTTT'
k = 5
d = 3
print (*frequentWordsWithMismatchesAndReverseComplements(s,k,d),sep = '\n')
print(*neighbors('ATCATCACT',2), sep = '\n')


def findPattern(pattern,seq, d=3):
    pos={}
    for i in range(0, len(seq)-len(pattern)+1):
        s=seq[i:(i+len(pattern))]
        nmismatch=0
        for j in range(len(s)):
            if(s[j]!=pattern[j]):
                nmismatch=nmismatch+1
                if(nmismatch>d):
                    break
        if(nmismatch<=d):
            pos[i]=s
    return(pos)

print (findPattern('ATTCTGGA', 'CGCCCGAATCCAGAACGCATTCCCCTGGCCTCCATTCTGGAACGGTACGGACGTCAATCAAAT',3))


#with open('data/stepic_1f.txt') as input_data:
	#pattern, dna, n = [line.strip() if index != 2 else int(line.strip()) for index, line in enumerate(input_data.readlines())]
print('COUNTING:')
dna = 'CATGCCATTCGCATTGTCCCAGTGA'
pattern ="CCC"
n = 2
approx_match = []
for i in range(len(dna)-len(pattern)+1):
	mismatch_count = 0
	for j in range(len(pattern)):
		if dna[i:i+len(pattern)][j] != pattern[j]:
			mismatch_count += 1

	if mismatch_count <= n:
		approx_match.append(str(i))

#print ' '.join(approx_match)
print (*approx_match,sep = ' ')
#with open('output/Assignment_01F.txt', 'w') as output_data:
	#output_data.write(' '.join(approx_match))
print("PROBLEM 5")
s = "ACGT"
k = 4
d = 3
def frequentWordsWithMismatches( s, k, d ):
    counts = {}
    for i in range(len(s)-k+1):
        for neighbor in neighbors(s[i:i+k],d):
            if neighbor not in counts:
                counts[neighbor] = 0
            counts[neighbor] += 1
    m = max(counts.values())
    return [kmer for kmer in counts if counts[kmer] == m]
print(*frequentWordsWithMismatches( s, k, d ),sep = '\n')



#with open('data/stepic_1e.txt') as input_data:
#	dna = input_data.read().strip()
dna ="CATTCCAGTACTTCATGATGGCGTGAAGA"
skew_value, min_skew, min_ind = 0, 1, []
for index, nucleotide in enumerate(dna):
	# Determine the skew value.
	if nucleotide == 'C':
		skew_value -= 1
	elif nucleotide == 'G':
		skew_value += 1
	# Check if it matches the current minimum, or is a new minimum.
	if skew_value == min_skew:
		min_ind.append(str(index+1))
	elif skew_value > min_skew:
		min_skew = skew_value
		min_ind = [str(index+1)]
print('SKEW PROBELM 4')

print (min_ind)


#%% -------------------- SKEW ---------------------
#Retornar una cadena numérica que vaya acumulando las G y restando las C
seq="CATTCCAGTACTTCATGATGGCGTGAAGA"
def skew(seq):
    count=0
    sk=[]
    sk.append(count)
    for i in range(len(seq)):
        if(seq[i]=='G'):
            count=count+1
        elif(seq[i]=='C'):
            count=count-1
        sk.append(count)
    return(sk)

#%% Pruebas

#%% -------------------- MINIMUM SKEW ---------------------
#Retornar la posición o posiciones donde el skew es mínimo
def minimumSkew(seq):
    sk=skew(seq)
    mins=[]
    msk=max(sk)
    for i in range(len(sk)):
        if(sk[i]==msk):
            mins.append(i)
    return(mins)

#%% Pruebas
#seq= 'CATTCCAGTACTTCATGATGGCGTGAAGA'
#sk=skew(seq)
print (minimumSkew(seq))