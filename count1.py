# region CHARGING STATION: The Frequency Array
import itertools
#counts how many times certain pattern is observed in the sequence
def pattern_count(text, pattern):
    count = 0
    for i in range(0, len(text) - len(pattern) + 1):
        if text[i : len(pattern) + i] == pattern:
            count += 1
    return count


text1 = "GACCATCAAAACTGATAAACTACTTAAAAATCAGT"
pattern1 = "AAA"
#print (pattern_count(text1, pattern1))

#FrequentWords is shown below.,http://ideone.com/fork/A2kahh
#https://ideone.com/fork/LrznON
def frequentwords(text, k):
    frequentpatterns, count = [], []
    for i in range(len(text) - k + 1):
        count.append(pattern_count(text, text[i:i+k]))
    maxCount = max(count)
    for i in range(len(count)):
        if count[i] == maxCount:
            frequentpatterns.append(text[i:i+k])
    return set(frequentpatterns)
text= "CGCCTAAATAGCCTCGCGGAGCCTTATGTCATACTCGTCCT"
k= 3

#print(frequentwords(text, k))

#function to reverse complement of the string
def reversecomplement(text):
    base = ['A', 'C', 'G', 'T']
   # print (list(reversed(base)))
    #print((base.index('C')))
    #print(list(reversed(base))[base.index('C')])
    result = [list(reversed(base))[base.index(ch)] for ch in reversed(text)]
    return ''.join(result)

text2='CCAGATC'


f = open('myfile','w')
f.write(reversecomplement(text2)) # python will convert \n to os.linesep
f.close() # you can omit in most cases as the destructor will call it
#print (reversecomplement(text2))


#pattern matching, all positions where the pattern starts
def substr(pattern, genome):
    k, result = len(pattern), []
    for i in range(len(genome) - k + 1):
        if genome[i:i+k] == pattern:
            result.append(i)
    #import string
    #string.join( result, '(*)' )
    return result

pattern ='CTTGATCAT'
genome ='aa'
#print(*substr(pattern, genome), sep=' ')

def patterntonumber(pattern):
    patterns = [''.join(x) for x in itertools.product(['A', 'C', 'G', 'T'], repeat=len(pattern))]
    return patterns.index(pattern)

def numbertopattern(number, k):
    patterns = [''.join(x) for x in itertools.product(['A', 'C', 'G', 'T'], repeat=k)]
    return patterns[number]

def findclump(genome1, k, l, t):
    result = set()
    for i in range(len(genome1) - l + 1):
        wordsrange = range(l - k + 1)
        index, count = [patterntonumber(genome1[i+j:i+j+k]) for j in wordsrange], [1 for j in wordsrange]
        index.sort()
        for j in range(1, l - k + 1):
            if index[j] == index[j - 1]:
                count[j] = count[j - 1] + 1
        for j in wordsrange:
            if count[j] >= t:
                result.add(numbertopattern(index[j], k))
    return result

def ClumpFinding(Text, k, L, t):
    """Clump Finding Problem: Find patterns forming clumps in a string.
        Input: A string Genome, and integers k (frequent k-mers), L (cluster/clumps length), and t (time/count).
        Output: All distinct k-mers forming (L, t)-clumps in Genome.
    """
    FrequentPatterns=[]
    for n in range (0,(len(Text)-L)):
        Text_L=Text[n:(n+L)]
        for i in range(0,(len(Text_L)-k)):
            Pattern = Text[i:(i+k)]
            count=PatternCount(Text_L, Pattern)
            if count == t:
                FrequentPatterns.append(Text[i:(i+k)]) #add Text(i, k) to FrequentPatterns
            FrequentPatterns=list(set(FrequentPatterns))  #remove duplicates from FrequentPatterns
    return FrequentPatterns


def PatternCount(Text, Pattern):
    """"The pattern of the motif for the mismatches,
    The input Text is a string and Pattern is also a string
    """
    count= 0
    for i in range (0,(len(Text)-len(Pattern))):
        if Text[i:(i+len(Pattern))] == Pattern:
            count=count+1
    return count

#example
Text="yyyyyyyyyyyyyyy"

k=10
L=514
t=20
print (ClumpFinding(Text, k, L, t))

genome1='CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA'
k=5
l=50
t=4
print( findclump(genome1, k, l, t))

genome2='tttttttttttttttttttttTTTTTTTTTTTT'
k=10
l=579
t=17
print( findclump(genome2, k, l, t))


# region CHARGING STATION: Finding Frequent Words by Sorting
#http://codereview.stackexchange.com/questions/37932/genome-string-clump-finding-problem
# region CHARGING STATION: Solving the Clump Finding Problem
# region CHARGING STATION: The Frequency Array
import itertools


def patterntonumber(pattern):
    patterns = [''.join(x) for x in itertools.product(['A', 'C', 'G', 'T'], repeat=len(pattern))]
    return patterns.index(pattern)

print(patterntonumber('ATGCAA'))

def numbertopattern(number, k):
    patterns = [''.join(x) for x in itertools.product(['A', 'C', 'G', 'T'], repeat=k)]
    return patterns[number]

def computefrequency(text, k):
    result = [0 for i in range(4 ** k)]
    for i in range(len(text) - k + 1):
        result[patterntonumber(text[i:i+k])] += 1
    return result


def fasterfrequentwords(text, k):
    frequentpatterns, frequencyarray = [], computefrequency(text, k)
    maxCount = max(frequencyarray)
    for i in range(4 ** k - 1):
        if frequencyarray[i] == maxCount:
            frequentpatterns.append(text[i:i+k])
    return frequentpatterns

def betterclumpfinding(genome, k, l, t):
    result, clump = set(), [0 for i in range(4 ** k)]
    frequencyarray = computefrequency(genome[:l], k)
    for i in range(4 ** k):
        if frequencyarray[i] >= t:
            clump[i] = 1
    for i in range(1, len(genome) - l + 1):
        frequencyarray[patterntonumber(genome[i-1:i+k-1])] -= 1
        j = patterntonumber(genome[i+l-k:i+l])
        frequencyarray[j] += 1
        if frequencyarray[j] >= t:
            clump[j] = 1
    for i in range(4 ** k):
        if clump[i] == 1:
            result.add(numbertopattern(i, k))
    return result
# endregion

#betterclumpfinding(genome, k, l, t)
genome3='CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA'
k3=5
l3=50
t3=4
print(betterclumpfinding(genome3, k3, l3, t3))

print('i do not know555')

import math
import numpy
#genome='CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA'
#k=5
#L=50
#t=4
def betterClumpFinding(genome,k,L,t):
    frequentPatterns = []
    clump = numpy.zeros( math.pow(4,k), dtype=numpy.int )
    text = genome[0:L]
    frequencyArray = computingFrequencies(text,k)
    for i in range(len(clump)):
        if frequencyArray[i] >= t:
            clump[i] = 1
    for i in range(1,len(genome)-L+1):
        firstPattern = genome[i-1:i-1+k]
        j = patternToNumber(firstPattern)
        frequencyArray[j] -= 1
        lastPattern = genome[i+L-k:i+L]
        j = patternToNumber(lastPattern)
        frequencyArray[j] += 1
        if frequencyArray[j] >= t:
            clump[j] = 1
    for i in range(len(clump)):
        if clump[i] == 1:
            pattern = numberToPattern(i,k)
            frequentPatterns.append(pattern)
    return frequentPatterns

def testClumpFinding():
    f = open("file1", r)
    lines = f.readlines()
    genome = lines[0].rstrip()
    k,L,t = map(int,lines[1].rstrip().split(" "))
    output = betterClumpFinding(genome,k,L,t)
    print(" ".join(map(str, output)))

def computingFrequencies(text, k):
    frequencyArray = numpy.zeros( math.pow(4,k), dtype=numpy.int )
    for i in range(len(text)-k+1):
        pattern = text[i:i+k]
        j = patternToNumber(pattern)
        frequencyArray[j] += 1

    return frequencyArray

def patternToNumber(pattern):
    if len(pattern) == 0:
        return 0
    prefix = pattern[:-1]
    symbol = pattern[-1:]
    return 4*patternToNumber(prefix) + symbolToNumber(symbol)

def symbolToNumber(symbol):
    return {
        'A':0,
        'C':1,
        'G':2,
        'T':3
    }[symbol]

def numberToPattern(number,k):
    if k == 1:
        return numberToSymbol(number)
    quotient, remainder = divmod(number,4)
    prefixPattern = numberToPattern(quotient,k-1)
    symbol = numberToSymbol(remainder)
    return prefixPattern+symbol

def numberToSymbol(number):
    return {
        0:'A',
        1:'C',
        2:'G',
        3:'T'
    }[number]

pattern ='ATGCAA'

print(patternToNumber(pattern))

seq = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
n = 4
k = 1
def countKmers(seq,n,k):
    kmers={}
    kfinals={}
    for i in range(len(seq)-k+1):
        kmer=seq[i:(i+k)]
        if(kmers.has_key(kmer)):
            kmers[kmer]+=1
        else:
            kmers[kmer]=1
    for i in kmers.keys():
        if(kmers[i]>=n):
            kfinals[i]=kmers[i]
    return kfinals
print (countKmers(seq,n,k))
