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
genome = 'AAAAATTTTTTTTGGGGGGGGCATA'
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
print('i do not know555')

# region CHARGING STATION: Finding Frequent Words by Sorting
#http://codereview.stackexchange.com/questions/37932/genome-string-clump-finding-problem
# region CHARGING STATION: Solving the Clump Finding Problem
# region CHARGING STATION: The Frequency Array
import itertools


def patterntonumber(pattern):
    patterns = [''.join(x) for x in itertools.product(['A', 'C', 'G', 'T'], repeat=len(pattern))]
    return patterns.index(pattern)


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
genome='CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA'
print(betterclumpfinding(genome, k, l, t))


# -*- coding: utf-8 -*-
"""
BIOINFORMÁTICA MUII

@author: rodri
"""
#%% ---------------- Lectura de archivos ---------------
import urllib2
f=urllib2.urlopen("http://vis.usal.es/rodrigo/documentos/bioinfo/avanzada/genomas/Vibrio_cholerae.txt", "r")
vc=f.readlines()    #toma todas las líneas del archivo y las mete en una lista
vc=vc[0]            #el fichero de genoma sólo tiene una línea
import string
vc=string.replace(vc, "\n", "")
print len(vc)          #(muy larga, claro)
f.close()
#%%
f=urllib2.urlopen("http://vis.usal.es/rodrigo/documentos/bie/E-coli.txt", "r")
ec=f.readlines()    #toma todas las líneas del archivo y las mete en una lista
ec=ec[0]            #el fichero de genoma sólo tiene una línea
import string
ec=string.replace(ec, "\n", "")
print len(ec)          #(muy larga, claro)
f.close()

#%%---------------------- COMPLEMENTARIA ----------------------
def complement(seq):
    comp=[]
    us=seq.upper()
    for s in us:
        if(s=='A'):    comp.append('T')
        elif(s=='T'):  comp.append('A')
        elif(s=='G'):  comp.append('C')
        elif(s=='C'):  comp.append('G')
    return('%s' % ''.join(map(str, comp)))
complement("ACGTTCGT")
#%%---------------------- INVERSA ----------------------
def reverse(seq):
    rev=[]
    for i in range(len(seq)-1,-1,-1):
        rev.append(seq[i])
    return('%s' % ''.join(map(str, rev)))
reverse("ACGTTCGT")
#%%----------- INVERSA COMPLEMENTARIA -----------
def revcomp(seq):
    return complement(reverse(seq))
print revcomp("ACGTTCGT")

print revcomp(vc[((int)(len(vc)*0.5)-10):((int)(len(vc)*0.5)+10)])

#%% ---------------- FRECUENCIA ---------------
def freq(seq):
    f={'A':0,'C':0,'G':0,'T':0}
    for x in seq:
        if(x in f.keys()):
            f[x]=f[x]+1
    for k in f.keys():
        f[k]=(float)(f[k])/len(seq)
    return f
print freq("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT")
print freq("ACGTTCGT")

freqEC=freq(ec)
#%%------------ GC CONTENT ---------------
freqEC["G"]+freqEC["C"]

#%%------------ GC CONTENT de distintos organismos
def GCcontent(url):
    import urllib2
    f=urllib2.urlopen(url, "r")
    lines=f.readlines()    #toma todas las líneas del archivo y las mete en una lista
    f.close()

    #f=open(path, "r")
    #lines=f.readlines()    #toma todas las líneas del archivo y las mete en una lista
    #f.close()
    gc=0
    totSeqs=0
    for line in lines:
        if(line[0]!='>'):
            f=freq(string.replace(line, "\n", ""))
            gc+=f['C']+f['G']
            totSeqs+=1
    return gc/totSeqs
print GCcontent("http://vis.usal.es/rodrigo/documentos/bie/C_albicans_SC5314.fasta")
#%%
import time
t0=time.clock()
print "Contenido: {}%".format(100*GCcontent("http://vis.usal.es/rodrigo/documentos/bie/Plasmodium_falciparum.fa"))
print "Tarda {} s".format((time.clock()-t0))	#resultado en segundos
#%%
t0=time.clock()
print "Contenido: {}%".format(100*GCcontent("http://vis.usal.es/rodrigo/documentos/bie/Streptomyces_coelicolor.fa"))
print "Tarda {} s".format((time.clock()-t0))	#resultado en segundos



#%% -------------------- SKEW ---------------------
#Retornar una cadena numérica que vaya acumulando las G y restando las C
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
seq="CATGGGCATCGGCCATACGCC"
sk=skew(seq)
sk

#%% Min skew
#Retornar la posición o posiciones donde el skew es mínimo
def minimumSkew(seq):
    sk=skew(seq)
    mins=[]
    msk=min(sk)
    for i in range(len(sk)):
        if(sk[i]==msk):
            mins.append(i)
    return(mins)

#%% Pruebas
seq="TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"
sk=skew(seq)
minimumSkew(seq)


#%% Gráficas
#Mostrando de manera gráfica el skew
import matplotlib.pyplot as plt
plt.plot(sk, label='skew')
plt.ylabel('skew')
plt.xlabel('position')
plt.show()

#%% Con Vibrio cholerae
sk=skew(vc)
print minimumSkew(vc)
plt.plot(sk, label='skew')
plt.ylabel('skew')
plt.xlabel('position')
plt.show()

#%% Con E coli
sk=skew(ec)
ms=minimumSkew(ec)
plt.plot(sk, label='skew')
plt.ylabel('skew')
plt.xlabel('position')
plt.show()

oricEC=ec[ms[0]:(ms[0]+500)]

#%% Ahora con Thermotoga!
import urllib2
f=urllib2.urlopen("http://vis.usal.es/rodrigo/documentos/bioinfo/avanzada/genomas/Thermotoga-petrophila.txt", "r")
tg=f.readlines()    #toma todas las líneas del archivo y las mete en una lista
tg=tg[0]            #el fichero de genoma sólo tiene una línea


#%%
sk=skew(tg)
print minimumSkew(tg)
plt.plot(sk, label='skew')
plt.title="Thermotoga petrophila"
plt.ylabel('skew')
plt.xlabel('position')
plt.show()


