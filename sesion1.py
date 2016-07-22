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

    
