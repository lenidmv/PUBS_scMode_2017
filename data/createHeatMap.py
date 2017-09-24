#!/usr/bin/python
import datetime
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pickle
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.ticker import FixedLocator

#load data
data = pickle.load( open( "allele_dic_with_WT.pkl", "rb" ) )
codondict = pickle.load(open("translate.pkl", "rb"))
aminoacidpos = pickle.load(open("aminotonumber.pkl", "rb"))
aminostring = []

errork = []
lengths = 0
#get amino acids
for key, value in sorted(aminoacidpos.iteritems(), key=lambda (k,v): (v,k)):
    aminostring.append(key)
#labels for y axis
#get rid of places where there aren't 3 aa in codons and get protein length

for key,value in data.iteritems():
    if len(value[1]) != 3:
        errork.append(key)
    if value[0] > lengths:
        lengths = value[0]

positions = range(1,lengths)
#create output array
numouts=lengths*21
x = np.repeat(0,numouts+1)
x = np.resize(x,(21,lengths+1))

#check errors
for item in errork:
    del data[item]


#translate to amino acids
counter = []
for k,v in sorted(data.values()):
    newcodon = ''
    #translate DNA to mRNA
    for i in range(len(v)):
        if v[i] == 'G':
            newcodon+='C'
        if v[i] == 'T':
            newcodon+='A'
        if v[i] == 'A':
            newcodon+='U'
        if v[i] == 'C':
            newcodon+='G'
    if newcodon in codondict:
        aa = codondict[newcodon]
    #print k, v
    temptup = (k,aa)
    counter.append(temptup)
#print sorted(counter)
counter = sorted(counter)
#print counter
#frequency variables for number amino acids seen, count number of appearences
#this could probably be condensed, but I'm being lazy
STOP, W, F, Y,L,I,M,V,C,A, G,P,S,T,N,Q,H,R,K,D,E = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
iterate = 1
for i in range(len(counter)):
    if counter[i][0] == iterate:
        if counter[i][1] == 'W':
            W+=1
            x[1][iterate] = W
        if counter[i][1] == 'F':
            F+=1
            x[2][iterate] = F
        if counter[i][1] == 'Y':
            Y+=1
            x[3][iterate] = Y
        if counter[i][1] == 'L':
            L+=1
            x[4][iterate]=L
        if counter[i][1] == 'I':
            I+=1
            x[5][iterate]=I
        if counter[i][1] == 'M':
            M+=1
            x[6][iterate]=M
        if counter[i][1] == 'V':
            V+=1
            x[7][iterate]=V
        if counter[i][1] == 'C':
            C+=1
            x[8][iterate]=C
        if counter[i][1] == 'A':
            A+=1
            x[9][iterate]=A
        if counter[i][1] == 'G':
            G+=1
            x[10][iterate]=G
        if counter[i][1] == 'P':
            P+=1
            x[11][iterate]=P
        if counter[i][1] == 'S':
            S+=1
            x[12][iterate]=S
        if counter[i][1] == 'T':
            T+=1
            x[13][iterate]=T
        if counter[i][1] == 'N':
            N+=1
            x[14][iterate]=N
        if counter[i][1] == 'Q':
            Q+=1
            x[15][iterate]=Q
        if counter[i][1] == 'H':
            H+=1
            x[16][iterate]=H
        if counter[i][1] == 'R':
            R+=1
            x[17][iterate]=R
        if counter[i][1] == 'K':
            K+=1
            x[18][iterate]=K
        if counter[i][1] == 'D':
            D+=1
            x[19][iterate] = D
        if counter[i][1] == 'E':
            E+=1
            x[20][iterate]=E
        if counter[i][1] == 'STOP':
            STOP+=1
            x[0][iterate]=STOP
    #is this skipped?
    elif counter[i][0] != iterate:
        iterate+=1
        STOP, W, F, Y,L,I,M,V,C,A, G,P,S,T,N,Q,H,R,K,D,E = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0

##ploting
fig, ax = plt.subplots()
heatmap = ax.pcolor(x, cmap='hot')

#legend
cbar = plt.colorbar(heatmap)
highest = np.amax(x)
lowest = np.amin(x)
cbar.ax.get_yaxis().set_ticks([])
#plot sidebar
for j, lab in enumerate([highest*1/6,highest*2/6, highest*4/6,">"+str(highest*5/6)]):
    cbar.ax.text(.5, (2 * j + 1) / 8.0, lab, ha='center', va='center')
cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_ylabel('# of barcodes', rotation=270)

majorLocator   = FixedLocator(np.linspace(0,lengths,5))
ax.set_yticks(np.arange(x.shape[0]), minor=False)
ax.invert_yaxis()

#lebels
row_labels = aminostring
#column_labels = list(range(0,77,5))
column_labels=positions
ax.set_xticklabels(column_labels, minor=True)
ax.set_yticklabels(row_labels, minor=False)
ax.set_xlabel('protein position')
ax.set_ylabel('amino acid')
ax.axis('tight')
fig.suptitle('Frequency of Amino Acids Labeled by Barcodes', fontsize=14, fontweight='bold')

plt.show()
