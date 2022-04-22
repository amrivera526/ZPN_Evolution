#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
from string import digits

import numpy
import scipy
#import Image
#import kiwisolver
import matplotlib #Had to reinstall Pillow 
import matplotlib.pyplot as plt
import seaborn

import re

import sys
import pandas
import random
import math

import sklearn
#from sklearn.model_selection import metrics
from sklearn import linear_model
from sklearn.model_selection import cross_validate
#from sklearn.cross_validation import cross_val_score
#from sklearn.model_selection import metrics

import time,datetime
from mpl_toolkits import mplot3d
from matplotlib.text import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
from matplotlib import cm

from operator import itemgetter
#import logomaker
import pickle

# In[2]:


sorted(numpy.arange(0.01,0.05,0.002),reverse=True)


# In[3]:


class Sequence:
    
    def __init__(self):
        self.seq=""
        self.name=""
        self.file=""
        self.alngroup=""
        
    def GiveName(self,N):
        self.name = N
    
    def GiveSeq (self,SEQ):
        self.seq = SEQ
    
    def GiveFile(self,F):
        self.file = F
        
    def GiveAln(self,A):
        self.alngroup = A
    


# In[4]:


def ReadIn(f,GeneTag=''):
    S=[]
    total = 0
    unique = 0
    

    blump = 0
    fin = open(f,'r')
    #line = fin.readline()
    for line in fin:
        if line[0] == '>':
            if blump > 0:
                AddHere = Sequence()
                AddHere.GiveSeq(nuc)
                AddHere.GiveName(NAME)
                AddHere.GiveFile(f)
                S.append(AddHere)
                #seqs[name] = nuc #Push previous sequence
            total+=1 
            nuc = ''
            blump = 1
            prename = line[1:].replace(' ','_').strip('\n').strip('\r')

           
            #N=name.split('_')
            #print(N)
            #N0=N[0]
            #N1=N[1]
            #N1="ENS"+N1#Should be more similar to prename from previous version
            if GeneTag !='':
                prename=GeneTag+'_'+prename
            NAME=prename


            #print(prename.translate(table))
            #if N1.translate(table)[:-1] in ensembl.keys():
             #   name = ensembl[N1.translate(table)[:-1]]+'_' + N0 + '_' + N1
            #else:
             #   name = N0 + '_' + N1
        else:
            nuc+=line.strip('\n').strip('\r')
    #if name not in seqs.keys():
        #seqs[name] = nuc #Push the final sequence
    AddHere = Sequence()
    AddHere.GiveSeq(nuc)
    AddHere.GiveName(NAME)
    AddHere.GiveFile(f)
    S.append(AddHere)
    fin.close()
    return(S)

   


# In[5]:


def SimpMult(A,B):
    X=0.0
    for i in numpy.arange(len(A)):
        X+=A[i]*B[i]
        
    return(X)
    


# In[6]:


numpy.exp(1)


# In[7]:


def LogitEst(V1,V2,c):
    polynom = SimpMult(V1,V2)+c
    Ans=1/(1+numpy.exp(-1*polynom))
    return(Ans)


# In[8]:


#codes = ['-','A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
#         'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

#def create_dict(codes):
#    char_dict = {}
#    for index, val in enumerate(codes):
#        char_dict[val] = index+1

#    return char_dict

#char_dict = create_dict(codes)

#char_dict


# In[9]:

'''
#Cite this function?
def draw_logo(all_scores, fontfamily='Arial', size=80):
    if fontfamily == 'xkcd':
        plt.xkcd()
    else:
        matplotlib.rcParams['font.family'] = fontfamily

    fig, ax = plt.subplots(figsize=(len(all_scores), 2.5))

    font = FontProperties()
    font.set_size(size)
    font.set_weight('bold')
    
    #font.set_family(fontfamily)

    ax.set_xticks(range(1,len(all_scores)+1))    
    ax.set_yticks(range(0,3))
    ax.set_xticklabels(range(1,len(all_scores)+1), rotation=90)
    ax.set_yticklabels(numpy.arange(0,3,1))    
    seaborn.despine(ax=ax, trim=True)
    
    trans_offset = transforms.offset_copy(ax.transData, 
                                          fig=fig, 
                                          x=1, 
                                          y=0, 
                                          units='dots')
   
    for index, scores in enumerate(all_scores):
        yshift = 0
        for base, score in scores:
            txt = ax.text(index+1, 
                          0, 
                          base, 
                          transform=trans_offset,
                          fontsize=80, 
                          color=COLOR_SCHEME[base],
                          ha='center',
                          fontproperties=font,

                         )
            txt.set_path_effects([Scale(1.0, score)])
            fig.canvas.draw()
            window_ext = txt.get_window_extent(txt._renderer)
            yshift = window_ext.height*score
            trans_offset = transforms.offset_copy(txt._transform, 
                                                  fig=fig,
                                                  y=yshift,
                                                  units='points')
        trans_offset = transforms.offset_copy(ax.transData, 
                                              fig=fig, 
                                              x=1, 
                                              y=0, 
                                              units='points')    
    plt.show()

'''
# In[ ]:





# In[10]:


#os.chdir(r'C:\Users\alber\OneDrive\Documents\Swanson_2017Onward\Sequences\ZP_Alignments')

#Aln = ReadIn('RepresentativeZPN_CryStrAlignedSimpler_202009141616.fasta')
#Aln=ReadIn('AlgoAln_ZPN_90c_20200625_75Train_202009241307.fasta')
#Aln=ReadIn('AlgoAln_ZPN_90c_20200625_75Train_202011021352.fasta')
#Aln = ReadIn('AlgoAln_ZPN_PSI90c_Train_202111031718.fasta')
#Aln = ReadIn('AlgoAln_ZPN_PSI90c_MAFFTProMals_Train_202112140138.fasta')

Aln = ReadIn('AlgoAln_ZPN_PSI100c_MAFFTProMals_Train_202203242127.fasta')

# In[ ]:





# In[11]:


#Basic mapper
residues = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
           'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

mapper = dict(zip(residues, [[1.0 if residues[i] == r else 0.0 for i in range(len(residues))] for r in residues]))
mapper['-'] = [0.0]*20


# In[12]:

'''
PCADat = open('AminoAcid_8PCACodes.txt','r')

mapper2={}
for pca in PCADat:
    SP = pca.split('\t')
    AAn=(SP[0][len(SP[0])-1])
    V =[]
    for i in numpy.arange(len(SP)):
        if i>0:
            Vi = float(SP[i].strip('\n'))
            V.append(Vi)
    #print(AAn,V)
    mapper2[AAn]=V

mapper2['-']=[0.0]*8 #Must make by 8. Previously had 5
PCADat.close()
mapper2
'''




# In[13]:


def ConvertRes(seq,mapperx):
    #print(seq)
    results = []
    for i in range(len(seq)):
        results += mapperx[seq[i]]
    return(results)


# In[14]:


random.seed(2112)
#Randomize BEts

NumPam = 20*len(Aln[0].seq)
print(NumPam)

InitBetas=numpy.random.normal(0,0.1,NumPam) #Let Params vary since logistic functions can take input from -inf to inf


# In[15]:


#len(InitBetas)

# Assign #labels


LastZPN = ["ZP1ZPN2","ZP2ZPN4","ZP3ZPN1","ZP4ZPN2","ZPAXLastZPN","UMODZPN","TECTAZPN","CUZD1ZPN","ZPDZPN"]

NameList=[]


for i in numpy.arange(len(Aln)):
    NS = Aln[i].name.split('_')
    #print(NS)
    NameList.append(NS[0]+NS[1])
    
#print(NameList)


GroundTruth = [1.0 if NameList[i] in LastZPN else 0.0 for i in range(len(NameList))]

GroundTruth


# In[ ]:





# In[ ]:





    


# In[ ]:





# In[16]:


CostLim=0.1
lamb = 0.1
#lamb=0
lamb=lamb/len(Aln)
Betas = InitBetas

#limit condition
ite=0
CFList=[]

CF=10

alph0=0.1
alph=alph0

b0=random.uniform(-0.1,0.1)


# In[17]:


AllDat=[]

for i in numpy.arange(len(Aln)):
    #print(ConvertRes(Aln[i].seq))
    XLine=ConvertRes(Aln[i].seq,mapper) #Decide which mapper you want
    #XLine.append(GroundTruth[i])
    #print(XLine)
    AllDat.append(XLine)

print(len(AllDat[0]))


# In[18]:


#AllDat[0]


# In[ ]:





# In[ ]:





# In[19]:


#AlVal = ReadIn('AlgoAln_ZPN_90c_20200625_25Test_202009241307.fasta')
#AlVal = ReadIn('AlgoAln_ZPN_90c_20200625_25Test_202011021352.fasta')
#AlVal = ReadIn('AlgoAln_ZPN_PSI90c_Test_202111031718.fasta')
#AlVal = ReadIn('AlgoAln_ZPN_PSI90c_MAFFTProMals_Test_202112140138.fasta')
AlVal = ReadIn('AlgoAln_ZPN_PSI100c_MAFFTProMals_Test_202203242127.fasta')

ValDat=[]
NL=[]

for j in numpy.arange(len(AlVal)):
    YLine = ConvertRes(AlVal[j].seq,mapper) #Decide which mapper you want
    ValDat.append(YLine)
    
    NLi = AlVal[j].name.split('_')
    NL.append(NLi[0]+NLi[1])
    
    
    
GTV = [1.0 if NL[i] in LastZPN else 0.0 for i in range(len(NL))]
#GTV
#print(numpy.ndarray(ValDat))


# In[20]:


current_time = time.time()


# In[21]:


#[1,0.5,0.2]/10


# In[ ]:





# In[22]:


#LRs = [0,0.25,0.5,0.75,1]
#LRs=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]

LR1=numpy.arange(0,1,0.02)
LR2=numpy.arange(0.4,0.8,0.02) #ZOomed in

LR3 = numpy.arange(0.25,0.75,0.01)
#LRs = [0.01*i for i in numpy.arange(101)]
#Alphs = [0,0.001,0.01,0.1,1]
Alphs=[1]
#Cs =[1000,100,10,1]
C1=[1,0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001] # Skip first 4 C's In the 0.01 to the 0.05 range
C2=sorted(numpy.arange(0.01,0.05,0.002),reverse=True) #Going down from 0.05 to 0.01

C3 = sorted(numpy.arange(0.01,0.1,0.002),reverse=True)

NumPam=[]
Acc=[]

GT = numpy.reshape(GroundTruth,(1,len(GroundTruth)))
GT=GT[0]



# In[ ]:



reg1 = linear_model.LogisticRegressionCV(penalty='elasticnet',solver='saga',l1_ratios=LR1,max_iter=10000,Cs=C1,scoring='neg_mean_squared_error').fit(AllDat,GT)
reg2 = linear_model.LogisticRegressionCV(penalty='elasticnet',solver='saga',l1_ratios=LR2,max_iter=10000,Cs=C2,scoring='neg_mean_squared_error').fit(AllDat,GT)
reg3 = linear_model.LogisticRegressionCV(penalty='elasticnet',solver='saga',l1_ratios=LR3,max_iter=10000,Cs=C3,scoring='neg_mean_squared_error').fit(AllDat,GT)


#pickle.dump(reg1, open('ZPNPSI90c_MAFFTPromalAF2_Reg1_NoZoom_202112140156.sav', 'wb'))
#pickle.dump(reg2, open('ZPNPSI90c_MAFFTPromalAF2_Reg2_NoZoom_202112140156.sav', 'wb'))
#pickle.dump(reg3, open('ZPNPSI90c_MAFFTPromalAF2_Reg3_NoZoom_202112140156.sav', 'wb'))

pickle.dump(reg1, open('ZPNPSI100c_MAFFTPromalAF2_Reg1_202203242152.sav', 'wb'))
pickle.dump(reg2, open('ZPNPSI100c_MAFFTPromalAF2_Reg2_202203242152.sav', 'wb'))
pickle.dump(reg3, open('ZPNPSI100c_MAFFTPromalAF2_Reg3_202203242152.sav', 'wb'))