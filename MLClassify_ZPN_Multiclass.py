#!/usr/bin/env python
# coding: utf-8

# In[37]:


import os
from string import digits

import numpy as np
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
import logomaker
import pickle

# In[2]:


sorted(np.arange(0.01,0.05,0.002),reverse=True)


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
    for i in np.arange(len(A)):
        X+=A[i]*B[i]
        
    return(X)
    


# In[6]:


np.exp(1)


# In[7]:


def LogitEst(V1,V2,c):
    polynom = SimpMult(V1,V2)+c
    Ans=1/(1+np.exp(-1*polynom))
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
    ax.set_yticklabels(np.arange(0,3,1))    
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
    for i in np.arange(len(SP)):
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


#Let Params vary since logistic functions can take input from -inf to inf


# In[15]:


#len(InitBetas)

# Assign #labels



#Aln = ReadIn('RepresentativeZPN_CryStrAlignedSimpler_202009141616.fasta')
#Aln=ReadIn('AlgoAln_ZPN_90c_20200625_75Train_202009241307.fasta')
#Aln=ReadIn('AlgoAln_ZPN_90c_20200625_75Train_202011021352.fasta')

#Aln = ReadIn('AlgoAln_ZPN_PSI90c_MAFFTProMals_Train_202112140138.fasta')


Aln = ReadIn('AlgoAln_ZPN_PSI100c_MAFFTProMals_Train_202203242127.fasta')
random.seed(2112)
#Randomize BEts

NumPam = 20*len(Aln[0].seq)
print(NumPam)

InitBetas=np.random.normal(0,0.1,NumPam)


LastZPN = ["ZP1ZPN2","ZP2ZPN4","ZP3ZPN1","ZP4ZPN2","ZPAXLastZPN","UMODZPN","TECTAZPN","CUZD1ZPN","ZPDZPN"]

FirstZPN=["ZP1ZPN1","ZP2ZPN1","ZP4ZPN1","ZPAXZPN1"]


NameList=[]


for i in np.arange(len(Aln)):
    NS = Aln[i].name.split('_')
    #print(NS)
    NameList.append(NS[0]+NS[1])
    
#print(NameList)

A2=[] #Start with this empty!

NL2=[]


#for j in np.arange(len(Aln)):
   # if NameList[j] in LastZPN:
   #     NL2.append(NameList[j])
   #     A2.append(Aln[j])
#NameList=NL2
#Aln=A2




GroundTruth = [0.0 if NameList[i] in FirstZPN else 2.0 if NameList[i] in LastZPN else 1.0 for i in range(len(NameList))]
#print(len(GroundTruth))
GroundTruth    


# In[16]:


sum(GroundTruth)


# In[ ]:





# In[ ]:





    


# In[ ]:





# In[17]:


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


# In[18]:


AllDat=[]

for i in np.arange(len(Aln)):
    #print(ConvertRes(Aln[i].seq))
    XLine=ConvertRes(Aln[i].seq,mapper) #Decide which mapper you want
    #XLine.append(GroundTruth[i])
    #print(XLine)
    AllDat.append(XLine)


# In[19]:


#AllDat[0]


# In[ ]:





# In[ ]:





# In[20]:


#AlVal = ReadIn('AlgoAln_ZPN_90c_20200625_25Test_202009241307.fasta')
#AlVal = ReadIn('AlgoAln_ZPN_90c_20200625_25Test_202011021352.fasta')

#AlVal = ReadIn('AlgoAln_ZPN_PSI90c_MAFFTProMals_Test_202112140138.fasta')
AlVal = ReadIn('AlgoAln_ZPN_PSI100c_MAFFTProMals_Test_202203242127.fasta')


ValDat=[]
NL=[]

for j in np.arange(len(AlVal)):
    YLine = ConvertRes(AlVal[j].seq,mapper) #Decide which mapper you want
    ValDat.append(YLine)
    
    NLi = AlVal[j].name.split('_')
    NL.append(NLi[0]+NLi[1])
    

A3=[]

NL3=[]
V3=[]



GTV = [0.0 if NL[i] in FirstZPN else 2.0 if NL[i] in LastZPN else 1.0 for i in range(len(NL))]
#GTV
#print(np.ndarray(ValDat))


# In[21]:


current_time = time.time()


# In[22]:


#[1,0.5,0.2]/10


# In[ ]:





# In[23]:


LR1=np.arange(0,1,0.02)
LR2=np.arange(0.4,0.8,0.02) #ZOomed in


#LRs = [0.01*i for i in np.arange(101)]
#Alphs = [0,0.001,0.01,0.1,1]
Alphs=[1]
#Cs =[1000,100,10,1]
#C1=[1,0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001] # Skip first 4 C's In the 0.01 to the 0.05 range
#C1=sorted(np.arange(0.01,0.05,0.002),reverse=True)
#Going down from 0.05 to 0.01
C1=[1,0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001] # Skip first 4 C's In the 0.01 to the 0.05 range
C2=sorted(np.arange(0.01,0.05,0.002),reverse=True) #Going down from 0.05 to 0.01

NumPam=[]
Acc=[]

GT = np.reshape(GroundTruth,(1,len(GroundTruth)))
GT=GT[0]


#reg1 = linear_model.LogisticRegressionCV(penalty='elasticnet',solver='saga',l1_ratios=LR1,max_iter=10000,Cs=C1,scoring='neg_mean_squared_error',multi_class='ovr').fit(AllDat,GT)
#reg2 = linear_model.LogisticRegressionCV(penalty='elasticnet',solver='saga',l1_ratios=LR2,max_iter=10000,Cs=C2,scoring='neg_mean_squared_error',multi_class='ovr').fit(AllDat,GT)
pickle.dump(reg1, open('ZPNPSI_MAFFTPromalAF2_Reg1_MultiClass_100C_20220424.sav', 'wb'))
pickle.dump(reg2, open('ZPNPSI_MAFFTPromalAF2_Reg2_MultiClass_100C_20220424.sav', 'wb'))


unreg = linear_model.LogisticRegressionCV(penalty='elasticnet',solver='saga',l1_ratios=[],max_iter=10000,Cs=C1,scoring='neg_mean_squared_error',multi_class='ovr').fit(AllDat,GT) #Put _ in multi_class
pickle.dump(unreg, open('ZPNPSI_MAFFTPromalAF2_UnReg_MultiClass_100C_20220424.sav', 'wb'))

reg = reg2

DatOne=reg.coefs_paths_[0.0]

DatX=reg.coefs_paths_[1.0]

DatM=reg.coefs_paths_[2.0]

ScoreOne=reg.scores_[0.0]
ScoreX=reg.scores_[1.0]
ScoreM=reg.scores_[2.0]

#print(reg.scores_[1.0]) #Key to the array is 1.0
AvgScoresOne=sum(reg.scores_[0.0])/len(reg.scores_[0.0])
StdScoresOne = np.std(reg.scores_[0.0])
AvgCoefsOne = sum(reg.coefs_paths_[0.0])/len(reg.coefs_paths_[0.0])

AvgScoresX=sum(reg.scores_[1.0])/len(reg.scores_[1.0])
StdScoresX = np.std(reg.scores_[1.0])
AvgCoefsX = sum(reg.coefs_paths_[1.0])/len(reg.coefs_paths_[1.0])

AvgScoresM=sum(reg.scores_[2.0])/len(reg.scores_[2.0])
StdScoresM = np.std(reg.scores_[2.0])
AvgCoefsM = sum(reg.coefs_paths_[2.0])/len(reg.coefs_paths_[2.0])


#print(AvgScores)
#fig = plt.figure()


BetaCutoff=0.0


Caxis=C1
Laxis=LR1


#Manual Reshape. Making all same size

CDat=[]
LDat=[]

ZOneDat=[]
ZOneStd=[]
POne=[]

ZXDat=[]
ZXStd=[]
PX=[]

ZMDat=[]
ZMStd=[]
PM=[]




for i in np.arange(len(AvgScoresOne)):
    #print(i)
    for j in np.arange(len(AvgScoresOne[i])):
        CDat.append(Caxis[i])
        LDat.append(Laxis[j])
        ZOneDat.append(AvgScoresOne[i][j])
        ZXDat.append(AvgScoresX[i][j])
        ZMDat.append(AvgScoresM[i][j])
        
        ParItOne = AvgCoefsOne[i][j]
        ParItX = AvgCoefsX[i][j]
        ParItM = AvgCoefsM[i][j]
        
        NuOne=0
        for p in ParItOne:
            if abs(p)>BetaCutoff:
                NuOne=NuOne+1
        POne.append(NuOne)
        
        NuX=0
        for p in ParItX:
            if abs(p)>BetaCutoff:
                NuX+=1
        PX.append(NuX)
        
        NuM=0
        for p in ParItM:
            if abs(p)>BetaCutoff:
                NuM+=1
        PM.append(NuM)
                
                
                
        
        ItVecOne=[]
        for k in np.arange(len(reg.scores_[0.0])):
            ItVecOne.append(reg.scores_[0.0][k][i][j])
            
        ZOneStd.append(np.std(ItVecOne))
        
        ItVecX=[]
        for k in np.arange(len(reg.scores_[1.0])):
            ItVecX.append(reg.scores_[1.0][k][i][j])
            
        ZXStd.append(np.std(ItVecX))
        
        ItVecM=[]
        for k in np.arange(len(reg.scores_[2.0])):
            ItVecM.append(reg.scores_[2.0][k][i][j])
            
        ZMStd.append(np.std(ItVecM))
        
        
        
        #print(j)

#print(np.log(Caxis))

#print(PDat)


MOne=max(ZOneDat)

BestInOne=0

NewCutOne=MOne

for i in np.arange(len(ZOneDat)):
    if ZOneDat[i] == MOne:
        #print("Matches")
        BestInOne=i
        NewCutOne = MOne - 1.96*ZOneStd[i]
        
        

MX=max(ZXDat)

BestInX=0

NewCutX=MX

for i in np.arange(len(ZXDat)):
    if ZXDat[i] == MX:
        #print("Matches")
        BestInX=i
        NewCutX = MX - 1.96*ZXStd[i]
        
        

MM=max(ZMDat)

BestInM=0

NewCutM=MM

for i in np.arange(len(ZMDat)):
    if ZMDat[i] == MM:
        #print("Matches")
        BestInM=i
        NewCutM = MM - 1.96*ZMStd[i]
        
 

ww_df = logomaker.get_example_matrix('ww_information_matrix',
                                     print_description=False)
       

#Look at params
#print(reg.scores_)
reg.coef_[0] #This is the main answer

#Make LogoDat for later

Val0One=[]
Pos0One=[]
ResN0One=[]



for i in np.arange(len(reg.coef_[0])):
    if i == len(reg.coef_[0])-1:
        continue   
    if abs(reg.coef_[0][i])>0:
        pos = math.floor(i/20) #Use 8 or 20
        AAid = i % 20 #Use 8 or 20. Or this is a PCA ID
        #AAnam = residues[AAid]
        print(reg.coef_[0][i],pos+1,residues[AAid])
        LetHeight=0
        RawBet = reg.coef_[0][i]
        
        if RawBet > 0:
            LetHeight=np.exp(RawBet)-1.0
        else:
            LetHeight=-1.0*(np.exp(-1.0*RawBet)-1.0)
        
        Val0One.append(LetHeight)
        Pos0One.append(pos+1)
        ResN0One.append(residues[AAid])
        #Last one is an intercept
        
print("GAP")
        
Val0X=[]
Pos0X=[]
ResN0X=[]

for i in np.arange(len(reg.coef_[1])):
    if i == len(reg.coef_[1])-1:
        continue   
    if abs(reg.coef_[1][i])>0:
        pos = math.floor(i/20) #Use 8 or 20
        AAid = i % 20 #Use 8 or 20. Or this is a PCA ID
        #AAnam = residues[AAid]
        print(reg.coef_[1][i],pos+1,residues[AAid])
        LetHeight=0
        RawBet = reg.coef_[1][i]
        
        if RawBet > 0:
            LetHeight=np.exp(RawBet)-1.0
        else:
            LetHeight=-1.0*(np.exp(-1.0*RawBet)-1.0)
        
        Val0X.append(LetHeight)
        Pos0X.append(pos+1)
        ResN0X.append(residues[AAid])
        #Last is an intercept
        
print("GAP")
        
Val0M=[]
Pos0M=[]
ResN0M=[]

for i in np.arange(len(reg.coef_[2])):
    if i == len(reg.coef_[2])-1:
        continue   
    if abs(reg.coef_[2][i])>0:
        pos = math.floor(i/20) #Use 8 or 20
        AAid = i % 20 #Use 8 or 20. Or this is a PCA ID
        #AAnam = residues[AAid]
        print(reg.coef_[2][i],pos+1,residues[AAid])
        LetHeight=0
        RawBet = reg.coef_[2][i]
        
        if RawBet > 0:
            LetHeight=np.exp(RawBet)-1.0
        else:
            LetHeight=-1.0*(np.exp(-1.0*RawBet)-1.0)
        
        Val0M.append(LetHeight)
        Pos0M.append(pos+1)
        ResN0M.append(residues[AAid])
        #Last is an intercept

LogoDat0One=[]
#loop over all positions
for a in np.arange(len(AllDat[0])/20):
    SiteDat0=[]
    #Loop over residues
    for b in residues:
        DatIt0=(b,0.0)
        for c in np.arange(len(Pos0One)):
            if Pos0One[c]==(a+1) and b==ResN0One[c]:
                DatIt0=(b,Val0One[c])
        SiteDat0.append(DatIt0)
    LogoDat0One.append(SiteDat0)
            
        
LogoDat0X=[]
#loop over all positions
for a in np.arange(len(AllDat[0])/20):
    SiteDat0=[]
    #Loop over residues
    for b in residues:
        DatIt0=(b,0.0)
        for c in np.arange(len(Pos0X)):
            if Pos0X[c]==(a+1) and b==ResN0X[c]:
                DatIt0=(b,Val0X[c])
        SiteDat0.append(DatIt0)
    LogoDat0X.append(SiteDat0)
            
LogoDat0M=[]
#loop over all positions
for a in np.arange(len(AllDat[0])/20):
    SiteDat0=[]
    #Loop over residues
    for b in residues:
        DatIt0=(b,0.0)
        for c in np.arange(len(Pos0M)):
            if Pos0M[c]==(a+1) and b==ResN0M[c]:
                DatIt0=(b,Val0M[c])
        SiteDat0.append(DatIt0)
    LogoDat0M.append(SiteDat0)
    
def MakeLogoDF (LogoDat):
    Logo_df = pandas.DataFrame(0,columns=ww_df.columns,index=np.arange(len(LogoDat)))
    for i in np.arange(len(LogoDat)):
        #Make Row
        NuRow=[]
        for j in np.arange(len(LogoDat[i])):
            NuRow.append(LogoDat[i][j][1])
            Logo_df.iloc[i,j] = LogoDat[i][j][1]
    return(Logo_df)


Logo_dfOne = MakeLogoDF(LogoDat0One)

Logo_dfX = MakeLogoDF(LogoDat0X)

Logo_dfM = MakeLogoDF(LogoDat0M)

AllAvgScores=[AvgScoresOne,AvgScoresX,AvgScoresM]
AllAvgCoefs=[AvgCoefsOne,AvgCoefsX,AvgCoefsM]
AllNewCut=[NewCutOne,NewCutX,NewCutM]

CinVec=[]
LinVec=[]

for i in np.arange(len(AllAvgScores)):
    AvgScores=AllAvgScores[i]
    AvgCoefs=AllAvgCoefs[i]
    NewCut=AllNewCut[i]
    #Xq=np.quantile(ZDat,0.95)
    #print(Xq)
    #NewCut Defined above
    #Allow other cutoff
    #Loop over new cutoff
    BestCin=0
    BestLin=0
    MinNP = len(reg.coefs_paths_[0][0][0][0])
    ModScore=-100*MinNP
    #NuPa=MinNP


    Cnu=[]
    LRnu=[]
    Scores=[]
    ParamV=[]


    SkipCon=0

    for a in np.arange(len(AvgScores)):
            for b in np.arange(len(AvgScores[a])):
                Take = AvgCoefs[a][b]
                SkipCon+=1
                #if SkipCon>10 : break
                #print(AvgScores[a][b])
                #cross_validate.cross_val_predict(reg,AlVal,GTV)
                #cross_validate(reg,AlVal,GTV)
                #OneOff=linear_model.LogisticRegression().fit()
                OneOff=reg
                OneOff.intercept=Take[len(Take)-1]
                OneOff.coef=Take[:-1]
                OneScore=OneOff.score(ValDat,GTV)


                #print(OneOff)
                #print(len(AvgCoefs[a][b]))

                if OneScore>= NewCut:

                    Take = AvgCoefs[a][b]
                    #print(len(Take),Take[0],Take[len(Take)-1])
                    NuPa = 0

                    #print(MinNP)
                    for t in Take:
                        if abs(t) > 0:
                            NuPa+=1

                    if NuPa < MinNP:
                        MinNP = NuPa
                        BestCin=a
                        BestLin=b
                        ModScore=OneScore

                    elif NuPa == MinNP:
                        if OneScore> ModScore:
                            #MinNP = NuPa
                            BestCin=a
                            BestLin=b
                            #Tie breaker if two have same number of parameters
                            ModScore=OneScore

                #print("accept",AvgScores[a][b],C1[a],LR1[b],NuPa,MinNP)
                #print(len(Take),Take[0],Take[len(Take)-1])
    print("Final accept",ModScore,C1[BestCin],LR1[BestLin],MinNP)
    CinVec.append(BestCin)
    LinVec.append(BestLin)

#Save selected model, but do it in a looped way

#AllLogoDat=['A','B','C']
AllLogoDat={}

for w in np.arange(len(AllAvgScores)):
    PosX=[]
    BetaY=[]
    ResI=[]
    ResN=[]

    #Sweep=reg3.coef_[0]

    AvgCoefs=sum(reg.coefs_paths_[w])/len(reg.coefs_paths_[w])

    Sweep = AvgCoefs[CinVec[w]][LinVec[w]]
    #print(len(AllDat[0])/20)


    for i in np.arange(len(Sweep)):
        if abs(Sweep[i])>0:
            pos = math.floor(i/20) #Use 8 or 20
            if pos == len(AllDat[0])/20:
                break
                #Let's me skip the intercept term
            AAid = i % 20 #Use 8 or 20. Or this is a PCA ID
            #AAnam = residues[AAid]
            #print(Sweep[i],pos+1,residues[AAid])
            PosX.append((pos+1))
            RawBet=Sweep[i]
            LetHeight=0

            if RawBet > 0:
                LetHeight=np.exp(RawBet)-1.0
            else:
                LetHeight=-1.0*(np.exp(-1.0*RawBet)-1.0)

            #Undo Exponential Transform:
            #LetHeight=RawBet


            BetaY.append(LetHeight)
            ResI.append(AAid)
            ResN.append(residues[AAid])


    LogoDat=[]
    #loop over all positions
    for a in np.arange(len(AllDat[0])/20):
        SiteDat=[]
        #Loop over residues
        for b in residues:
            DatIt=(b,0.0)
            for c in np.arange(len(PosX)):
                if PosX[c]==(a+1) and b==ResN[c]:
                    DatIt=(b,BetaY[c])
            SiteDat.append(DatIt)
        LogoDat.append(SiteDat)
#print(AllLogoDat)
    AllLogoDat[w] = LogoDat # Must be in the loop!
#print(AllLogoDat)

#print(AllLogoDat)



LogOne= MakeLogoDF(AllLogoDat[0])
LogX= MakeLogoDF(AllLogoDat[1])
LogM= MakeLogoDF(AllLogoDat[2])

#Convert values #Start with ZP Module


#Wrong Alignment Examples



MouseZP3_UnAlig = 'VKVECLEAELVVTVSRDLFGTGKLVQPGDLTLGSEGCQPRVSVDTDVVRFNAQLHECSSRVQMTKDALVYSTFLLHDPRPVSGLSILRTNRVEVPIECRYPR'
#RatZP3_Alig = '-VE-VECKE-----A---E---------------LVVT----VRRDL-FGTG-K-L----------------------------------------VQP-------GD-L---TL-G--SE--------G-CQ--------------P-LVAV----D-T--------DVVRLNAQLH-E----------------CSSG-------VQVT----------------------E-D----A-L-VYS--TF-LLHD-PRPV------------------------N-GLSILRTN--RV-E--V------PIE-CRY-P-R-------'
#RatZP3_Alig = '--VKVECLEA-ELVVT--VS---R------D--------------------L-F---GT----------------------------GK-LV--------------Q-PG-D---LT-L--G----SEG-C-------------------------Q--PRVS-VD---------T--D-VV--RFN-----A-Q-LH---------ECSSR--V--Q--M---T--KD-------------------AL---VYST-FL--LH-D----P------RPVS-----------G-----L---SI-LRTNR--VEVP-----IECRY-PR'
#This is the renewed alignment

RatZP3_Alig = '-VKVECLEA-E--L-------VV---------T---------V------S-RDLFG----T-------------------------------------GK-----------LV----Q--P---G------------------D---LT-L------------GSE-------------------------------------------------------------------------G-CQP--R-----V-S----VD--T--D---V-------------VRFNAQL------------HEC-------SSRV--Q---M---------T-----------------------K-D----AL-----VY-STF----LL---H---------------D--PR----P-V--S--G-LS------I------L---RTNR-----VEVP-----IECR-YPR'






MouseZP2ZPN1_UnAlig='GTLICDKDEVRIEFSSRFDMEKWNPSVVDTLGSEILNCTYALDLERFVLKFPYETCTIKVVGGYQVNIRVGDTTTDVRYKDDMYHFFCPAIQ'



#RatZP2ZPN_Alig='-GT-LICDK-----D---E---------------VRVE----FSSRF-D----M-E-KW-----------------------------------N-PSL--VD-TFGN-E---IS--------------------------------N-CTYA----L-DLEK-------FILKFPYE-T----------------C----------TIKV----------------------I-G----G-Y-QVN--IR-V--Q-DTNA------D-----------------V------SYK--DD-V--H------HFF-CPA-I-Q-------'
#RatZP2ZPN_Alig = '--GTLICDKD-EVRIE--FS---S------R--------------------F-D---M------------------------------------------------E-KW-N---PS-V--V----D-----------------------------T--LGSEILNC-TYAL--DL--E-RF--VLK-----F-P-YE---------TCTI------K--V---V--G-------------------------GYQV-NI--RV-G----D------TTT-----------------------D-VRYKD--DMYH-----FFCPA-IQ'
#New PSI-BLAST based alignment 

RatZP2ZPN_Alig = '-GTLICDKD-E--V-------RI---------E---------F------S-SRF---------------------------------------------------------DM----E--K---W------------------N---PS-V------------VDT-------------------------L---------G----SE-----------------------I------LN-CTY--A-----L-D-----------L-E-R-------------FVLKFPY------------ETC-------T--I--K---V---------V-------------------------------G-----GY-QVN----IR---V---------------G--DT----T----------T------D------V---RYKD-----DMYH-----FFCP-AIQ'

#This is really a mouse.

#Identify gap only sites in both

BothGap=[]
FineSite=[]

for i in np.arange(len(RatZP3_Alig)):
    if RatZP3_Alig[i]=='-' and RatZP2ZPN_Alig[i] == '-' :
        BothGap.append(i)
    else:
        #print(sum(Logo_df.iloc[i,])) #Not I lose one of my parameters because it falls in a gap
        FineSite.append(i)

        
CleanLogoOne=LogOne.iloc[FineSite,:].reset_index(drop=True)
CleanLogoX=LogX.iloc[FineSite,:].reset_index(drop=True)
CleanLogoM = LogM.iloc[FineSite,:].reset_index(drop=True)

CleanLogoOne.to_csv('CleanLogoOne_100C_20220424.csv',index = False)
CleanLogoX.to_csv('CleanLogoX_100C_20220424.csv',index = False)
CleanLogoM.to_csv('CleanLogoM_100C_20220424.csv',index = False)








    
# %%
