#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
import operator
from operator import itemgetter
import logomaker
import pickle


# In[2]:


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
    


# In[3]:


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

   


# In[4]:


def SimpMult(A,B):
    X=0.0
    for i in np.arange(len(A)):
        X+=A[i]*B[i]
        
    return(X)
    
def LogitEst(V1,V2,c):
    polynom = SimpMult(V1,V2)+c
    Ans=1/(1+np.exp(-1*polynom))
    return(Ans)


# In[5]:


residues = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
           'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

mapper = dict(zip(residues, [[1.0 if residues[i] == r else 0.0 for i in range(len(residues))] for r in residues]))
mapper['-'] = [0.0]*20


# In[6]:


def ConvertRes(seq,mapperx):
    #print(seq)
    results = []
    for i in range(len(seq)):
        results += mapperx[seq[i]]
    return(results)


# In[7]:


os.chdir(r'C:\Users\alber\OneDrive\Documents\Swanson_2017Onward\Sequences\ZP_Alignments')


LastZPN = ["ZP1ZPN2","ZP2ZPN4","ZP3ZPN1","ZP4ZPN2","ZPAXLastZPN","UMODZPN","TECTAZPN","CUZD1ZPN","ZPDZPN"]


#ValDat = ReadIn('AlgoAln_ZPN_PSI90c_Test_202111031718.fasta')
#AlVal = ReadIn('AlgoAln_ZPN_90c_20200625_25Test_202009241307.fasta')
#AlVal = ReadIn('AlgoAln_ZPN_90c_20200625_25Test_202011021352.fasta')
#AlVal = ReadIn('AlgoAln_ZPN_PSI90c_Test_202111031718.fasta')
AlVal = ReadIn('AlgoAln_ZPN_PSI100c_MAFFTProMals_Test_202203242127.fasta')


ValDat=[]
NL=[]

for j in np.arange(len(AlVal)):
    YLine = ConvertRes(AlVal[j].seq,mapper) #Decide which mapper you want
    ValDat.append(YLine)
    
    NLi = AlVal[j].name.split('_')
    NL.append(NLi[0]+NLi[1])
    
    
    
GTV = [1.0 if NL[i] in LastZPN else 0.0 for i in range(len(NL))]
#GTV
#print(np.ndarray(ValDat))


# In[8]:


20*len(AlVal[0].seq)


# In[9]:


#F1 = open('ZPNPSI_Reg1_NoZoomLocal_202112021554.sav','rb')
#F2 = open('ZPNPSI_Reg2_ZoomLocal_202112021554.sav','rb')

#R1 = pickle.load(open('ZPNPSI_Reg1_NoZoom_202112021554.sav','rb'))
#R2 = pickle.load( open('ZPNPSI_Reg2_Zoom_202112021554.sav','rb'))
#R3 = pickle.load(open('ZPNPSI_Reg3_Zoom_202112061916.sav','rb'))


R1 = pickle.load(open('ZPNPSI100c_MAFFTPromalAF2_Reg1_202203242152.sav','rb'))

R2 = pickle.load(open('ZPNPSI100c_MAFFTPromalAF2_Reg2_202203242152.sav','rb'))

R3 = pickle.load(open('ZPNPSI100c_MAFFTPromalAF2_Reg3_202203242152.sav','rb'))



#F1.close()
#F2.close()



# C:\Users\alber\Anaconda3\envs\myenv8\lib\site-packages\sklearn\utils\deprecation.py:143: FutureWarning: The sklearn.linear_model.logistic module is  deprecated in version 0.22 and will be removed in version 0.24. The corresponding classes / functions should instead be imported from sklearn.linear_model. Anything that cannot be imported from sklearn.linear_model is now part of the private API.
#   warnings.warn(message, FutureWarning)
# C:\Users\alber\Anaconda3\envs\myenv8\lib\site-packages\sklearn\base.py:329: UserWarning: Trying to unpickle estimator LogisticRegressionCV from version 0.21.2 when using version 0.23.2. This might lead to breaking code or invalid results. Use at your own risk.
#   warnings.warn(

# In[10]:



#F1.close()
#F2.close()
(R3.scores_[1.0].shape)


# In[11]:


#reg = R1

reg = R2

#reg = R3

AvgScores=sum(reg.scores_[1.0])/len(reg.scores_[1.0])

StdScores = np.std(reg.scores_[1.0])



AvgCoefs = sum(reg.coefs_paths_[1.0])/len(reg.coefs_paths_[1.0])



# In[12]:


#fig = plt.figure()




LR1=np.arange(0,1,0.02)
LR2=np.arange(0.4,0.8,0.02) #ZOomed in

LR3 = np.arange(0.25,0.75,0.01)
#LRs = [0.01*i for i in np.arange(101)]
#Alphs = [0,0.001,0.01,0.1,1]
Alphs=[1]

C1=[1,0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001] # Skip first 4 C's In the 0.01 to the 0.05 range
C2=sorted(np.arange(0.01,0.05,0.002),reverse=True)

C3 = sorted(np.arange(0.01,0.1,0.002),reverse=True)

#Fix when using R2 or R3
C1=C2
LR1=LR2


#C1 = C1
#LR1=LR1


Caxis=C1
Laxis=LR1

#Caxis = C2
#Laxis = LR2


#Manual Reshape. Making all same size

CDat=[]
LDat=[]
ZDat=[]
ZStd=[]
PDat=[]




for i in np.arange(len(AvgScores)):
    #print(i)
    for j in np.arange(len(AvgScores[i])):
        CDat.append(Caxis[i])
        LDat.append(Laxis[j])
        ZDat.append(AvgScores[i][j])
        
        ParIt = AvgCoefs[i][j]
        
        NuParIt=0
        for p in ParIt:
            if abs(p)>0:
                NuParIt=NuParIt+1
        PDat.append(NuParIt)
                
        
        ItVec=[]
        for k in np.arange(len(reg.scores_[1.0])):
            ItVec.append(reg.scores_[1.0][k][i][j])
            
        ZStd.append(np.std(ItVec))
        
        #print(j)

#print(np.log(Caxis))

#print(PDat)


# In[13]:


M=max(ZDat)

BestIn=0

NewCut=M

for i in np.arange(len(ZDat)):
    if ZDat[i] == M:
        #print("Matches")
        BestIn=i
        NewCut = M - 1.96*ZStd[i]
        
#print(M,ZStd[BestIn],NewCut)
if NewCut==0:
    NewCut = np.quantile(ZDat,0.95)

    
print(NewCut)
    
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})


Lam=[1/i for i in C1]

Cgood, LRgood = np.meshgrid(Caxis,Laxis)






#Lamgood, LRgood = np.meshgrid(Lam,LR1)
#print(Cgood)

#ax.scatter3D(CDat, LDat, AvgScores,c=AvgScores<-0.4)
#ax.scatter3D(LDat,np.log(CDat),ZDat)
#ax.scatter3D(LDat,CDat,ZDat,c=AvgScores<=NewCut,cmap=cm.seismic) #Linear for this one #THis works!!!
#C~0.026, LR ~ 0.78

#ax.plot_surface(Cgood, LRgood,AvgScores)
#When I input an array, I have to keep track of the correct side based on the 


#ax.plot_surface(LRgood,Cgood,AvgScores,cmap=cm.coolwarm) 
#ax.plot_surface(LRgood,Lamgood,AvgScores,cmap=cm.coolwarm) 

#ax.plot_surface(LRgood,Cgood,AvgScores)


#ax.plot_surface(Cgood,LRgood,AvgScores,cmap=cm.coolwarm)

print(Cgood.shape,LRgood.shape,AvgScores.shape)

ax.plot_surface(Cgood,LRgood,np.transpose(AvgScores),cmap=cm.coolwarm)

LamDat = [1/x for x in CDat]

ax.plot3D(CDat,LDat,ZDat)

#ax.plot_surface(np.transpose(Cgood),np.transpose(LRgood),AvgScores,cmap=cm.coolwarm)

#plot_surface(C1,LR1,AvgScores)
plt.show()


# In[14]:


AvgScores.shape


# In[15]:


len(C1)


# In[16]:


len(LR1)


# In[17]:


len(Lam)


# In[18]:


Cgood.shape


# In[19]:


LRgood.shape


# In[20]:


#Xq=np.quantile(ZDat,0.95)
#print(Xq)



#NewCut Defined above

#Allow other cutoff

#Loop over new cutoff
BestCin=0
BestLin=0
MinNP = len(reg.coefs_paths_[1.0][0][0][0])
ModScore=-100*MinNP
#NuPa=MinNP


Cnu=[]
LRnu=[]
Scores=[]
ParamV=[]


SkipCon=0

OptAcc=0

for a in np.arange(len(AvgScores)):
        for b in np.arange(len(AvgScores[a])):
            Take = AvgCoefs[a][b]
            SkipCon+=1
            #if SkipCon>10 : break
            #print(AvgScores[a][b])
            #cross_validate.cross_val_predict(reg,AlVal,GTV)
            #cross_validate(reg,AlVal,GTV)
            #OneOff=linear_model.LogisticRegression().fit()
            OneOff=reg #This maintains settings from before
            OneOff.intercept=Take[len(Take)-1]
            OneOff.coef=Take[:-1]
            OneScore=OneOff.score(ValDat,GTV)
            AccMod=OneOff
            AccMod.scoring='accuracy'
            AccScore = AccMod.score(ValDat,GTV)
            #Get accuracy
            
            
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
                    OptAcc=AccScore
                    
                elif NuPa == MinNP:
                    if OneScore> ModScore:
                        #MinNP = NuPa
                        BestCin=a
                        BestLin=b
                        #Tie breaker if two have same number of parameters
                        ModScore=OneScore
                        OptAcc=AccScore
                    
            print("accept",AvgScores[a][b],C1[a],LR1[b],NuPa,MinNP)
            #print(len(Take),Take[0],Take[len(Take)-1])
print("Final accept",ModScore,C1[BestCin],LR1[BestLin],MinNP)
print("Optimized Accuracy",OptAcc)
      
#12/15/21 100% Accuracy


# In[21]:


print(C1,BestCin,LR1,BestLin)


# In[22]:


len(LR1)


# In[23]:


Lamb1=[]
ScoreC=[]
StdC=[]

LRVec=[]
ScoreL=[]
StdL=[]

NumParC=[]

NumParL=[]


for i in np.arange(len(ZDat)):
    if LDat[i] == LR1[BestLin]:
        Lamb1.append(1/CDat[i])
        ScoreC.append(ZDat[i])
        StdC.append(ZStd[i])
        
        NumParC.append(PDat[i])
        
        
    if CDat[i] == C1[BestCin]:
        LRVec.append(LDat[i])
        ScoreL.append(ZDat[i])
        StdL.append(ZStd[i])
        
        NumParL.append(PDat[i])


plt.scatter(Lamb1,ScoreC)


# In[24]:





plt.errorbar(Lamb1,ScoreC,c="black",yerr=StdC) #Changed to black


for i in np.arange(len(NumParC)):
    #if i in [1,2,3,4,5,6,7,8,9]: continue
    
    plt.text(Lamb1[i],0.2,str(NumParC[i]),fontsize=12)


plt.show()


# In[146]:


fig0, ax1 = plt.subplots(2,gridspec_kw={'height_ratios':[1,2]}) #Re-do with different plots

ax1[0].plot(Lamb1,NumParC,c='black')
ax1[0].set_ylabel("Parameter Number")
ax1[0].set_title('Degree of Regularization')

ax1[0].axes.xaxis.set_visible(False)
#Putting parameters on top

MinLamb1 = min(Lamb1)
MaxLamb1 = max(Lamb1)

ax1[0].set_xlim(MinLamb1,MaxLamb1)
ax1[1].set_xlim(MinLamb1,MaxLamb1)

#ax2=ax1.twinx() #Shares same x-axis
#ax1[0].set_xlabel("Strength of Regularization")

#ax1.set_ylim(-1.4,0.1)
#ax2.set_ylim(14,19.5)
#ax2.set_ylabel("Number of Parameters",color="blue")
#ax2.tick_params(axis='y',labelcolor="blue")
ax1[1].set_ylabel("Score")


#ax1[1].errorbar(Lamb1,ScoreC,c="black",yerr=StdC) #pre 1/2022 this was error bars



LowC = [ScoreC[i]-StdC[i] for i in np.arange(len(StdC))]
HighC = [ScoreC[i]+StdC[i] for i in np.arange(len(StdC))]


ax1[1].plot(Lamb1,ScoreC,c="black",linewidth=5) #1/25/2022

ax1[1].fill_between(Lamb1,LowC,HighC,facecolor = '#dddddd')

#ax1[1].set_xlabel("Degree of Regularization")

ax1[1].hlines(y=(ScoreC[0]-StdC[0]),xmin=MinLamb1,xmax=100,linewidth=4, color='b', linestyle = 'dotted')

#ax1[0].axvspan(38.46,40, alpha=0.5, color='#222222')
#ax1[1].axvspan(38.46,40, alpha=0.5, color='#222222')

ax1[0].axvspan(24.5,25.5, alpha=0.5, color='green')
ax1[1].axvspan(24.5,25.5, alpha=0.5, color='green')


ax1[0].set_xlim(21,30)
ax1[1].set_xlim(21,30)

ax1[0].set_ylim(20,45)
ax1[1].set_ylim(-0.03,0)


fig0.set_size_inches(4,4)

fig0.tight_layout()
#fig0.savefig('OneStandardErrorRule_100C_202204061618.png',dpi=600)


# In[26]:


#Get location on graph

cutwant = LowC[0]
#38.46 to 41.67

for i in np.arange(len(ScoreC)):
    if ScoreC[i]>cutwant:
        print(Lamb1[i],'good')
    else:
        print(Lamb1[i],'bad')
       


# In[27]:


1/C1[BestCin]


# In[28]:


Lam


# In[29]:


plt.errorbar(LRVec,ScoreL,c='red',yerr=StdL)


# In[145]:


LRCons = [LR1[BestLin] for x in np.arange(len(Lamb1))]
fig,ax = plt.subplots(subplot_kw={"projection": "3d"})

#ax.zaxis.set_label_position("bottom")
#ax.zaxis.tick_bottom()


#ax.invert_zaxis()
#axo=ax.twinx()

#ax=plt.axes(projection='3d')


#Slice and Surface placement 


#Cgood, LRgood = np.meshgrid(C1,LR1)
Lamgood, LRgood = np.meshgrid(Lam,LR1)
#print(Cgood)

#ax.scatter3D(CDat, LDat, AvgScores,c=AvgScores<-0.4)
#ax.scatter3D(LDat,np.log(CDat),ZDat)
#ax.scatter3D(LDat,CDat,ZDat,c=AvgScores<=NewCut,cmap=cm.seismic) #Linear for this one
#C~0.026, LR ~ 0.78

#ax.plot_surface(Cgood, LRgood,AvgScores)
#When I input an array, I have to keep track of the correct side based on the 


#ax.plot_surface(LRgood,Cgood,AvgScores,cmap=cm.coolwarm) 

ax.plot_surface(LRgood,Lamgood,AvgScores.T,color='lightgray') #May need to transpose

#ax.scatter3D(LRCons,Lamb1,ScoreC,c=ScoreC,cmap=cm.coolwarm) #This does color gradient
ax.scatter3D(LRCons,Lamb1,ScoreC,color='black')
#LRCons is the best L1-Ratio

YPlane=np.arange(21,30)
ZPlane=np.arange(-0.03,0,0.005)

#YPlane=np.arange(0,120)
#ZPlane=np.arange(-1,0.2,0.01)

yy,zz = np.meshgrid(YPlane,ZPlane)
X1=[LR1[BestLin] for x in np.arange(len(YPlane))]
X2=[LR1[BestLin] for x in np.arange(len(ZPlane))]

xx1, xx2= np.meshgrid(X1,X2)
#ax.plot_surface(xx1,yy,zz,color=[1.0,0.0,0.0,0.5])



xx, yy = np.meshgrid(LRCons,YPlane)

normal =[0,1,0]

normal=[]
point=[LR1[BestLin] , Lamb1[0],ScoreC[0]]

print(point)

zz= xx

#ax.plot_surface(xx,yy,zz,color='red')


ax.set_xlabel("Regularization Ratio")



ax.set_ylabel("Degree of Regularization")
#ax.set_zlabel("Score")


#ax.set_zticks([])
ax.axes.zaxis.set_visible('False')




#ax.set_ylim(21,30)

#ax.set_xlim(20,45)
#ax.set_zlim(-0.03,0)



#axz=ax.twinx()
#axz.set_zlabel("Doof")

#ax.plot_surface(LRgood,Cgood,AvgScores)
fig.set_size_inches(6,5)
fig.tight_layout()


#ax.set_title("Model Scores",fontsize=12)

#ax.plot_surface(Cgood,LRgood,AvgScores,cmap=cm.coolwarm)
#fig.savefig('ModelSelection_Waterfall_202204061615_NoSq.png',dpi=600)
plt.show()


# In[31]:


print(LRgood.shape,Lamgood.shape,AvgScores.shape)


# In[32]:


#Lamgood


# In[33]:


plt.errorbar(LDat,ZDat,yerr=ZStd,c='red')
plt.show()


# In[34]:


fig0, ax1 = plt.subplots(2,gridspec_kw={'height_ratios':[1,2]}) #Re-do with different plots



ax1[0].plot(LRVec,NumParL,c='gray')
ax1[0].set_ylabel("Parameter Number")

ax1[0].axes.xaxis.set_visible(False)
#Putting parameters on top


#ax1[0].set_xlim(20,40)
#ax1[1].set_xlim(20,40)

#ax2=ax1.twinx() #Shares same x-axis
#ax1[0].set_xlabel("Strength of Regularization")

#ax1.set_ylim(-1.4,0.1)
#ax2.set_ylim(14,19.5)
#ax2.set_ylabel("Number of Parameters",color="blue")
#ax2.tick_params(axis='y',labelcolor="blue")
ax1[1].set_ylabel("Score")
ax1[1].errorbar(LR1,ScoreL,c="black",yerr=StdL) 
ax1[1].set_xlabel("Regulatization Type")

ax1[1].hlines(y=(ScoreL[0]-StdL[0]),xmin=0.4,xmax=0.8,linewidth=2, color='r') #Change x limits #Make this line dotted



fig0.set_size_inches(4,4)

fig0.tight_layout()
#fig0.savefig('OneStandardErrorRule_202106162338.png',dpi=600)


# In[35]:


#Make Penalty1, Pentalty2

P1 = []
P2 = []


for a in np.arange(len(AvgScores)):
    
    for b in np.arange(len(AvgScores[a])):
        
        cost1 = LR1[b] * (1 / C1[a])
        cost2 = (1.0-LR1[b]) * (1 / C1[a])
        
        P1.append(cost1)
        P2.append(cost2)


# In[36]:


P1use = np.reshape(np.array(P1),(len(C1),len(LR1)))
P2use = np.reshape(np.array(P2),(len(C1),len(LR1)))
Zuse = np.reshape(np.array(ZDat),(len(C1),len(LR1)))




# In[37]:


fig,ax = plt.subplots(subplot_kw={"projection": "3d"})

ax.plot_surface(P1use,P2use,Zuse)


# In[102]:


###BREAK##



PosX=[]
BetaY=[]
ResI=[]
ResN=[]

#Sweep=reg3.coef_[0]

AvgCoefs=sum(reg.coefs_paths_[1.0])/len(reg.coefs_paths_[1.0])

Sweep = AvgCoefs[BestCin][BestLin]
#print(len(ValDat[0])/20)


for i in np.arange(len(Sweep)):
    if abs(Sweep[i])>0:
        pos = math.floor(i/20) #Use 8 or 20
        if pos == len(ValDat[0])/20:
            break
            #Let's me skip the intercept term
        AAid = i % 20 #Use 8 or 20. Or this is a PCA ID
        #AAnam = residues[AAid]
        print(Sweep[i],pos+1,residues[AAid])
        PosX.append((pos+1))
        RawBet=Sweep[i]
        LetHeight=0
        
        if RawBet > 0:
            LetHeight=np.exp(RawBet)-1.0
        else:
            LetHeight=-1.0*(np.exp(-1.0*RawBet)-1.0)
            
        #Undo Exponential Transform:
        LetHeight=RawBet
        
        
        BetaY.append(LetHeight)
        ResI.append(AAid)
        ResN.append(residues[AAid])


# In[103]:


def GetLogoVecs(Sweep):
    
    ###BREAK##
    #print(Sweep)


    PosX=[]
    BetaY=[]
    ResI=[]
    ResN=[]

    for i in np.arange(len(Sweep)):
        if abs(Sweep[i])>0:
            pos = math.floor(i/20) #Use 8 or 20
            if pos == len(ValDat[0])/20:
                break
                #Let's me skip the intercept term
            AAid = i % 20 #Use 8 or 20. Or this is a PCA ID
            #AAnam = residues[AAid]
            print(Sweep[i],pos+1,residues[AAid])
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
    return PosX, BetaY, ResI, ResN


# In[104]:


PosXU, BetaYU, ResIU, ResNU = GetLogoVecs(R1.coef_[0])


# In[105]:


LogoDat=[]
#loop over all positions
for a in np.arange(len(ValDat[0])/20):
    SiteDat=[]
    #Loop over residues
    for b in residues:
        DatIt=(b,0.0)
        for c in np.arange(len(PosX)):
            if PosX[c]==(a+1) and b==ResN[c]:
                DatIt=(b,BetaY[c])
        SiteDat.append(DatIt)
    LogoDat.append(SiteDat)
            


# In[106]:


logomaker.get_example_matrix('ww_information_matrix',print_description=False)    


# In[107]:


def GetLogoPanda(PosX,BetaY,ResI,ResN):
    LogoDat=[]
#loop over all positions
    for a in np.arange(len(ValDat[0])/20):
        SiteDat=[]
        #Loop over residues
        for b in residues:
            DatIt=(b,0.0)
            for c in np.arange(len(PosX)):
                if PosX[c]==(a+1) and b==ResN[c]:
                    DatIt=(b,BetaY[c])
            SiteDat.append(DatIt)
        LogoDat.append(SiteDat)
        
        
    ww_df = logomaker.get_example_matrix('ww_information_matrix',print_description=False)    
    Logo_df = pandas.DataFrame(0,columns=ww_df.columns,index=np.arange(len(LogoDat)))


    for i in np.arange(len(LogoDat)):
        #Make Row
        NuRow=[]
        for j in np.arange(len(LogoDat[i])):
            NuRow.append(LogoDat[i][j][1])
            Logo_df.iloc[i,j] = LogoDat[i][j][1]
    return(Logo_df)


# In[108]:


UnReg = GetLogoPanda(PosXU, BetaYU, ResIU, ResNU)


# In[109]:


# load ww information matrix
fig, ax = plt.subplots()
fig.set_size_inches(7,2)
ww_df = logomaker.get_example_matrix('ww_information_matrix',
                                     print_description=False)

# create Logo object
#logomaker.Logo(ww_df)


#print(ww_df.columns)
print(ww_df.index)

Logo_df = pandas.DataFrame(0,columns=ww_df.columns,index=np.arange(len(LogoDat)))


for i in np.arange(len(LogoDat)):
    #Make Row
    NuRow=[]
    for j in np.arange(len(LogoDat[i])):
        NuRow.append(LogoDat[i][j][1])
        Logo_df.iloc[i,j] = LogoDat[i][j][1]
    #print(LogoDat[i])
    #print(NuRow)
    #Row_df=pandas.DataFrame(np.transpose(NuRow),columns=residues)
    #Logo_df.append(NuRow,ignore_index=True)
    #if i >10:
       # break
        

#print(ww_df.head())
#print(Logo_df)
logomaker.Logo(Logo_df,color_scheme='hydrophobicity',ax=ax)
#fig.savefig('OneRow_PsuedoLogo_202104141508.png',dpi=600)


# In[110]:


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

#Renewed alignment


# In[111]:


##Identify gap only sites in both

BothGap=[]
FineSite=[]

for i in np.arange(len(RatZP3_Alig)):
    if RatZP3_Alig[i]=='-' and RatZP2ZPN_Alig[i] == '-' :
        BothGap.append(i)
    else:
        #print(sum(Logo_df.iloc[i,])) #Not I lose one of my parameters because it falls in a gap
        FineSite.append(i)

        
CleanLogo=Logo_df.iloc[FineSite,:].reset_index(drop=True) 
CleanUnReg=UnReg.iloc[FineSite,:].reset_index(drop=True) 
CleanLogo


# In[112]:


fig3, ax3 = plt.subplots(2)
fig3.set_size_inches(7,2)



#logomaker.Logo(CleanUnReg, ax =ax3[0], color_scheme='skylign_protein')
#logomaker.Logo(CleanLogo, ax =ax3[1], color_scheme='skylign_protein')
logomaker.Logo(CleanLogo, ax =ax3[1], color_scheme='black')

ax3[0].set_ylim(-1.5,1.5)
ax3[1].set_ylim(-1.5,1.5)
ax3[0].set_xticks([])
fig3.tight_layout()

fig3.subplots_adjust(hspace=0.05)
#fig3.savefig('Clean_PsuegoLogo_202104141738.png',dpi=600)


# In[113]:


fig3, ax3 = plt.subplots(2)
fig3.set_size_inches(7,2)



logomaker.Logo(CleanUnReg, ax =ax3[0], color_scheme='skylign_protein')
logomaker.Logo(CleanLogo, ax =ax3[1], color_scheme='skylign_protein')
#logomaker.Logo(CleanLogo, ax =ax3[1], color_scheme='black')


fig3.tight_layout()

fig3.subplots_adjust(hspace=0.05)
#fig3.savefig('TestLogo_202112082210.png',dpi=600)


# In[114]:


R1


# In[115]:


#Turn R1 into good LogoDat

print(len(R1.coef_[0]),len(AlVal[0].seq)*20 )


# In[ ]:





# In[116]:


#CleanUnReg.to_csv('CleanUnReg_LogoDat_202112131147x.csv')
#CleanLogo.to_csv('CleanReg_LogoDatR?_202112131147.csv')


# In[117]:


#Add Colors


ModSeq=''
FreeSeq=''



for i in np.arange(len(CleanLogo)):
    Take=CleanLogo.iloc[i,]
    if (sum(Take!=0.0))==0.0:
        #Remember that the sum of Trues (non-zero values) is 0
        ModSeq+='-'
        FreeSeq+='-'
    else:
        #print(Take.index(min(Take)))
        VminI, Vmin = min(enumerate(Take), key=operator.itemgetter(1)) #With enumerate I output both the index and value
        
        VmaxI, Vmax = max(enumerate(Take), key=operator.itemgetter(1))
        #print(Take)
        #print(Vmin,Vmax)
        if Vmin <0:
            FreeSeq+=residues[VminI] #Testing
            #print(VminI)
        else:
            FreeSeq+='-'
        if Vmax > 0:
            ModSeq+=residues[VmaxI]
            #print(VmaxI)
        else:
            ModSeq+='-'
    
print(ModSeq)
print(FreeSeq)





# In[118]:


aa_pi = {'G':5.97, 'A':6.00, 'V':5.96, 'L':5.98,
        'I':6.02, 'M':5.74, 'P':6.30, 'F':5.48,
        'W':5.89, 'N':5.41, 'Q':5.65, 'S':5.68,
        'T':5.60, 'Y':5.66, 'C':5.07, 'D':2.77,
        'E':3.22, 'K':9.74, 'R':10.76, 'H':7.59}

aa_sorted = sorted(aa_pi, key=lambda x:aa_pi[x], reverse=True)
colors = cm.turbo(np.linspace(0.1,0.9,20))
#print(colors)
print([x[:-1] for x in colors]) 
#colors=[x[:-1] for x in colors] #Remove last entry which is always a 1?
aa_colors = dict(zip(aa_sorted,colors))


# In[119]:


def UnTransform(ProcDat):
    
    NewDat = ProcDat.copy(deep=True)
    
    for i in np.arange(len(NewDat)):
        for j in np.arange(len(NewDat.iloc[i,:])):
            #print(ProcDat.iloc[i,j])
            
            if NewDat.iloc[i,j] > 0.0 :
                NewDat.iloc[i,j] = np.log(1.0+NewDat.iloc[i,j])
            elif NewDat.iloc[i,j] < 0.0:
                NewDat.iloc[i,j] = -1.0*np.log(1.0-NewDat.iloc[i,j])
    
    return(NewDat)
   
RawCleanUnReg = UnTransform(CleanUnReg)
RawCleanLogo = UnTransform(CleanLogo)


# In[144]:


fig3, ax3 = plt.subplots(2)
fig3.set_size_inches(7,2)

#TestSeq="VKV"+('-'*(len(CleanUnReg)-3))
TestSeq="L"*len(CleanUnReg)



#logomaker.Logo(CleanUnReg, ax =ax3[0], color_scheme=aa_colors)
#L1=logomaker.Logo(CleanUnReg, ax =ax3[0], color_scheme='dimgrey')

L1=logomaker.Logo(CleanUnReg, ax =ax3[0], color_scheme='lightgrey')
#L1.style_glyphs_in_sequence(sequence=ModSeq,color=[1,0.69,0.52]) 
#L1.style_glyphs_in_sequence(sequence=FreeSeq,color=[0.611,1,0.611])
L1.style_glyphs_in_sequence(sequence=ModSeq,color=[1,0.58,0.18]) 
L1.style_glyphs_in_sequence(sequence=FreeSeq,color=[0.18,0.827,0.424])





L2=logomaker.Logo(CleanLogo, ax =ax3[1], color_scheme=aa_colors)
#L2.style_glyphs_in_sequence(sequence=ModSeq,color=[1,0.69,0.52]) 
#L2.style_glyphs_in_sequence(sequence=FreeSeq,color=[0.611,1,0.611]) 
L2.style_glyphs_in_sequence(sequence=ModSeq,color=[1,0.58,0.18]) 
L2.style_glyphs_in_sequence(sequence=FreeSeq,color=[0.18,0.827,0.424])



ax3[0].set_ylim(-.99,0.99)
ax3[1].set_ylim(-.99,0.99)
ax3[0].set_xticks([])
ax3[1].set_xticks([])


fig3.tight_layout()
#Add X and Y labels


fig3.subplots_adjust(hspace=0.05)
fig3.savefig('Logo_UnRegvsReg_202204051426.png',dpi=600)


# In[121]:


BasicCleanLogo = CleanLogo.copy(deep=True)

for i in np.arange(len(BasicCleanLogo)):
    for j in np.arange(len(BasicCleanLogo.iloc[i,:])):
        
        if BasicCleanLogo.iloc[i,j] > 0.0:
            
            BasicCleanLogo.iloc[i,j] = BasicCleanLogo.iloc[i,j] + 1.0 
            


# In[122]:


CleanLogo


# In[123]:


fig3, ax3 = plt.subplots(3)
fig3.set_size_inches(9,3)

#TestSeq="VKV"+('-'*(len(CleanUnReg)-3))
TestSeq="L"*len(CleanUnReg)



#logomaker.Logo(CleanUnReg, ax =ax3[0], color_scheme=aa_colors)
#L1=logomaker.Logo(CleanUnReg, ax =ax3[0], color_scheme='dimgrey')

L1=logomaker.Logo(RawCleanUnReg, ax =ax3[0], color_scheme='lightgrey')
#L1.style_glyphs_in_sequence(sequence=ModSeq,color=[1,0.69,0.52]) 
#L1.style_glyphs_in_sequence(sequence=FreeSeq,color=[0.611,1,0.611])
L1.style_glyphs_in_sequence(sequence=ModSeq,color=[1,0.58,0.18]) 
L1.style_glyphs_in_sequence(sequence=FreeSeq,color=[0.18,0.827,0.424])





L2=logomaker.Logo(RawCleanLogo, ax =ax3[1], color_scheme=aa_colors)
#L2.style_glyphs_in_sequence(sequence=ModSeq,color=[1,0.69,0.52]) 
#L2.style_glyphs_in_sequence(sequence=FreeSeq,color=[0.611,1,0.611]) 
L2.style_glyphs_in_sequence(sequence=ModSeq,color=[1,0.58,0.18]) 
L2.style_glyphs_in_sequence(sequence=FreeSeq,color=[0.18,0.827,0.424])


L3=logomaker.Logo((BasicCleanLogo), ax =ax3[2], color_scheme=aa_colors) #CleanLogo+1 is easy when all positive Not True!!
#L2.style_glyphs_in_sequence(sequence=ModSeq,color=[1,0.69,0.52]) 
#L2.style_glyphs_in_sequence(sequence=FreeSeq,color=[0.611,1,0.611]) 
L3.style_glyphs_in_sequence(sequence=ModSeq,color=[1,0.58,0.18]) 
L3.style_glyphs_in_sequence(sequence=FreeSeq,color=[0.18,0.827,0.424])


a#x3[0].set_ylim(-0.8,1.0)
#ax3[1].set_ylim(-0.7,0.84)
ax3[0].set_xticks([])
ax3[1].set_xticks([])
fig3.tight_layout()
#Add X and Y labels


fig3.subplots_adjust(hspace=0.05)
#fig3.savefig('Logo_R1vR3_202112141409X.png',dpi=600)


# In[124]:


fig3, ax3 = plt.subplots(2)
fig3.set_size_inches(8,3)

#TestSeq="VKV"+('-'*(len(CleanUnReg)-3))
TestSeq="L"*len(CleanUnReg)



#logomaker.Logo(CleanUnReg, ax =ax3[0], color_scheme=aa_colors)
#L1=logomaker.Logo(CleanUnReg, ax =ax3[0], color_scheme='dimgrey')

L1=logomaker.Logo(RawCleanUnReg, ax =ax3[0], color_scheme='lightgrey')
#L1.style_glyphs_in_sequence(sequence=ModSeq,color=[1,0.69,0.52]) 
#L1.style_glyphs_in_sequence(sequence=FreeSeq,color=[0.611,1,0.611])
L1.style_glyphs_in_sequence(sequence=ModSeq,color=[1,0.58,0.18]) 
L1.style_glyphs_in_sequence(sequence=FreeSeq,color=[0.18,0.827,0.424])





L2=logomaker.Logo(RawCleanLogo, ax =ax3[1], color_scheme=aa_colors)
#L2.style_glyphs_in_sequence(sequence=ModSeq,color=[1,0.69,0.52]) 
#L2.style_glyphs_in_sequence(sequence=FreeSeq,color=[0.611,1,0.611]) 
L2.style_glyphs_in_sequence(sequence=ModSeq,color=[1,0.58,0.18]) 
L2.style_glyphs_in_sequence(sequence=FreeSeq,color=[0.18,0.827,0.424])


#ax3[0].set_ylim(-0.8,1.0)
#ax3[1].set_ylim(-0.7,0.84)
ax3[0].set_xticks([])
#ax3[1].set_xticks([])
fig3.tight_layout()
#Add X and Y labels


fig3.subplots_adjust(hspace=0.05)
ax3[1].set_xticks([])
#fig3.savefig('Logo_R1vR3_202201272253UnTran.png',dpi=600)


# In[125]:


print(max([max(RawCleanLogo.iloc[i,:]) for i in np.arange(len(RawCleanLogo))]))

print(max([max(RawCleanUnReg.iloc[i,:]) for i in np.arange(len(RawCleanUnReg))]))


# In[126]:



    
MaxReg = [max(CleanLogo.iloc[i,:]) for i in np.arange(len(CleanLogo))]

MaxUnreg = [max(CleanUnReg.iloc[i,:]) for i in np.arange(len(CleanUnReg))]


# In[127]:


print(max(MaxReg),max(MaxUnreg))


# In[128]:


print(min([min(CleanLogo.iloc[i,:]) for i in np.arange(len(CleanLogo))]),min([min(CleanUnReg.iloc[i,:]) for i in np.arange(len(CleanUnReg))]))


# In[129]:



pos=0

PosMap_Mod={}

for i in np.arange(len(RatZP3_Alig)):
    if RatZP3_Alig[i]=='-':
        continue
    else:
        #print(i,pos)
        PosMap_Mod[i+1]=(pos+1)
        pos+=1
        


# In[130]:


PosMap_Mod


# In[131]:


pos=0

PosMap_Free={}

for i in np.arange(len(RatZP2ZPN_Alig)):
    if RatZP2ZPN_Alig[i]=='-':
        continue
    else:
        #print(i,pos)
        PosMap_Free[i+1]=(pos+1)#Add 46 to match pdb 
        pos+=1
#PosMap_Free


# ### PosMap_Free

# In[132]:


for i in np.arange(len(BetaY)):
    
    #print(PosX[i])
    
    if BetaY[i] > 0:
        print('Mod',PosMap_Mod[PosX[i]])
    else:
        print('Free',PosMap_Free[PosX[i]])
    


# In[133]:


BetaY


# In[134]:


PosX


# In[135]:


print(MouseZP3_UnAlig,'\n',RatZP3_Alig)


# In[136]:


CleanLogo


# In[137]:


ww_df.columns


# In[138]:


residues


# In[139]:


len(LR3)


# In[ ]:





# In[ ]:





# In[ ]:




