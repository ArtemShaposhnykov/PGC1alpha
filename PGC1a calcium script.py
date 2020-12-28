# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 23:46:15 2020

@author: Arty
"""

import numpy as np
import matplotlib.pyplot as plt

F = np.genfromtxt('C:\\Users\\Arty\\tricalcalcium\\2TT.csv', delimiter=',')

a= len(F)    #counter of the cells traces
b= len (F[1]) # counter of the time point in the single trace

pik = [[0]*b for i in range(a)] 
pik = np.array(pik, dtype=float)
pik1 = [[0]*b for i in range(a)]
pik1 = np.array(pik, dtype=float)
bl = [[0]*b for i in range(a)]
bl = np.array(pik, dtype=float)
tau = [[0]*6 for i in range(a)]
  #tau == [[0]all spike;[1] spike of alive cells;
  # [2] frequence [3]real tau]
tau = np.array(tau, dtype = float)

u = [1]*a
z = [1]*a
uR = [0]*a

#find maximum and position of the maximum
for x in range (0, a):
    #print(x)
    z[x] = np.amax(F[x])
    
 #normalization to Fmax   
for x in range (0, a):
    for j in range (0, b):
        F[x,j] = F[x,j]/z[x]

#filling of the pik array with avarage value of the neibor +-5wells
for x in range (0, a):
    for j in range (7, b-7):
        pik[x,j]= ((F[x,j-7]+F[x,j-6]+F[x,j-5]+ F[x,j-4]+F[x,j-3]+
         F[x,j-2]+F[x,j-1]+F[x,j]+F[x,j+1]+F[x,j+2]+F[x,j+3]+F[x,j+4]
         +F[x,j+5]+F[x,j+6]+F[x,j+7])/15)

for x in range (0, a):
    for j in range (7, b-7):
        pik1[x,j]=((F[x,j-7]+F[x,j-6]+F[x,j-5]+ F[x,j-4]+F[x,j-3]+
         F[x,j-2]+F[x,j-1]+F[x,j]+F[x,j+1]+F[x,j+2]+F[x,j+3]+F[x,j+4]
         +F[x,j+5]+F[x,j+6]+F[x,j+7])/15)
        bl[x,j]=((F[x,j-7]+F[x,j-6]+F[x,j-5]+ F[x,j-4]+F[x,j-3]+
         F[x,j-2]+F[x,j-1]+F[x,j]+F[x,j+1]+F[x,j+2]+F[x,j+3]+F[x,j+4]
         +F[x,j+5]+F[x,j+6]+F[x,j+7])/15)
#pick detection
for x in range (0, a):
    for j in range (7, b-7):
        if (F[x,j]>F[x,j-1] and F[x,j]>F[x,j+1] 
        and F[x,j]>F[x,j-2] and F[x,j]>F[x,j+2]
        and F[x,j]>F[x,j-3] and F[x,j]>F[x,j+3] 
        and F[x,j]>F[x,j-4] and F[x,j]>F[x,j+4] 
        and F[x,j]>F[x,j-5] and F[x,j]>F[x,j+5]): 
       # and (((F[x,j]-pik[x,j])/pik[x,j])>0.1)):
            pik[x,j] = F[x,j]
           
          # tau[x,0]=tau[x,0]+1
        else:
            pik[x,j] = 0

for x in range (0, a):
    for j in range (7, b-7):
        if (((F[x,j]-pik1[x,j])/pik1[x,j] )>0.25   or (F[x,j]>0.6)): 
            pik1[x,j] = F[x,j]
        else:
            pik1[x,j] = 0
       
for x in range (0, a):
    for j in range (7, b-7):
        if (pik1[x,j] == pik[x,j]): 
            pik[x,j] = pik[x,j]
        else:
            pik[x,j] = 0  
            
n=1
z=1
pik1 = [[0]*b for i in range(a)]
pik1 = np.array(pik, dtype=float)
for x in range (0,a):
    for j in range (7, b-7):
        if pik[x,j] !=0:
            for g in range (j+1, b-3):
                if ((0.1*(F[x,j]-bl[x,j])+bl[x,j])>0)and(F[x,g] > (0.1*(F[x,j]-bl[x,j])+bl[x,j]) and pik[x,g] == 0) :
                    z=z+1
                else:
                    pik1[x,j]=z
                    z=1
                    n=n+1
                    break
# z, n - parameters to define of the part of the ca trace for the analysis
z=800
n=1800
for x in range (0,a):
    for j in range (z, n):
        if pik1[x,j] != 0:
            tau[x,0]=tau[x,0]+1

for x in range (0,a):
    for j in range (z, n):
        tau[x,5]=tau[x,5]+F[x,j]
#summ of all real tau            
for x in range (0,a):
    for j in range (z, n):
        if pik1[x,j] != 0:
            tau[x,1]=tau[x,1]+pik1[x,j] 
#real tau
for x in range (0,a):
    if tau[x,0] != 0:
        tau[x,2]=tau[x,1]/tau[x,0]
        
for x in range (0,a):
    for j in range (z, n):
        tau[x,3]=tau[x,3]+pik[x,j] 

for x in range (0,a):
    for j in range (z, n):
        if pik1[x,j] != 0:
           tau[x,4]=tau[x,4]+bl[x,j]

# visual control, parameters need to be adjusted
k=0
plt.xlim(1000,3000)
plt.plot(pik[1], 'bo')
plt.plot(F[1])
#plt.plot(F[19])
#plt.plot(F[10])
plt.show()
for g in range (4,4):
    plt.plot(F[g])
plt.show()        
        
        
