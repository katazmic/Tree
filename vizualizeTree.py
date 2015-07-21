# -*- coding: utf-8 -*-

from __future__ import division                                                                                                                                                    

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

import math
import scipy as sp
import numpy as np


data = np.loadtxt('rsltscmp.txt');
#data = np.loadtxt('rsltscmpEx.txt');
input = np.loadtxt('inputsPlot.txt');


N = (int) (input[0])

k = []
for i in data[1]:
    k.append(i);

i=1
while i<len(k):
    if k[i]==k[i-1]:
        k.pop(i);
        i=i-1
    i=i+1
        



max_range = 0
i= 0
while data[1][-1]==data[1][-1-i]:
    max_range = max_range+1;
    i = i+1    



X = range(max_range)
X, kn = np.meshgrid(X, k)


av = data[-300]
for i in data[-300:-1]:
    av =av + i

av= data[-1]
Mat = np.zeros((len(k),max_range))


l= 0
for i in range(10):#len(k)):
    h=0;
    j=0;
    while j <max_range:
        h=h+1
        while(j<h*(max_range/2**i) and j<max_range ):
            Mat[i][j] = av[l]
            j=j+1    
        l=l+1

i=i+1
while i < len(k):
    j=0;
    while j <max_range:
        Mat[i][j] = av[l]
        l=l+1
        j=j+1
    i=i+1

plt.plot(Mat[10],'.r')
plt.plot(Mat[11],'.b')
plt.plot(Mat[12],'.g')
plt.show()


########## plotting loglog 2D surfaces #########



X = range(max_range)
X, lkn = np.meshgrid(X, np.log(k))

lMat = np.zeros((len(k),max_range))

for i in range(len(k)):
    for j in range(max_range):
        lMat[i][j] = np.log(Mat[i][j])




#fig = plt.figure()
#ax = fig.gca(projection='3d')


#surf = ax.plot_surface(X, lkn, lMat,cmap=cm.coolwarm,rstride=1, cstride=1,  antialiased=False)

#ax.set_xlabel('X') 
#ax.set_ylabel('log(k)') 
#ax.set_zlabel('$\phi_{k,n}^2$')

#plt.show()

################


X = range(max_range)
lkn, X = np.meshgrid(np.log(k), X)

lMat = np.zeros((max_range,len(k)))
avSp = np.zeros(len(k))
for i in range(len(k)):
    for j in range(max_range):
        lMat[j][i] = np.log(Mat[i][j])
        avSp[i] = avSp[i] + Mat[i][j]

plt.loglog(k,k*avSp/len(k))
plt.xlabel('k')
plt.ylabel('$E(k) = k\phi^2$')
plt.show()




fig = plt.figure()
ax = fig.gca(projection='3d')


surf = ax.plot_surface(lkn, X, lMat,cmap=cm.coolwarm,rstride=1, cstride=1,  antialiased=False)

ax.set_ylabel('X') 
ax.set_xlabel('log(k)') 
ax.set_zlabel('$\phi_{k,n}^2$')

plt.show()



