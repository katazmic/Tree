# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import math
import scipy as sp
import numpy as np



data = np.loadtxt('rslts_tree12D.txt');


N = 30;
n1 = 15000;
Nav = 400;
k = data[n1*N:(n1+1)*N,0];
ph = data[n1*N:(n1+1)*N,1]; 

for n in range(Nav):
    ph = ph + data[(n1+1+n)*N:(n1+2+n)*N,1];

dataT2 = np.loadtxt('rslts_tree22D.txt');


k = dataT2[n1*N:(n1+1)*N,0];
phT2 = dataT2[n1*N:(n1+1)*N,1];

for nt in range(Nav):
    phT2 = phT2 + dataT2[(n1+1+nt)*N:(n1+2+nt)*N,1];

dataT3 = np.loadtxt('rslts_tree32D.txt');


k = dataT3[n1*N:(n1+1)*N,0];
phT3 = dataT3[n1*N:(n1+1)*N,1];

for nt in range(Nav):
    phT3 = phT3 + dataT3[(n1+1+nt)*N:(n1+2+nt)*N,1];

dataT4 = np.loadtxt('rslts_tree42D.txt');


k = dataT4[n1*N:(n1+1)*N,0];
phT4 = dataT4[n1*N:(n1+1)*N,1];

for nt in range(Nav):
    phT4 = phT4 + dataT4[(n1+1+nt)*N:(n1+2+nt)*N,1];


dataT5 = np.loadtxt('rslts_tree52D.txt');


k = dataT5[n1*N:(n1+1)*N,0];
phT5 = dataT5[n1*N:(n1+1)*N,1];

for nt in range(Nav):
    phT5 = phT5 + dataT5[(n1+1+nt)*N:(n1+2+nt)*N,1];

dataT6 = np.loadtxt('rslts_tree62D.txt');


k = dataT6[n1*N:(n1+1)*N,0];
phT6 = dataT6[n1*N:(n1+1)*N,1];

for nt in range(Nav):
    phT6 = phT6 + dataT6[(n1+1+nt)*N:(n1+2+nt)*N,1];

dataT7 = np.loadtxt('rslts_tree72D.txt');


k = dataT7[n1*N:(n1+1)*N,0];
phT7 = dataT7[n1*N:(n1+1)*N,1];

for nt in range(Nav):
    phT7 = phT7 + dataT7[(n1+1+nt)*N:(n1+2+nt)*N,1];


dataT8 = np.loadtxt('rslts_tree82D.txt');


k = dataT8[n1*N:(n1+1)*N,0];
phT8 = dataT8[n1*N:(n1+1)*N,1];

for nt in range(Nav):
    phT8 = phT8 + dataT8[(n1+1+nt)*N:(n1+2+nt)*N,1];

dataT9 = np.loadtxt('rslts_tree92D.txt');


k = dataT9[n1*N:(n1+1)*N,0];
phT9 = dataT9[n1*N:(n1+1)*N,1];

for nt in range(Nav):
    phT9 = phT9 + dataT9[(n1+1+nt)*N:(n1+2+nt)*N,1];




dataT10 = np.loadtxt('rslts_tree102D.txt');

k = dataT10[n1*N:(n1+1)*N,0];
phT10 = dataT10[n1*N:(n1+1)*N,1];

for nt in range(Nav):
    phT10 = phT10 + dataT10[(n1+1+nt)*N:(n1+2+nt)*N,1];




dataT11 = np.loadtxt('rslts_tree112D.txt');


k = dataT11[n1*N:(n1+1)*N,0];
phT11 = dataT11[n1*N:(n1+1)*N,1];

for nt in range(Nav):
    phT11 = phT11 + dataT11[(n1+1+nt)*N:(n1+2+nt)*N,1];







plt.loglog(k,phT2/ph,'+-b')
plt.loglog(k,phT3/ph,'+-g')
plt.loglog(k,phT4/ph,'+-y')
plt.loglog(k,phT5/ph,'+-r')
plt.loglog(k,phT6/ph,'-m')
plt.loglog(k,phT7/ph,'-c')
plt.loglog(k,phT8/ph,'-r')
plt.loglog(k,phT9/ph,'-b')
plt.loglog(k,phT10/ph,'-y')
plt.loglog(k,phT11/ph,'-g')


plt.loglog(k,ph/ph,'o')#,k,1.e5*k**(-8./3.),k,1.e7*k**(-4.),'-r')


plt.show()

plt.loglog(k,ph,'r',k,1.e5*k**(-8./3.),k,1.e7*k**(-4.))
plt.loglog(k,phT11,'-b')
plt.show()
