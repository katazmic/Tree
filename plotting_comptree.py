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


dataT5 = np.loadtxt('rslts_tree52D.txt');


k = dataT5[n1*N:(n1+1)*N,0];
phT5 = dataT5[n1*N:(n1+1)*N,1];

for nt in range(Nav):
    phT5 = phT5 + dataT5[(n1+1+nt)*N:(n1+2+nt)*N,1];

plt.loglog(k,phT5,'+b')

plt.loglog(k,ph,'o',k,1.e5*k**(-8./3.),k,1.e7*k**(-4.),'-r')


plt.show()
