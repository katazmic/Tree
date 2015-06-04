# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import math
import scipy as sp
import numpy as np



data = np.loadtxt('rslts.txt');

input = np.loadtxt('inputsPlot.txt');

N = (int) (input[0]);
n1 = (int) (input[1]);
Nav = (int) (input[2]);
k = data[n1*N:(n1+1)*N,0];
ph = data[n1*N:(n1+1)*N,1]; 


for n in range(Nav-1):
    ph = ph + data[(n1+n)*N:(n1+1+n)*N,1];


plt.loglog(k,ph,'o',k,2e5*k**(-8./3.),k,5e6*k**(-4.),'-r')
plt.xlabel('$k_n$')
plt.ylabel('$\Phi_n^2$')

plt.show()
