# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 16:49:47 2025

@author: Usuario
"""

import numpy as np
import matplotlib.pyplot as plt



h = 1
eigenvalue = 30
eta_beta = 0.0
beta_mu = 0.034
mu = 0.5

for beta_mu in np.linspace(0.034,0.06,50):
    
    beta = beta_mu*mu
    eta = eta_beta*beta
    
    p = (beta * eigenvalue + eta - mu)/((1 - mu )*(beta * eigenvalue + eta ))
    q = 1 - beta * eigenvalue * p - eta * p
    
    W = ((1-q)+mu*q*(q*(1-mu))**h)/(1-q*(1-mu))
    
    r = (W*p-p**2)/(p-p**2)
    plt.plot(beta_mu,r,marker='o')
plt.show()
    