# -*- coding: utf-8 -*-
"""
Author: Tomasz Neska
Date: 06/03/2020
"""
# initialisation
import string
from math import *
import numpy as np
import matplotlib.pyplot as plt
import random
import time
import cmath

plt.rcParams.update({'font.size': 14})
plt.style.use('default')
figure = plt.figure()
plt.rcParams.update({'errorbar.capsize': 2})

m = 1.2
k = 0.8
b = 0.08

T = 221
h = 0.01
nsteps = int(np.rint(T/h))

def solver(x_init,v_init):
    A = 0
    B = 0
    temp = (b**2/(4*m**2))
    if (k/m) > temp: # imaginary
        marker = 3
        K = complex(-b/2*m, np.sqrt((k/m) - b**2 / (4*m**2)))
        a = K.real()
        b = np.abs(K.imag())
        #initial conditions
        A = (v_init-x_init*b)/(a-b)
        B = (v_init-x_init*a)/(b-a)
        return a,b,marker,A,B
    elif (k/m) = temp:
        marker = 2 # repeated real solutions
        K = -b / 2m
        # initial conditions
        A = x_init
        return K,0,marker,A,B
    elif (k/m) < temp:
        marker = 1
        a = -b/2*m + np.sqrt((k/m) - b**2 / (4*m**2))
        b = -b/2*m - np.sqrt((k/m) - b**2 / (4*m**2))
        # initial conditions
        A = v_init / b
        B = x_init
        return a,b,marker,A,B
    
def Position_t(t):
    k_1, k_2, marker,A,B = solver()
    if marker == 1:
        return A * np.exp(k_1*t) + B* np.exp(k_2*t)
    elif marker == 2:
        return A * np.exp(k_1)
    elif marker == 3:
        return np.exp(k_1) * (A*np.sin(k_2*t) + B*np.cos(k_2*t))

def initialCond(x_init, v_init, marker):
    

def Euler(steps,x_init, v_init):
    x = np.zeros(steps)
    v = np.zeros(steps)
    
    x[0] = x_init
    v[0] = v_init
    for i in range(steps-1):
        a = -(k/m)*x[i] - (b/m)*v[i]
        x[i+1] = x[i] + v[i]*h
        v[i+1] = v[i] + a*h
    
    plt.plot(x,'b')
    plt.plot(v,'r')
    return x,v
    
def Energy_t(x,v):
    return 1/2 * k * x**2 + 1/2 * m * v**2 

def plot(x,y):
    figure = plt.figure()
    axes_1 = figure.add_subplot(111)
    axes_1.plot(x,y)
    plt.xlabel("x")
    plt.ylabel("v")

def main():
    x,v = Euler(nsteps,2.3,0.3)
    E = Energy_t(x,v)
    plot(x,E)
    plot(v,E)
    
main()
plt.show()