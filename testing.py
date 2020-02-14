# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 12:35:35 2020

@author: d06939tn
"""
import numpy as np
import matplotlib.pyplot as plt
import math as m

def search(arr, x):
    for i in range(len(arr)): 
        if arr[i] == x: 
            return i 
    return -1

def extractPeaks(x,array,n):
    # n is the expected number of peaks
    peaks = []
    dips = []
    for counter in range(0,n,1):
        maximumIndex = search(array,np.max(array))
        peaks.append([x[maximumIndex],np.max(array)])
        minimumIndex = search(array,np.min(array))
        dips.append([x[minimumIndex],np.min(array)])
        array = np.delete(array,minimumIndex)
        array = np.delete(array,maximumIndex)
        
    return [peaks,dips,array]
        
array = np.array([0,1,3,2,2,1,1,2,200,2,-200,1,1,2,1,100,1,-100,2,2,1,1,1,1])
x = np.arange(1,len(array)+1)
plt.plot(x,array)

array = extractPeaks(x,array,2)[2]
x = np.arange(1,len(array)+1)

peaks = extractPeaks(x,array,2)[0]
dips  = extractPeaks(x,array,2)[1]

print(peaks)
print(dips)

plt.plot(x,array)

plt.show()