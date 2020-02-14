# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
fileName  = 'h2o.dat'

# functions

def splitArrays2(channelNumber,array):
    allArrays = []
    allArrays.append([channelNumber[200:500],array[200:500]])
    allArrays.append([channelNumber[500:650],array[500:650]])
    allArrays.append([channelNumber[650:750],array[650:750]])
    allArrays.append([channelNumber[720:1000],array[720:1000]])
    for entry in allArrays:
        plt.plot(entry[0],entry[1])

    # colour codes the plots
    return allArrays

def findStatPoints(x,y):
    peaks = np.where((y[1:-1] > y[0:-2]) * (y[1:-1] > y[2:]))[0] + 1
    dips = np.where((y[1:-1] < y[0:-2]) * (y[1:-1] < y[2:]))[0] + 1

    
    return [[x[peaks], y[peaks]],[x[dips], y[dips]]]

def convertToNormalArray(numpyArray):
    temp = []
    for entry in numpyArray:
        temp.append(entry)
    return temp

def findZeros(channelNumber,arrays):
    zeroes = []
    for array in arrays:
        # work with one array at a time
        #print("separated arrays")
        #print(array)
        points = findStatPoints(array[0],array[1])
        line = [ channelNumber[points[1][0][0]:points[0][0][0]] ,arrays[points[1][0][0]:points[0][0][0]] ]
        # line[0] are the x vales while
        # line[1] are the y values
        plt.plot(line[0],line[1])
        # get a line and interpolate its zero point. append to the zeroes array
        
    return zeroes
    # return
        


def findZero(array):
    maxima, minima = findStatPoints(array[0],array[1])
        

h20data = np.genfromtxt(fileName, delimiter = ',')
channelNumber  = np.arange(1,len(h20data)+1)
waveNumber = [12683.782,12685.540,12685.769,12687.066]     # water
plt.plot(x,h20data)

peaksArrays = splitArrays(channelNumber, h20data) # data arrays for the peaks
zeroes = findZeros(channelNumber,peaksArrays)
###################################



'''

h20dataPeaks = []




fit_parameters = np.polyfit(channelNumber,h20dataPeaks,1)
fit_m = fit_parameters[0]
fit_c = fit_parameters[1]
print('Linear np.polyfit of y = m*x + c')
print('Gradient  m = {:04.2f}'.format(fit_m)),
print('Intercept c = {:04.2f}'.format(fit_c))
'''


plt.show()
