'''
Description:
Author: Tomasz Neska
Date: 12.02.2020
'''
# initialization
import string
from math import *
import numpy as np
import matplotlib.pyplot as plt
import random

plt.rcParams.update({'font.size': 14})
plt.style.use('default')

# constants
c = 299792458  # speed of light in ms^-1
h = 6.626E-34  # Planck's constant in Js
water_file = 'h2o.dat'
experiment_file = 'c2h2.dat'


# functions
def getData(filename):
    data = np.genfromtxt(filename, delimiter=',')
    temp = []
    for counter in range(0, len(data), 1):
        if data[counter] == np.NaN:
            temp.append(counter)  # appends the index values that contain invalid entries

    for entry in temp:  # deletes the data that is not valid
        data = np.delete(data, entry)

    return data


def plotDataSimple(title, xAxisTitle, yAxisTitle, x, y, error_y):
    print("printing data")
    plt.rcParams["figure.figsize"] = (6, 3)
    plt.figure()
    plt.plot(x, y, "ro")
    plt.errorbar(x, y, error_y, )
    plt.xlabel(xAxisTitle)
    plt.ylabel(yAxisTitle)
    plt.title(title)
    plt.show()

############################ uselesss function ####################################################################
def findStatPoints(x, y):
    # this function finds the maxima and minima points by comparing the neightbours
    peaks = np.where((y[1:-1] > y[0:-2]) * (y[1:-1] > y[2:]))[0] + 1
    dips = np.where((y[1:-1] < y[0:-2]) * (y[1:-1] < y[2:]))[0] + 1

    return [[x[peaks], y[peaks]], [x[dips], y[dips]]]
####################################################################################################################

def splitArrays(channelNumber, array):
    # splits the array into the peaks for the h2o data set
    allArrays = []
    allArrays.append([channelNumber[200:500], array[200:500]])
    allArrays.append([channelNumber[500:650], array[500:650]])
    allArrays.append([channelNumber[650:750], array[650:750]])
    allArrays.append([channelNumber[720:1000], array[720:1000]])
    for entry in allArrays:
        plt.plot(entry[0], entry[1])

    # "colour codes" the plots
    return allArrays


def search(arr, x):
    # linear search function
    for i in range(len(arr)):
        if arr[i] == x:
            return i
    return -1

def findZeroes(array):
    x = array[0]
    array = array[1] # splits the variable into two arrays
    maximum = [x[search(array, np.max(array))], np.max(array)] # finds the maximum point
    minimum = [x[search(array, np.min(array))], np.min(array)] # finds the minimum point

def main():
    h2odata = getData(water_file)
    channel_number = np.arange(1, len(h2odata) + 1)  # get the x-coordiantes
    wave_number_water = [12683.782, 12685.540, 12685.769, 12687.066]
    peaks_arrays = splitArrays(channel_number, h2odata)  # data arrays for the peaks
    for peak in peaks_arrays: # peak has both x and y coordinates
        print(peak)
        #zeroes = findZeroes(peak):

main()
