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
water_file = 'h20.dat'
experiment_file = 'c2h2.dat'


# functions
def getData(filename):
    data = np.genfromtxt(filename, delimiter=',')
    temp = []
    for counter in range(0, len(data), 1):
        if data[counter] == np.NaN:
            temp.append(counter)

    for entry in temp:  # deletes the data that is invalid
        data = np.delete(data, entry)

    # add any other validation if necessary
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

def plotDataAdvanced():
    # this function

def main():
    data = getData(water_file)


# def main

main()
