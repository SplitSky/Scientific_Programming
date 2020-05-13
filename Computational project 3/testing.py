import string
from math import *
import numpy as np
import matplotlib.pyplot as plt
import random
import time
import cmath
import json
from scipy import optimize
from mpl_toolkits.mplot3d import Axes3D


def generateArray(first, last, N):
    array = []
    for counter in range(0, N, 1):
        temp = random.randint(first, last)
        array.append(temp)

    return array


def testRandomGenerator_withCounting(array, first, last, N):
    counters = []
    for i in range(0, last - first):
        counters.append([i, 0])

    counters = np.array(counters)
    counters += first
    # add first to all entries in array

    for i in range(0, len(counters)):
        for number in array:
            if number == counters[i][0]:
                counters[i][1] += 1

    temp = []
    x = []
    for entry in counters:
        x.append(entry[0])
        temp.append(entry[1])

    return x, temp


def testRandomGenerator(first, last, N):
    x = np.arange(first, last+2) # too get all the "bins" for the histogram
    print(x)
    array = generateArray(first, last, N)
    print(array)
    plt.hist(array, x)
    mean = np.mean(array)
    std = np.std(array)
    print("The mean of the generated array is: {0:9f}".format(mean))
    print("The standard deviation is: {0:9f}".format(std))
    print("The square root of the mean is: {0:9f}".format(np.sqrt(mean)))



first = 0
last = 100
N = 10000
a = generateArray(first, last, N)
testRandomGenerator(first, last, N)

plt.show()
