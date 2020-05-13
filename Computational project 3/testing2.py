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
from numba import jit
jit(nopython=True)

a = 33*0.6
bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
random_set = np.random.random(100000000)
random_set *= 100
random_set = -1 * 1/a * np.log(random_set)
bin_frequencies, bin_locations = np.histogram(random_set, bins=bins)
plt.hist(random_set, bins=bins)
std = np.std(bin_frequencies)
mean = np.mean(bin_frequencies)
plt.show()
print("The standard deviation is: {0:9f}".format(std))
print("The mean is: {0:9f}".format(mean))
