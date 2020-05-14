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

a = np.array(bins)
a = a.tolist()
print(a)

#################################################
figure2 = plt.figure()
ax11 = figure2.add_subplot(111)

ax11.scatter(array2[2][0], array2[2][1])
ax11.plot(array2[2][0] * self.gradients[2][0][0])
ax11.errorbar(array2[2][0], array2[0][1], yerr=array2[0][2])