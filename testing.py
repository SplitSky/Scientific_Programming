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


def extractPeaks(x, array, n):
    # n is the expected number of peaks
    peaks = []
    dips = []
    for counter in range(0, n, 1):
        maximumIndex = search(array, np.max(array))
        peaks.append([x[maximumIndex], np.max(array)])
        minimumIndex = search(array, np.min(array))
        dips.append([x[minimumIndex], np.min(array)])
        array = np.delete(array, minimumIndex)
        array = np.delete(array, maximumIndex)

    return [peaks, dips, array]


x = [651, 652, 653, 654, 655, 656, 657, 658, 659, 660, 661, 662, 663,
     664, 665, 666, 667, 668, 669, 670, 671, 672, 673, 674, 675, 676,
     677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689,
     690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702,
     703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715,
     716, 717, 718, 719, 720, 721, 722, 723, 724, 725, 726, 727, 728,
     729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741,
     742, 743, 744, 745, 746, 747, 748, 749, 750]

array = np.array([5.1003, 5.0788, 5.0669, 5.0645, 5.0693, 5.0287, 4.9904, 4.9713,
                  4.9403, 4.8495, 4.7038, 4.5127, 4.2929, 4.1066, 3.9824, 4.0301,
                  4.3168, 4.7946, 5.3655, 5.8289, 6.1156, 6.1801, 6.0439, 5.8385,
                  5.6522, 5.5064, 5.3679, 5.2819, 5.2293, 5.1648, 5.1553, 5.1481,
                  5.1409, 5.1242, 5.1194, 5.1171, 5.0764, 5.0884, 5.0956, 5.0836,
                  5.0645, 5.0693, 5.0812, 5.0621, 5.0597, 5.0478, 5.0597, 5.0597,
                  5.0573, 5.0526, 5.0454, 5.0717, 5.0717, 5.0645, 5.0621, 5.0884,
                  5.0741, 5.0836, 5.0812, 5.0741, 5.0788, 5.0741, 5.0764, 5.0836,
                  5.0645, 5.0812, 5.0764, 5.0645, 5.0597, 5.0549, 5.0597, 5.0526,
                  5.0573, 5.0526, 5.0478, 5.0502, 5.0597, 5.0621, 5.0621, 5.0597,
                  5.0597, 5.0597, 5.0621, 5.0573, 5.043, 5.0549, 5.0573, 5.0597,
                  5.0526, 5.0358, 5.0382, 5.043, 5.0502, 5.0382, 5.043, 5.0478,
                  5.0526, 5.0502, 5.0478, 5.043])
x = np.array(x)
array = array - 5
plt.plot(x, array)
plt.plot(x, array)
#

temp = np.min(array)
index2 = search(array, temp)
minimum = [x[index2], temp]

temp = np.max(array)
index1 = search(array, temp)
maximum = [x[index1], temp]

cut_array = array[index2:index1]
x_temp = x[index2:index1]

plt.plot(x_temp, cut_array)

temp = np.linspace(x_temp[0], x_temp[len(x_temp) - 1], 100)
interpolated_y = np.interp(temp, x_temp, cut_array)

fitted_parameters = np.polyfit(temp, interpolated_y, 1)
zero = -1 * fitted_parameters[1] / fitted_parameters[0]
plt.plot(temp, fitted_parameters[0] * temp + fitted_parameters[1], "g")
print(zero)

zero = np.rint(zero)
zero = int(zero)

plt.plot(zero, 0, "or")
print(zero)
plt.plot(x, x * 0)

#
plt.show()
