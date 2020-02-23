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


def splitArrays(channelNumber, array, marker):
    # splits the array into the peaks for the h2o data set and the c2h2 data set
    allArrays = []
    if marker:
        allArrays.append([channelNumber[200:500], array[200:500]])
        allArrays.append([channelNumber[500:650], array[500:650]])
        allArrays.append([channelNumber[650:750], array[650:750]])
        allArrays.append([channelNumber[720:1000], array[720:1000]])
    else:
        allArrays.append([channelNumber[:150], array[:150]])
        allArrays.append([channelNumber[350:550], array[350:550]])
        allArrays.append([channelNumber[800:900], array[800:900]])
        allArrays.append([channelNumber[1200:1300], array[1200:1300]])
        allArrays.append([channelNumber[1580:1650], array[1580:1650]])

    # "colour codes" the plots
    return allArrays


def search(arr, x):
    # linear search function
    for i in range(len(arr)):
        if arr[i] == x:
            return i
    return -1


def findZero(x, array):
    temp2 = array - 5
    # finds and returns the channel number at which the section of the array passes a zero

    index2 = search(temp2, np.min(temp2))
    index1 = search(temp2, np.max(temp2))

    cut_array = temp2[index2:index1]
    x_temp = x[index2:index1]

    # plt.plot(x_temp, cut_array+5)

    temp = np.linspace(x_temp[0], x_temp[len(x_temp) - 1], 100)
    interpolated_y = np.interp(temp, x_temp, cut_array)

    fitted_parameters = np.polyfit(temp, interpolated_y, 1)
    zero = -1 * fitted_parameters[1] / fitted_parameters[0]

    return int(np.rint(zero))  # returns channel number

def linearFit(x,y):
    # Perform a linear fit and get the errors
    fit_parameters, fit_errors = np.polyfit(x, y, 1, cov=True)
    fit_m = fit_parameters[0]
    fit_c = fit_parameters[1]
    # Here, we (rather laboriously) explicitly define some variables so you can see
    # exactly which matrix element is which
    variance_m = fit_errors[0][0]
    variance_c = fit_errors[1][1]
    sigma_m = np.sqrt(variance_m)
    sigma_c = np.sqrt(variance_c)
    print('Linear np.polyfit of y = m*x + c')
    print('Gradient  m = {:04.10f} +/- {:04.10f}'.format(fit_m, sigma_m))
    print('Intercept c = {:04.10f} +/- {:04.10f}'.format(fit_c, sigma_c))
    return [[fit_m,sigma_m], [fit_c, sigma_c]]

def quadraticFit(x,y):
    # Perform a quadratic fit and get the errors
    fit_parameters, fit_errors = np.polyfit(x, y, 2, cov=True)
    fit_a = fit_parameters[0]
    fit_b = fit_parameters[1]
    fit_c = fit_parameters[2]
    # Here, we (rather laboriously) explicitly define some variables so you can see
    # exactly which matrix element is which
    variance_a = fit_errors[0][0]
    variance_b = fit_errors[1][1]
    variance_c = fit_errors[2][2]
    sigma_a = np.sqrt(variance_a)
    sigma_b = np.sqrt(variance_a)
    sigma_c = np.sqrt(variance_c)
    print('Quadratic fit of y = a*x**2 + b*x + c')
    print('Quadratic term  a = {:04.10f} +/- {:04.10f}'.format(fit_a, sigma_a))
    print('Linear term     b = {:04.10f} +/- {:04.10f}'.format(fit_b, sigma_b))
    print('Intercept       c = {:04.10f} +/- {:04.10f}'.format(fit_c, sigma_c))
    return [[fit_a,sigma_a],[fit_b,sigma_b],[fit_c, sigma_c]]

def main():
    h2o_data = getData(water_file)

    channel_number = np.arange(1, len(h2o_data) + 1)  # get the x-coordinates

    # plt.plot(channel_number, h2o_data)

    wave_number_water = [12683.782, 12685.540, 12685.769, 12687.066]
    peaks_arrays = splitArrays(channel_number, h2o_data, True)  # data arrays for the peaks
    zeroes = []  # the points where the graph crosses the 5V point

    for peak in peaks_arrays:
        temp = findZero(peak[0], peak[1])
        zeroes.append([temp, h2o_data[temp - 1]])

    temp = []
    for counter in range(0, len(zeroes)):
        temp.append(zeroes[counter][0])

    temp = np.array(temp)
    wave_number_water = np.array(wave_number_water)

    fitting_parameter_calibration = np.polyfit(temp, wave_number_water, 1)  # calibration fitting
    '''
    enter the channel number and get the wave number
    '''
    c2h2_data = getData(experiment_file)

    channel_number = np.arange(1, len(c2h2_data) + 1)  # get the x-coordinates

    peaks_arrays = splitArrays(channel_number, c2h2_data, False)
    zeroes = []

    for peak in peaks_arrays:
        temp = findZero(peak[0], peak[1])
        zeroes.append([temp, c2h2_data[temp - 1]])
    temp = []
    for counter in range(0, len(zeroes)):
        temp.append(zeroes[counter][0])

    temp = np.array(temp)
    wave_numbers_c2h2 = fitting_parameter_calibration[0] * temp + fitting_parameter_calibration[1]
    # values of the wave number from calibrated spectrum
    energy_c2h2 = c * h * wave_numbers_c2h2

    x = np.array([1., 2, 3, 4.1, 5, 6, 6.5, 8, 9.3, 10.7])
    y = 2*x**2 + 5*x + 7
    results = quadraticFit(x,y)
    #####################################################################################################
    '''
    m = np.array([3, 4, 5, 6, 7])
    fit_parameters, fit_errors = np.polyfit(m, energy_c2h2, 1, cov=True)
    fit_m = fit_parameters[0]
    fit_c = fit_parameters[1]

    variance_m = fit_errors[0][0]
    variance_c = fit_errors[1][1]
    sigma_m = np.sqrt(variance_m)
    sigma_c = np.sqrt(variance_c)
    print("Linear fit of the quantum number against the energy values of the spectrum for h2c2 is:")
    print("Gradient: m = {:04.10f} +/- {:04.10f}".format(fit_m, sigma_m))
    print('Intercept c = {:04.10f} +/- {:04.10f}'.format(fit_c, sigma_c))

    # linear fitting of the c2h2 data
    '''

main()
plt.show()
