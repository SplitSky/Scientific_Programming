'''
Description:
Author: Tomasz Neska
Date: 12.02.2020
This program summarises and analises two spectra (water and c2h2). All the plots are saved in the script repository
'''
# initialization
import string
from math import *
import numpy as np
import matplotlib.pyplot as plt
import random

plt.rcParams.update({'font.size': 14})
plt.style.use('default')
figure = plt.figure()
plt.rcParams.update({'errorbar.capsize': 2})

# constants
c = 299792458  # speed of light in ms^-1
h = 6.626E-34  # Planck's constant in Js
k = 1.381E-23  # Boltzman's constant in J K^-1
m_e = 9.109E-31  # mass of an electron
h_bar = h / (2 * np.pi)
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


def plotData(title, xAxisTitle, yAxisTitle, x, y, error_y, label):
    figure = plt.figure()
    axes_1 = figure.add_subplot(121)
    axes_1.plot(x, y, "b+", label=label)
    axes_1.errorbar(x, y, error_y, fmt="b+")
    plt.xlabel(xAxisTitle)  #
    plt.ylabel(yAxisTitle)  # edit from axes
    plt.title(title)
    y_weights = (1 / error_y) * np.ones(np.size(y))
    y_errors = error_y * np.ones(np.size(y))
    fit_parameters, fit_errors = np.polyfit(x, y, 1, cov=True, w=y_weights)

    y_fitted = np.polyval(fit_parameters, x)
    axes_1.plot(x, y_fitted)
    axes_2 = figure.add_subplot(122)
    axes_2.set_xlabel(xAxisTitle)
    axes_2.set_ylabel(yAxisTitle)
    axes_2.errorbar(x, y - y_fitted, yerr=y_errors, fmt='b+')

    plt.savefig(title + ".png")


def plotQuadraticFitData(title, xAxisTitle, yAxisTitle, x, y, error_y, label):
    figure = plt.figure()
    axes_1 = figure.add_subplot(121)
    axes_1.plot(x, y, "b+", label=label)
    axes_1.errorbar(x, y, error_y)
    plt.ylabel(yAxisTitle)  # edit from axes
    plt.xlabel(xAxisTitle)
    plt.title(title)
    y_weights = (1 / error_y) * np.ones(np.size(y))
    y_errors = error_y * np.ones(np.size(y))
    fit_parameters, fit_errors = np.polyfit(x, y, 2, cov=True, w=y_weights)
    y_fitted = np.polyval(fit_parameters, x)
    axes_2 = figure.add_subplot(122)
    axes_1.plot(x, y_fitted)
    axes_2.set_xlabel(xAxisTitle)
    axes_2.set_ylabel(yAxisTitle)
    axes_2.errorbar(x, y - y_fitted, yerr=y_errors, fmt='b+')

    plt.savefig(title + ".png")


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

    temp = np.linspace(x_temp[0], x_temp[len(x_temp) - 1], 100)
    interpolated_y = np.interp(temp, x_temp, cut_array)

    fitted_parameters = np.polyfit(temp, interpolated_y, 1)
    zero = -1 * fitted_parameters[1] / fitted_parameters[0]

    #
    sigma = np.abs(array[index2] - array[index1])
    #

    return int(np.rint(zero)), sigma  # returns channel number and difference for temp calculation


def linearFit(x, y, ey):
    # Perform a linear fit and get the errors
    fit_parameters, fit_errors = np.polyfit(x, y, 1, cov=True, w=ey)
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

    return [[fit_m, sigma_m], [fit_c, sigma_c]]


def quadraticFit(x, y, ey):
    # Perform a quadratic fit and get the errors
    print(x)
    print(y)
    print(ey)

    fit_parameters, fit_errors = np.polyfit(x, y, 2, cov=True, w=ey)
    fit_a = fit_parameters[0]  # quadratic term
    fit_b = fit_parameters[1]  # linear term
    fit_c = fit_parameters[2]  # constant term

    variance_a = fit_errors[0][0]
    variance_b = fit_errors[1][1]
    variance_c = fit_errors[2][2]
    sigma_a = np.sqrt(variance_a)
    sigma_b = np.sqrt(variance_b)
    sigma_c = np.sqrt(variance_c)
    print('Quadratic fit of y = a*x**2 + b*x + c')
    print('Quadratic term  A = {:04.10f} +/- {:0.9}'.format(fit_a, sigma_a))
    print('Linear term     B = {:04.10f} +/- {:0.9}'.format(fit_b, sigma_b))
    print('Intercept       C = {:04.10f} +/- {:0.9}'.format(fit_c, sigma_c))
    return [[fit_a, sigma_a], [fit_b, sigma_b], [fit_c, sigma_c]]


def getZeroes(channel_number, data, marker):
    peaks_arrays = splitArrays(channel_number, data, marker)  # data arrays for the peaks
    zeroes = []  # the points where the graph crosses the 5V point
    differences = []
    # uses the maximum and minimum values to find the zero points

    for peak in peaks_arrays:
        temp, temp2 = findZero(peak[0], peak[1])
        zeroes.append([temp, data[temp - 1]])
        differences.append(temp2)

    temp = []
    for counter in range(0, len(zeroes)):
        temp.append(zeroes[counter][0])

    return np.array(temp), np.array(differences)


def getChiSqrt(fit_y, y, ey):
    # all arrays are numpy arrays
    # returns the chi squared value
    chi_sqrt = ((y - fit_y) / (ey)) ** 2
    return np.sum(chi_sqrt)


def getPeakWidth(x, array):
    temp2 = array - 5
    # finds and returns the channel number at which the section of the array passes a zero

    index2 = search(temp2, np.min(temp2))
    index1 = search(temp2, np.max(temp2))

    cut_array = temp2[index2:index1]
    x_temp = x[index2:index1]

    temp = np.linspace(x_temp[0], x_temp[len(x_temp) - 1], 100)
    interpolated_y = np.interp(temp, x_temp, cut_array)

    fitted_parameters = np.polyfit(temp, interpolated_y, 1)
    zero = -1 * fitted_parameters[1] / fitted_parameters[0]

    return int(np.rint(zero))  # returns channel number and the sigma


def main():
    h2o_data = getData(water_file)
    channel_number = np.arange(1, len(h2o_data) + 1)  # get the x-coordinates

    wave_number_water = [12683.782, 12685.540, 12685.769, 12687.066]
    wave_number_water = np.array(wave_number_water)
    wave_number_water = wave_number_water * 100  # unit conversion

    water_zeroes, water_zeroes_sigma = getZeroes(channel_number, h2o_data, True)

    print("calibration fitting for the wave number against channel number plot")
    fit_param_cal = linearFit(water_zeroes, wave_number_water, np.ones(np.size(water_zeroes)))
    print(" ")
    '''
    enter the channel number and get the wave number
    '''
    c2h2_data = getData(experiment_file)
    channel_number = np.arange(1, len(c2h2_data) + 1)  # get the x-coordinates
    c2h2_zeroes, c2h2_zeroes_sigma = getZeroes(channel_number, c2h2_data, False)

    # values of the wave number from calibrated spectrum
    # linear fitting of the c2h2 data #
    error_arr = np.ones(np.size(c2h2_zeroes))
    m = np.array([3, 4, 5, 6, 7])  # quantum number m
    print("Linear fitting for the c2h2 zeroes data against quantum number")
    results = linearFit(m, c2h2_zeroes, error_arr)
    chi_sqrt = getChiSqrt(m * results[0][0] + results[1][0], c2h2_zeroes, error_arr)
    print("The value of chi squared is: {:0.16}".format(chi_sqrt))
    print(" ")

    print("Quadratic fitting for the c2h2 zeroes data against quatum number")
    error_arr = np.ones(np.size(c2h2_zeroes))  # re declare remove later
    results_2 = quadraticFit(m, c2h2_zeroes, error_arr)
    chi_sqrt = getChiSqrt(results_2[0][0] * m ** 2 + results_2[1][0] * m + results_2[2][0], c2h2_zeroes, error_arr)
    print("The value of chi squared is: {:0.16}".format(chi_sqrt))
    print(" ")

    # this variable redeclaration is made to make the code readable

    B = results_2[1][0]
    B_sigma = results_2[1][1]
    C = results_2[0][0]
    C_sigma = results_2[0][1]
    F_sigma = fit_param_cal[0][1]
    F = fit_param_cal[0][0]

    I0 = (h_bar) ** 2 / ((F * h * c) * (B - C))
    I1 = (h_bar) ** 2 / ((F * h * c) * (B + C))

    mom_in_err0 = np.sqrt(I0 ** 2 * ((F_sigma / F) ** 2 + (np.sqrt(B_sigma ** 2 + C_sigma ** 2) / (B + C)) ** 2))
    mom_in_err1 = np.sqrt(I0 ** 2 * ((F_sigma / F) ** 2 + (np.sqrt(B_sigma ** 2 + C_sigma ** 2) / (B - C)) ** 2))

    delta_I = np.abs(h_bar ** 2 / (F * h * c) * (2 * C / B ** 2))
    delta_I_err = delta_I ** 2 * ((F_sigma / F) ** 2 + (2 * B_sigma / B) ** 2 + (C_sigma / C) ** 2)
    delta_I_err = np.sqrt(delta_I_err)

    print("Moment of inertia I_0 is: {:0.9} +/- {:0.9}".format(I0, mom_in_err0))
    print("Moment of inertia I_1 is: {:0.9} +/- {:0.9}".format(I1, mom_in_err1))
    print("The difference in moment of inertia is: {:0.9} +/- {:0.9}".format(delta_I, delta_I_err))

    # temperature calculations
    # water temperature
    temperature_error = 0

    # c2h2_zeroes, c2h2_zeroes_sigma
    temp = np.sqrt(m_e * c ** 2 / k) * (water_zeroes_sigma / water_zeroes)
    temperature_error = np.abs(np.std(temp))
    print(temperature_error)
    temp = np.mean(temp)
    print("The temperature of water is: {0:.09} +/- ".format(temp) + str(temperature_error))

    # c2h2 temperature
    temp = np.sqrt(m_e * c ** 2 / k) * (c2h2_zeroes_sigma / c2h2_zeroes)
    temperature_error = np.abs(np.std(temp))
    temp = np.mean(temp)
    print("The temperature of c2h2 is: {0:0.09} +/- ".format(temp) + str(temperature_error))

    # plots
    plotData("C2H2 data linear fit", "Quantum number m", "Channel number", m, c2h2_zeroes, error_arr, "C2H2")
    plotQuadraticFitData("C2H2 data polynomial fit", "Quantum number m", "Channel number", m, c2h2_zeroes, error_arr,
                         "C2H2")

    error_arr = np.delete(error_arr, 0)
    plotData("Water data", "Wave number", "Channel number", wave_number_water, water_zeroes, error_arr, "H2O")

    # finish the lines by adding or removing different stuff
    # make sure the graphs are consistent


main()
plt.show()
