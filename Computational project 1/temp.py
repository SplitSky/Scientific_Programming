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
    '''
    variable name               type                    description
    temp                        list                    stores the data entries
    data                        numpy array             stores the complete set of data
    '''
    data = np.genfromtxt(filename, delimiter=',')
    temp = []
    for counter in range(0, len(data), 1):
        if data[counter] == np.NaN:
            temp.append(counter)  # appends the index values that contain invalid entries

    for entry in temp:  # deletes the data that is not valid
        data = np.delete(data, entry)

    return data


def plotData(title, xAxisTitle, yAxisTitle, x, y, error_y, label):
    '''
    variable name               type                    description
    figure                      object                  stores the plot object
    axes_1                      object                  stores the subplot object (child of figure)
    axes_2                      object                  stores the subplot object (child of figure)
    y_weights                   numpy array             stores the weights for the residual plot
    y_errors                    numpy array             stores the array of errors for the plot considered
    fit_parameters              list                    stores the fitting parameters for the plot
    fit_errors                  list                    stores the covariance matrix of the errors
    y_fitted                    numpy array             stores the fitted y-value used for plotting
    title                       string                  the title of the plot
    xAxisTitle                  string                  x-axis title
    yAxisTitle                  string                  y-axis title
    x                           numpy array             the x-value of the plot
    y                           numpy array             the y-values of the plot
    error_y                     numpy array             the errors of the y-values
    '''
    figure = plt.figure()
    axes_1 = figure.add_subplot(121)
    axes_1.plot(x, y, "b+", label=label)
    axes_1.errorbar(x, y, error_y, fmt="b+")
    plt.xlabel(xAxisTitle)
    plt.ylabel(yAxisTitle)
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
    '''
    variable name               type                    description
    figure                      object                  stores the plot object
    axes_1                      object                  stores the subplot object (child of figure)
    axes_2                      object                  stores the subplot object (child of figure)
    y_weights                   numpy array             stores the weights for the residual plot
    y_errors                    numpy array             stores the array of errors for the plot considered
    fit_parameters              list                    stores the fitting parameters for the plot
    fit_errors                  list                    stores the covariance matrix of the errors
    y_fitted                    numpy array             stores the fitted y-value used for plotting
    title                       string                  the title of the plot
    xAxisTitle                  string                  x-axis title
    yAxisTitle                  string                  y-axis title
    x                           numpy array             the x-value of the plot
    y                           numpy array             the y-values of the plot
    error_y                     numpy array             the errors of the y-values
    '''
    figure = plt.figure()
    axes_1 = figure.add_subplot(121)
    axes_1.plot(x, y, "b+", label=label)
    axes_1.errorbar(x, y, error_y)
    plt.ylabel(yAxisTitle)
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
    '''
    variable name               type                    description
    allArrays                   list                    stores the sliced arrays
    channelNumber               numpy array             stores the x-values of the data
    array                       numpy array             stores the complete set of data
    marker                      boolean                 determines which slicing to use
    '''
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

    # slices the plots
    return allArrays


def search(arr, x):
    '''
    variable name               type                    description
    i                           integer                 counter
    arr                         list                    array
    x                           float                   the value being searched
    '''
    # linear search function
    for i in range(len(arr)):
        if arr[i] == x:
            return i
    return -1


def findZero(x, array):
    '''
    variable name               type                    description
    temp2                       numpy array             stores the data set
    index1                      int                     stores the index of the maximum
    index2                      int                     stores the index of the minimum
    cut_array                   list                    stores the sliced array between the peaks
    x_temp                      list                    stores the corresponding sliced x-coordinates
    temp                        list                    stores the x-values for interpolation
    interpolated_y              list                    stores the interpolated values of y
    fitted_parameters           list                    stores the parameters from the linear fitting
    zero                        float                   stores the value at which the curve crosses the 0V point
    sigma                       float                   stores the width of the peak
    '''

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
    sigma = np.abs(array[index2] - array[index1])

    return int(np.rint(zero)), sigma  # returns channel number and difference for temperature calculation


def linearFit(x, y, ey):
    '''
    variable name               type                    description
    x                           numpy array             stores the x values of the data
    y                           numpy array             stores the y values of the data
    ey                          numpy array             stores the errors on the y values
    fit_parameters              list                    stores the fitting parameters for the plot
    fit_errors                  list                    stores the covariance matrix of the errors
    variance_m                  float                   stores the variance of the plot gradient
    variance_c                  float                   stores the variance of the plot intercept
    sigma_m                     float                   stores the error on the gradient
    sigma_c                     float                   stores th error on the intercept
    '''
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
    '''
    variable name               type                    description
    x                           numpy array             stores the x values of the data
    y                           numpy array             stores the y values of the data
    ey                          numpy array             stores the errors on the y values
    fit_parameters              list                    stores the fitting parameters for the plot
    fit_errors                  list                    stores the covariance matrix of the errors
    variance_a                  float                   stores the variance of the quadratic term
    variance_b                  float                   stores the variance of the linear term
    variance_c                  float                   stores the variance of the constant term
    sigma_a                     float                   stores the error on the quadratic term
    sigma_b                     float                   stores th error on the linear term
    sigma_c                     float                   the error on the constant term
    fit_a                       float                   stores the value of the quadratic  term
    fit_b                       float                   stores the value of the linear term
    fit_c                       float                   stores the value of the constant term
    '''
    # Perform a quadratic fit and get the errors

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
    '''
    variable name               type                    description
    :param channel_number:      numpy array             stores the channel numbers
    :param data:                numpy array             stores the complete set of data
    :param marker:              boolean                 used for splitting arrays
    zeroes                      list                    stores the points at which the data crosses the 5V line
    differences                 list                    stores the peak widths
    temp                        float/list              temporary variable
    temp2                       float                   stores the peak width
    '''
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
    '''
    variable name               type                    description
    fit_y                       numpy array             stores the fitted y-values
    y                           numpy array             stores the y values of the data
    ey                          numpy array             stores the errors on the y values

    '''
    # all arrays are numpy arrays
    # returns the chi squared value
    chi_sqrt = ((y - fit_y) / (ey)) ** 2
    return np.sum(chi_sqrt)


def main():
    '''
    variable name               type                    description
    h2o_data                    numpy array             holds the water data
    channel_number              numpy array             holds the channel numbers
    wave_number_water           numpy array             wave numbers for water
    water_zeroes                numpy array             water 5V intersection points
    water_zeroes_sigma          numpy array             the peak widths for water
    fit_param_cal               numpy array             fit parameters of the calibration fit
    chi_sqrt                    float                   chi squared
    c2h2_data                   numpy array             the ethyne data
    c2h2_zeroes                 numpy array             the ethyne 5V intersection points
    m                           numpy array             quantum number j
    results                     numpy array             linear fit parameters of ethyne
    error_arr                   numpy array             the error array
    results_2                   numpy array             quadratic fit of ethyne
    B                           float                   fitting parameter
    B_sigma                     float                   fitting parameter error
    C                           float                   fitting parameter
    C_sigma                     float                   fitting parameter error
    F                           float                   calibration constant
    F_sigma                     float                   calibration constant error
    I0                          float                   moment of inertia
    I1                          float                   moment of inertia
    mom_in_err0                 float                   moment of inertia error
    mom_in_err1                 float                   moment of inertia error
    delta_I                     float                   the change in the moments of inertia
    delta_I_err                 float                   the change in the moments of inertia error
    temp                        float                   temperatures
    temperature_error           float                   the error in temperatures
    '''

    h2o_data = getData(water_file)
    channel_number = np.arange(1, len(h2o_data) + 1)  # get the x-coordinates

    wave_number_water = [12683.782, 12685.540, 12685.769, 12687.066]  # y-coordinate for calibration
    wave_number_water = np.array(wave_number_water)
    wave_number_water = wave_number_water * 100  # unit conversion to SI

    water_zeroes, water_zeroes_sigma = getZeroes(channel_number, h2o_data, True)

    print("calibration fitting for the wave number against channel number plot")
    fit_param_cal = linearFit(wave_number_water, water_zeroes, np.ones(np.size(water_zeroes)))
    chi_sqrt = getChiSqrt(fit_param_cal[0][0] * wave_number_water + fit_param_cal[1][0], water_zeroes,
                          np.ones(np.size(water_zeroes)))
    print("The value of chi squared is: {:0.9}".format(chi_sqrt))
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

    print("Quadratic fitting for the c2h2 zeroes data against quantum number")
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
    F = 1 / fit_param_cal[0][0]

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
    temp = np.mean(temp)
    print("The temperature of water is: {0:.09} +/-  ".format(temp) + str(temperature_error) + " K")

    # c2h2 temperature
    temp = np.sqrt(m_e * c ** 2 / k) * (c2h2_zeroes_sigma / c2h2_zeroes)
    temperature_error = np.abs(np.std(temp))
    temp = np.mean(temp)
    print("The temperature of c2h2 is: {0:0.09} +/-  ".format(temp) + str(temperature_error) + " K")

    # plots
    plotData("C2H2 data linear fit", "Quantum number m", "Channel number", m, c2h2_zeroes, error_arr, "C2H2")
    plotQuadraticFitData("C2H2 data polynomial fit", "Quantum number m", "Channel number", m, c2h2_zeroes, error_arr,
                         "C2H2")

    error_arr = np.delete(error_arr, 0)
    plotData("Water data", "Wave number / $m^{-1}$", "Channel number", wave_number_water, water_zeroes, error_arr,
             "H2O")

    # finish the lines by adding or removing different stuff
    # make sure the graphs are consistent


main()
plt.show()
