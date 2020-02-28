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
figure = plt.figure()

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


def linearFit(x, y, ey):
    #plt.plot(x, y)  #######################################################################################
    # Perform a linear fit and get the errors
    fit_parameters, fit_errors = np.polyfit(x, y, 1, cov=True, w= ey)
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

    # error = np.std(np.abs(y - (fit_m * x + fit_c)))  # returns the estimated error for the value # incorrect but still may be useful
    return [[fit_m, sigma_m], [fit_c, sigma_c]]


def quadraticFit(x, y, ey):
    # Perform a quadratic fit and get the errors
    print(x)
    print(y)
    print(ey)

    plt.plot(x,y,"b+")


    fit_parameters, fit_errors = np.polyfit(x, y, 2, cov=True, w= ey)
    fit_a = fit_parameters[0] # quadratuc term
    fit_b = fit_parameters[1] # linear term
    fit_c = fit_parameters[2] # constant term

    plt.plot(x, fit_a*x**2 + fit_b*x + fit_c)
    print(fit_errors)

    variance_a = fit_errors[0][0]
    variance_b = fit_errors[1][1]
    variance_c = fit_errors[2][2]
    sigma_a = np.sqrt(variance_a)
    sigma_b = np.sqrt(variance_b)
    sigma_c = np.sqrt(variance_c)
    print('Quadratic fit of y = a*x**2 + b*x + c')
    print('Quadratic term  A = {:04.10f} +/- {:04.10f}'.format(fit_a, sigma_a))
    print('Linear term     B = {:04.10f} +/- {:04.10f}'.format(fit_b, sigma_b))
    print('Intercept       C = {:04.10f} +/- {:04.10f}'.format(fit_c, sigma_c))
    return [[fit_a, sigma_a], [fit_b, sigma_b], [fit_c, sigma_c]]


def function(m, c, x):
    return m * x + c


def addErrors(type, a, error_a, b, error_b, c):
    # both functions return the absolute error
    if type:  # this type is adding in quadrature
        return np.sqrt(error_a ** 2 + error_b ** 2)
    else:  # this is for adding percentage error
        return c * np.sqrt((error_b / b) ** 2 + (error_a / a) ** 2)


def getZeroes(channel_number, data, marker):
    peaks_arrays = splitArrays(channel_number, data, marker)  # data arrays for the peaks
    zeroes = []  # the points where the graph crosses the 5V point
    # uses the maximum and minimum values to find the zero points
    for peak in peaks_arrays:
        temp = findZero(peak[0], peak[1])
        zeroes.append([temp, data[temp - 1]])

    temp = []
    for counter in range(0, len(zeroes)):
        temp.append(zeroes[counter][0])


    return np.array(temp)


def plot_residuals(x, y, err_y):
    y_sigma = err_y
    # Create array of weights (1/sigma) for y values, with same size as data array y
    y_weights = (1 / y_sigma) * np.ones(np.size(y))
    # Create array of error values for the error bar plot below; each element is y_sigma
    y_errors = y_sigma * np.ones(np.size(y))
    # Perform fit using np.polyfit
    fit_parameters, fit_errors = np.polyfit(x, y, 1, cov=True, w=y_weights)
    # Create set of fitted y values using fit parameters from np.polyfit, and original x values
    y_fitted = np.polyval(fit_parameters, x)
    # Make plot of the residuals, using the 'errorbar' plotting

    axes_1 = figure.add_subplot(111)
    axes_1.errorbar(x, y - y_fitted, yerr=y_errors, fmt='oy')

def getChiSqrt(fit_y,y,ey):
    # all arrays are numpy arrays
    # returns the chi squared value
    chi_sqrt = ( (y - fit_y) / (ey))**2
    return np.sum(chi_sqrt)

def main():

    ##          error = 1
    h2o_data = getData(water_file)
    channel_number = np.arange(1, len(h2o_data) + 1)  # get the x-coordinates
    # plt.plot(channel_number, h2o_data)

    wave_number_water = [12683.782, 12685.540, 12685.769, 12687.066]
    wave_number_water = np.array(wave_number_water)
    wave_number_water = wave_number_water*100 # unit conversion

    water_zeroes = getZeroes(channel_number, h2o_data, True)

    print("calibration fitting for the wave number against channel number plot")
    fit_param_cal = linearFit(water_zeroes, wave_number_water,np.ones(np.size(water_zeroes)))
    print(" ")
    '''
    enter the channel number and get the wave number
    '''
    c2h2_data = getData(experiment_file)
    channel_number = np.arange(1, len(c2h2_data) + 1)  # get the x-coordinates
    c2h2_zeroes = getZeroes(channel_number, c2h2_data, False)

    # values of the wave number from calibrated spectrum
    # linear fitting of the c2h2 data #
    error_arr = np.ones(np.size(c2h2_zeroes))
    m = np.array([3, 4, 5, 6, 7])
    print("Linear fitting for the c2h2 zeroes data against quantum number")
    results = linearFit(m, c2h2_zeroes,error_arr)
    chi_sqrt = getChiSqrt(m*results[0][0]+results[1][0],c2h2_zeroes,error_arr)
    print("The value of chi squared is: {:04.10f}".format(chi_sqrt))
    print(" ")

    print("Quadratic fitting for the c2h2 zeroes data against quatum number")
    error_arr = np.ones(np.size(c2h2_zeroes)) # re declare remove later
    results_2 = quadraticFit(m, c2h2_zeroes,error_arr)
    chi_sqrt = getChiSqrt(results_2[0][0]*m**2 + results_2[1][0]*m + results_2[2][0], c2h2_zeroes, error_arr)
    print("The value of chi squared is: {:04.10f}".format(chi_sqrt))
    print(" ")

    A = results_2[2][0]
    A_sigma = results_2[2][1]
    B = results_2[1][0]
    B_sigma = results_2[1][1]
    C = results_2[0][0]
    C_sigma = results_2[0][1]
    F_sigma = fit_param_cal[0][1]
    F = fit_param_cal[0][0]

    I0 = (h/ (2*np.pi))**2 / ((F*h*c)*(B-C))
    I1 = (h/ (2*np.pi))**2 / ((F*h*c)*(B+C))

    mom_in_err0 = np.sqrt(I0**2 * ((F_sigma/F)**2 + (np.sqrt(B_sigma**2 + C_sigma**2)/ (B + C) )**2))
    mom_in_err1 = np.sqrt(I0**2 * ((F_sigma/F)**2 + (np.sqrt(B_sigma**2 + C_sigma**2)/ (B - C) )**2))

    delta_I = np.abs( (h)**2/(F*h*c) * (2*C/B**2) )
    delta_I_err = delta_I **2 * ( (F_sigma/F)**2 + (2*B_sigma/B)**2 + (C_sigma/C)**2 )
    delta_I_err = np.sqrt(delta_I_err)

    print("moment of inertia I_0 is: {:0.9} +/- {:0.9}".format(I0, mom_in_err0))
    print("moment of inertia I_1 is: {:0.9} +/- {:0.9}".format(I1, mom_in_err1))
    print("the difference in moment of inertia is: {:0.9} +/- {:0.9}".format(delta_I, delta_I_err))





main()
plt.show()
