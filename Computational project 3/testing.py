# -*- coding: utf-8 -*-
"""
Author: Tomasz Neska
Date: 13/05/2020
Description: Project 3 - Monte Carlo techniques - Penetration of neutrons through shielding
"""
# initialisation
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

plt.rcParams.update({'font.size': 14})
plt.style.use('default')
figure = plt.figure()
plt.rcParams.update({'errorbar.capsize': 2})


def linearFit(x, flux, meanFreePath, plot, ey):
    '''
    :param x:
    :param flux:
    :param meanFreePath:
    :param plot:
    :return: gradient and the error

    variable name               type                    description
    temp                        float                   stores the temporary values
    flux                        float                   stores the parameters being evaluated
    x                           numpy array             stores the x-axis values for plotting
    y                           numpy array             stores the y-axis values for plotting
    fit_m                       float                   stores the gradient
    fit_c                       float                   stores the intercept
    variance_m                  float                   stores the variance of the gradient
    variance_c                  float                   stores the variance of the intercept
    sigma_m                     float                   stores the standard error of the gradient
    sigma_c                     float                   stores the standard error of the intercept
    figure                      object                  stores the figure object
    axes                        object                  stores the subplot object
    '''

    temp = []
    for counter in range(0, len(flux), 1):
        if (flux[counter] == 0):
            temp.append(counter)
    flux = np.delete(flux, temp)
    x = np.delete(x, temp)
    y = np.log(flux)

    y_weights = (1 / error_y) * np.ones(np.size(y)) #########################################

    fit_parameters, fit_errors = np.polyfit(x, y, 1, w=ey, cov=True)
    fit_m = fit_parameters[0]
    fit_c = fit_parameters[1]
    variance_m = fit_errors[0][0]
    variance_c = fit_errors[1][1]
    sigma_m = np.sqrt(variance_m)
    sigma_c = np.sqrt(variance_c)

    if plot:
        figure = plt.figure()
        axes = figure.add_subplot(111)
        axes.plot(x, fit_m * x + fit_c)
        axes.scatter(x, y)
        axes.set_xlabel("Number")
        axes.set_ylabel("Frequency")

        print('Linear np.polyfit of y = m*x + c')
        print('Gradient  m = {:04.10f} +/- {:04.10f}'.format(fit_m, sigma_m))
        print('Intercept c = {:04.10f} +/- {:04.10f}'.format(fit_c, sigma_c))
        print("mean free path -1/lambda =" + str(-1 / meanFreePath))

    return fit_m, sigma_m


class LCG:  # taken directly from the course material
    """
    A general linear congruential generator
    """

    def __init__(self, m, a, c):
        self.m = m
        self.a = a
        self.c = c
        self.seed = 0
        self.this_sample = self.seed  # Set initial sequence value to be the seed
        # Can return the original seed value if we want to!

    def sample(self):
        # Generate the sample (between 0 and m)
        self.this_sample = (self.a * self.this_sample + self.c) % self.m
        # Return the sample (between 0 and 1)
        return self.this_sample / self.m

    # Allow the seed value to be set explicitly
    def set_seed(self, seed_val):
        self.seed = seed_val
        self.this_sample = self.seed


class Random_Generator(object):
    def __init__(self):
        '''
        This is the initialisation function
            variable name               type                    description
            self.array                  list                    stores the values obtained from random number generator methods
            self.x                      list                    stores x-values for plotting
            self.expDist                nested list             stores the exponential distibution
            self.length                 numpy array             stores the length of an evaluated vector
            self.vectors                float                   store the complete set of vectors
        '''
        self.array = []
        self.x = []
        self.expDist = []
        self.length = []
        self.vectors = []

    def badGenerateArray(self, N):
        # generates numbers between 0 and 1.

        random_set = np.array([])
        # This is a pretty terrible way to generate 1000 random numbers, but it explicitly shows that we can
        # generate 1 random number at a time and put 1000 of them into an array
        for i in range(N):
            random_set = np.append(random_set, np.random.random())
        self.array = random_set

    def weightedAverage(self, y, ey):
        y = np.array(y)
        ey = np.array(ey)
        weighted_mean = np.sum(y / ey ** 2) / np.sum(1 / ey ** 2)
        error = np.sqrt(1 / np.sum(1 / ey ** 2))
        return weighted_mean, error

    def randssp(self, p, q):
        # q - how many sets of random numbers
        # p - how many numbers in each set

        global m, a, c, x

        try:
            x
        except NameError:
            m = pow(2, 31)
            a = pow(2, 16) + 3
            c = 0
            x = 123456789

        try:
            p
        except NameError:
            p = 1
        try:
            q
        except NameError:
            q = p

        r = np.zeros([p, q])

        for l in range(0, q):
            for k in range(0, p):
                x = np.mod(a * x + c, m)
                r[k, l] = x / m

        return r

    def betterGenerateArray(self, N):
        random_set = np.random.uniform(0, 1, N)
        self.array = random_set

    def testRandomGenerator(self, N):
        # plots the distribution of the random numbers generated
        figure = plt.figure()
        axes = figure.add_subplot(131)
        axes2 = figure.add_subplot(132)
        axes3 = figure.add_subplot(133)

        bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        # generate numbers bad
        self.badGenerateArray(N)
        bin_frequencies, bin_locations = np.histogram(self.array, bins=bins)
        axes.hist(self.array, bins=bins)
        axes.set_xlabel("Number")
        axes.set_ylabel("Frequency")
        mean = np.mean(bin_frequencies)
        std = np.std(bin_frequencies)

        # generate numbers good
        self.betterGenerateArray(N)
        bin_frequencies, bin_locations = np.histogram(self.array, bins=bins)
        axes2.hist(self.array, bins=bins)
        axes2.set_xlabel("Number")
        axes2.set_ylabel("Frequency")
        mean2 = np.mean(bin_frequencies)
        std2 = np.std(bin_frequencies)

        # pseudo random
        temp = self.randssp(1, N)
        temp = temp[0]
        bin_frequencies, bin_locations = np.histogram(self.array, bins=bins)
        axes3.hist(temp, bins=bins)
        axes3.set_xlabel("Number")
        axes3.set_ylabel("Frequency")
        mean3 = np.mean(bin_frequencies)
        std3 = np.std(bin_frequencies)

        # compare number generators
        # binomial calculation
        probablity = 1 / (len(bins) - 1)
        binomial_mean = N * probablity
        binomial_std = np.sqrt(binomial_mean * (1 - binomial_mean / N))
        print("The expected value from the binomial distribution is: {:0.9} with fluctuation of  {:0.9}".format(
            binomial_mean, binomial_std))
        print("The first generator gave a mean of {:0.9} +/- {:0.9}".format(mean, std))
        print("The second generator gave a mean of {:0.9} +/- {:0.9}".format(mean2, std2))
        print("The pseudo random generator gave a mean of {:0.9} +/- {:0.9}".format(mean3, std3))

    def spectraPhen(self):
        # Samples using the inbuilt generator
        nsamples = 2000
        x_np = np.zeros(nsamples)
        y_np = np.zeros(nsamples)
        z_np = np.zeros(nsamples)
        for i in range(nsamples):
            x_np[i] = 2 * np.random.random() - 1
            y_np[i] = 2 * np.random.random() - 1
            z_np[i] = 2 * np.random.random() - 1

        fig = plt.figure()
        ax1 = fig.add_subplot(111, projection='3d')
        ax1.scatter(x_np, y_np, z_np)

        # Makes samples using RANDU algorithm
        randu = LCG(2 ** 31, 65539, 0)
        nsamples = 5000
        x_np = np.zeros(nsamples)
        y_np = np.zeros(nsamples)
        z_np = np.zeros(nsamples)
        randu.set_seed(28538)
        for i in range(nsamples):
            x_np[i] = 2 * randu.sample() - 1
            y_np[i] = 2 * randu.sample() - 1
            z_np[i] = 2 * randu.sample() - 1

        fig = plt.figure()
        ax1 = fig.add_subplot(111, projection='3d')
        ax1.scatter(x_np, y_np, z_np)

    def plotDistribution3D(self, x, y, z):
        fig = plt.figure()
        ax1 = fig.add_subplot(111, projection='3d')
        ax1.scatter(x, y, z)

    def quickHist(self, data, step):
        bins = np.arange(0, 1, step)
        figure = plt.figure()
        axes = figure.add_subplot(111)
        axes.set_ylabel("Frequency")
        axes.set_xlabel("Number")
        axes.hist(data, bins=bins)

    def generate_vectors(self, accuracy, plot, r=1):
        # generates isotropic vectors distributed on a surface of a sphere
        # range of u is -1 to 1
        # range of theta is 0 to 2 pi
        self.betterGenerateArray(accuracy)
        u = 2 * self.array - 1
        self.betterGenerateArray(accuracy)
        theta = 2 * np.pi * self.array

        x = r * np.cos(theta) * np.sqrt(1 - u ** 2)
        y = r * np.sin(theta) * np.sqrt(1 - u ** 2)
        z = r * u
        self.length = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        if plot:
            self.plotDistribution3D(x, y, z)
        self.vectors = [x, y, z]

    def generateExponentialDist(self, N, plot, meanFreePath):
        step = 0.1
        bins = np.arange(0, 1, step)
        self.betterGenerateArray(N)
        self.expDist = -1 * meanFreePath * np.log(self.array)
        histogram = np.histogram(self.expDist)
        y = histogram[0]
        x = histogram[1][:len(histogram[1]) - 1] + step / 2
        y = np.log(y)
        # eliminate zeroes
        temp = []
        for counter in range(0, len(y), 1):
            if y[counter] == -np.inf:
                temp.append(counter)
        y = np.delete(y, temp)
        x = np.delete(x, temp)

        fit_parameters, fit_errors = np.polyfit(x, y, 1, cov=True)
        fit_m = fit_parameters[0]
        fit_c = fit_parameters[1]
        variance_m = fit_errors[0][0]
        variance_c = fit_errors[1][1]
        sigma_m = np.sqrt(variance_m)
        sigma_c = np.sqrt(variance_c)

        if plot:
            figure = plt.figure()
            axes = figure.add_subplot(111)
            axes.plot(x, fit_m * x + fit_c)
            axes.scatter(x, y)
            axes.set_xlabel("Number")
            axes.set_ylabel("Frequency")

            print('Linear np.polyfit of y = m*x + c')
            print('Gradient  m = {:04.10f} +/- {:04.10f}'.format(fit_m, sigma_m))
            print('Intercept c = {:04.10f} +/- {:04.10f}'.format(fit_c, sigma_c))
            print("mean free path -1/lambda =" + str(-1 / meanFreePath))

        self.expDist = [x, y]

    def generateDistVectors(self, meanFreePath, N):
        r = np.random.uniform(0, 1, N)
        r = -1 * meanFreePath * np.log(r)  # makes the lengths exponentially distributed

        number = np.random.uniform(0, 1, N)

        u = 2 * number - 1

        number = np.random.uniform(0, 1, N)

        theta = 2 * np.pi * number

        x = r * np.cos(theta) * np.sqrt(1 - u ** 2)
        y = r * np.sin(theta) * np.sqrt(1 - u ** 2)
        z = r * u

        self.vectors = [x, y, z]
        self.plotDistribution3D(x, y, z)

    def genDistVector(self, meanFreePath):
        r = np.random.uniform(0, 1, 1)
        r = -1 * meanFreePath * np.log(r)  # makes the lengths exponentially distributed

        self.betterGenerateArray(1)
        u = 2 * self.array - 1
        self.betterGenerateArray(1)
        theta = 2 * np.pi * self.array

        x = r * np.cos(theta) * np.sqrt(1 - u ** 2)
        y = r * np.sin(theta) * np.sqrt(1 - u ** 2)
        z = r * u

        return [x.tolist()[0], y.tolist()[0], z.tolist()[0]]

    def getLength(self, vector):
        return np.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)

    def generateDirection(self, r, vector):
        length = self.getLength(vector)  # normalizes the vector -> length 1
        vector[0] = vector[0] / length
        vector[1] = vector[1] / length
        vector[2] = vector[2] / length

        for counter in range(0, 3, 1, ):  # extends the vector appropriately
            vector[counter] *= r

        return vector


class Experiment(Random_Generator):
    def __init__(self, material1, N, thickness, thickness2, name, material2):
        self.data = []
        self.absArea = material1[0]
        self.scatterArea = material1[1]
        self.density = material1[2]  # number density
        self.thickness = thickness
        self.thickness2 = thickness2
        self.N = N
        self.meanFreePath = 1 / (self.density * self.absArea + self.density * self.scatterArea)
        print(self.meanFreePath)
        print("meme")
        self.history = []
        self.gradients = [[], [], []]  # stores gradients form many experiments
        # absorbed, reflected, transmitted
        self.name = name
        self.material2 = material2
        self.abs_err = []
        self.pass_err = []
        self.ref_err = []

        super().__init__()

    def randomWalk(self, T):
        vector = self.genDistVector(self.meanFreePath)
        # T -> thickness
        SigmaT = 0
        SigmaT += self.density * self.absArea
        SigmaT += self.density * self.scatterArea
        self.meanFreePath = 1 / SigmaT

        prob_absorption = self.density * self.absArea / (self.density * self.absArea + self.density * self.scatterArea)
        is_absorbed = 0
        i = 0
        x = 0
        self.history = [[], [], []]

        while is_absorbed == 0:
            if i == 0:
                vector[0] += self.getLength(vector)  ## adds the vectors
            else:
                vector = self.addVectors(self.genDistVector(self.meanFreePath), vector)

            self.history[0].append(vector[0])
            self.history[1].append(vector[1])
            self.history[2].append(vector[2])
            x = vector[0]
            if (x < 0):
                return "reflected", self.history
            elif (x > T):
                return "passed", self.history

            # evaluate what happens to the particle
            probability = np.random.uniform(0, 1, 1)[0]
            if probability <= prob_absorption:
                # gets absorbed
                is_absorbed = 1
                return "absorbed", self.history
            else:
                # scatters
                i += 1

    def compressData(self, array):
        # compresses experiment data
        for entry in array:
            absorbed = entry[0]
            passed = entry[1]
            reflected = entry[2]

        self.abs_err.append(np.std(absorbed))
        self.pass_err.append(np.std(passed))
        self.ref_err.append(np.std(reflected))

        return [np.mean(absorbed), np.mean(passed), np.mean(reflected)]

    def addVectors(self, vector1, vector2):
        temp = [0, 0, 0]
        for counter in range(0, 3, 1):
            temp[counter] = vector1[counter] + vector2[counter]

        return temp

    def plotRandomWalk(self):
        x = self.history[0]
        y = self.history[1]
        z = self.history[2]
        fig = plt.figure()
        ax1 = fig.add_subplot(111, projection='3d')
        ax1.plot(x, y, z)

    def experiment(self, N, thickness):
        # performs an experiment N times
        histories = []
        results = [0, 0, 0]
        for i in range(0, N, 1):
            result, temp = self.randomWalk(thickness)
            histories.append(result)

        for entry in histories:
            # counts the outcomes
            if entry == "absorbed":
                results[0] += 1
            elif entry == "passed":
                results[1] += 1
            elif entry == "reflected":
                results[2] += 1  #
        self.data = results

    def experimentWoodcock(self, N, T1, T2):
        # performs an experiment N times
        histories = []
        results = [0, 0, 0]
        for i in range(0, N, 1):
            result, temp = self.woodcockMethod()
            histories.append(result)

        for entry in histories:
            # counts the outcomes
            if entry == "absorbed":
                results[0] += 1
            elif entry == "passed":
                results[1] += 1
            elif entry == "reflected":
                results[2] += 1  #
        self.data = results

    def thicknessPlot(self, N, min_T, max_T, step):
        results = []
        thickness = np.arange(min_T, max_T, step=step)
        for entry in thickness:
            temp = []
            for i in range(10):
                self.experiment(N, entry)
                temp.append(self.data)
                print(temp)

            temp = self.compressData(temp)
            results.append(temp)

        self.data = results
        self.thickness = thickness  # array storing thickness

        absorbed = []
        transmitted = []
        reflected = []

        for entry in self.data:
            absorbed.append(entry[0])
            transmitted.append(entry[1])
            reflected.append(entry[2])

        figure = plt.figure()
        # plotting for absorption
        flux = np.array(absorbed)
        marker = [0, 0, 0]
        temp = 0
        for entry in absorbed:
            if entry == 0:
                temp += 1

        if not temp == len(absorbed):

            # printing area
            print("start")
            print(len(flux))
            print(len(thickness))
            print(len(self.ref_err))
            print(len(self.ref_err))
            print(len(self.pass_err))

            print("end")

            self.gradients[0].append(linearFit(thickness, flux, self.meanFreePath, False, self.abs_err))
            ax1 = figure.add_subplot(131)
            ax1.scatter(thickness, flux, label="Absorption")
            ax1.plot(thickness, flux)
            ax1.errorbar(thickness, flux, yerr=self.abs_err, fmt='b+')
        else:
            marker[0] = 1

        temp = 0
        for entry in reflected:
            if entry == 0:
                temp += 1

        if not temp == len(reflected):
            # plotting for reflection
            flux = np.array(reflected)



            self.gradients[1].append(linearFit(thickness, flux, self.meanFreePath, False, self.ref_err))
            ax2 = figure.add_subplot(132)
            ax2.scatter(thickness, flux, label="Reflection")
            ax2.errorbar(thickness, flux, yerr=self.ref_err, fmt='b+')
        else:
            marker[1] = 1

        temp = 0
        for entry in transmitted:
            if entry == 0:
                temp += 1
        if not temp == len(transmitted):
            # plotting for transmission
            flux = np.array(transmitted)
            self.gradients[2].append(linearFit(thickness, flux, self.meanFreePath, False, self.pass_err))
            ax3 = figure.add_subplot(133)
            ax3.scatter(thickness, flux, label="Transmission")
            ax3.plot(thickness, flux)
            ax3.errorbar(thickness, flux, yerr=self.ref_err, fmt='b+')
        else:
            marker[2] = 1

        if np.array(marker).sum() != 3:
            figure.legend()
            print(self.gradients)
            temp = -1 * 1 / self.gradients[2][0][0]
            temp2 = self.gradients[2][0][1] / self.gradients[2][0][0]
            temp2 *= temp
            temp2 = np.abs(temp2)
            print(
                "The attenuation length for " + str(self.name) + " is " + str(temp) + str(" +/- ") + str(temp2) + " cm")
        else:
            print("no data obtained")

    def thing(self):
        x = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        i = 0
        for item in self.gradients[0]:
            print(item[0])
            plt.plot(x, item[0] * x, label=str(i))
            i += 1
            plt.legend()

    def percentageAbsorption(self, N):
        # check variation of error with number of neutrons used
        thickness = 10  # cm
        results = []
        for counter in range(50):
            self.experiment(N, thickness)
            results.append(self.data)

        absorbed = []
        transmitted = []
        reflected = []

        for entry in results:
            absorbed.append(entry[0])
            transmitted.append(entry[1])
            reflected.append(entry[2])

        mean_absorbed = np.mean(absorbed) / N * 100
        mean_transmitted = np.mean(transmitted) / N * 100
        mean_reflected = np.mean(reflected) / N * 100

        std_abs = np.std(absorbed) / np.mean(absorbed) * mean_absorbed
        std_tra = np.std(transmitted) / np.mean(transmitted) * mean_transmitted
        std_ref = np.std(reflected) / np.mean(reflected) * mean_reflected

        print("----------------------------------------------------------------------------------------")
        print("Transmission through a fixed Thickness - " + str(self.name))
        print("Thickness : 10 cm ")
        print("Total Neutrons: " + str(N))
        print("Neutrons Reflected: " + str(np.mean(reflected)))
        print("Neutrons Transmitted: " + str(np.mean(transmitted)))
        print("Neutrons Absorbed: " + str(np.mean(absorbed)))
        print("Percentage Transmitted: " + str(mean_transmitted))

        print(
            'The average percentage of neutrons that were absorbed is: {:04.10f} % +/- {:04.10f}'.format(mean_absorbed,
                                                                                                         std_abs))
        print('The average percentage of neutrons that were transmitted is: {:04.10f} % +/- {:04.10f}'.format(
            mean_transmitted, std_tra))
        print('The average percentage of neutrons that were reflected is: {:04.10f} % +/- {:04.10f}'.format(
            mean_reflected, std_ref))

    def woodcockMethod(self):
        # Material 1 is always the first block encountered
        T1 = self.thickness
        T2 = self.thickness2 + T1
        marker = 0
        Sigma1 = 1 / self.meanFreePath
        Sigma2 = self.material2[0] * self.material2[2] + self.material2[1] * self.material2[2]
        SigmaT = 0
        if Sigma1 > Sigma2:
            marker = 2
            SigmaT = Sigma1
            prob_fictitious = 1 - Sigma2 / SigmaT
        else:
            marker = 1
            SigmaT = Sigma2
            prob_fictitious = 1 - Sigma1 / SigmaT

        self.meanFreePath = 1 / SigmaT

        # T -> thickness
        prob_absorption1 = self.density * self.absArea / (self.density * self.absArea + self.density * self.scatterArea)
        prob_absorption2 = self.material2[0] * self.material2[2] / (
                self.material2[0] * self.material2[2] + self.material2[1] * self.material2[2])

        is_absorbed = 0
        i = 0
        x = 0
        self.history = [[], [], []]

        r = np.random.uniform(0, 1)
        r = -1 * self.meanFreePath * np.log(r)
        vector = [r, 0, 0]

        while is_absorbed == 0:
            print(i)
            print(vector)
            if i == 0:
                vector[0] += self.getLength(vector)
            else:
                vector = self.addVectors(self.genDistVector(self.meanFreePath), vector)

            self.history[0].append(vector[0])
            self.history[1].append(vector[1])
            self.history[2].append(vector[2])
            x = vector[0]
            if (x < 0):
                return "reflected", self.history
            elif (x > T2):
                return "passed", self.history

            # evaluate what happens to the particle
            probability = np.random.uniform(0, 1)

            # check the region the particle is in:
            if (x < T1) and (x > 0):  # it's in region 1
                if marker == 1:
                    if probability > prob_fictitious:
                        r = np.random.uniform(0, 1)
                        r = -1 * self.meanFreePath * np.log(r)
                        direction = self.generateDirection(r, vector)
                        vector = self.addVectors(direction, vector)  # makes the step
                        self.history[0].append(vector[0])
                        self.history[1].append(vector[1])
                        self.history[2].append(vector[2])
                        i += 1
                    else:
                        probability = np.random.uniform(0, 1, 1)[0]
                        if probability >= prob_absorption1:
                            # gets absorbed
                            is_absorbed = 1
                            return "absorbed", self.history
                        else:
                            # scatters
                            i += 1
                else:
                    # normal behaviour in region 1
                    probability = np.random.uniform(0, 1, 1)[0]
                    if probability >= prob_absorption1:
                        is_absorbed = 1
                        return "absorbed", self.history
                    else:  # scatters
                        i += 1

            elif (x < T2) and (x >= T1):  # it's in region 2
                if marker == 2:
                    if probability > prob_fictitious:
                        r = np.random.uniform(0, 1)
                        r = -1 * self.meanFreePath * np.log(r)
                        direction = self.generateDirection(r, vector)
                        vector = self.addVectors(direction, vector)  # makes the step
                        self.history[0].append(vector[0])
                        self.history[1].append(vector[1])
                        self.history[2].append(vector[2])
                        i += 1
                    else:
                        probability = np.random.uniform(0, 1)
                        if probability >= prob_absorption2:
                            is_absorbed = 1
                            return "absorbed", self.history
                        else:
                            # scatters
                            i += 1
                else:
                    # normal behaviour in region 2
                    probability = np.random.uniform(0, 1, 1)[0]
                    if probability >= prob_absorption2:
                        is_absorbed = 1
                        return "absorbed", self.history
                    else:  # scatters
                        i += 1


def main():
    water = [0.6652, 103.0, 1.00]

    water_meanFreePath = 1 / (water[0] * water[2] + water[1] * water[2])

    lead = [0.158E-24, 11.221E-24, 1.985E46]
    graphite = [0.0045E-24, 4.74E-24, 6.828E46]
    water = [0.6652E-24, 103.0E-24, 3.343E22]

    gen1 = Experiment(water, 1000, 10, 10, "water",
                      lead)  # absorbArea, scatterArea, density, N, thickness,thickness2, name, material2
    gen2 = Experiment(lead, 1000, 10, 10, "graphite", graphite)
    gen3 = Experiment(graphite, 1000, 10, 10, "water", water)
    # gen1.experimentWoodcock(1000, 40, 40)
    # demonstrate the generator
    # gen1.testRandomGenerator(10000)
    # gen1.spectraPhen()
    # gen1.generateExponentialDist(10000, True, water_meanFreePath)
    # gen1.generateDistVectors(water_meanFreePath, 10000) # generates a cluster of vectors
    # gen1.generate_vectors(1000, True, 1) # sphere distribution
    # experiment 1
    # 10cm material - water. Only absorption
    # water = [0.6652E-24, 0, 3.343E22]
    # gen1 = Experiment(water, 1000, 10, 10, "water", lead)
    gen1.thicknessPlot(100, 0, 10, 1)
    gen1.plotRandomWalk()


main()
plt.show()
