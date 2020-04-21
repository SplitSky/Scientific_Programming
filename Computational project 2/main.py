# -*- coding: utf-8 -*-
"""
Author: Tomasz Neska
Date: 06/03/2020
Description: Project 2 - utilising three different methods it evaluates the effects of time step and accuracy on the
behaviour of the simple harmonic oscillator
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

plt.rcParams.update({'font.size': 14})
plt.style.use('default')
figure = plt.figure()
plt.rcParams.update({'errorbar.capsize': 2})


class SHO(object):
    def __init__(self, time_step, max_time, b=0.0, m=0.0, k=0.0, init_x=0.0, init_v=0.0, fileNameSave="data.txt",
                 fileNameLoad="data.txt"):
        self.fileNameSave = fileNameSave
        self.fileNameLoad = fileNameLoad
        self.b = b
        self.m = m
        self.k = k
        self.h = time_step
        self.init_v = init_v
        self.init_x = init_x
        self.no_steps = int(np.rint(max_time / time_step))
        self.natural_angular_frequency = np.sqrt(self.k / self.m)
        if self.b != 0:
            self.gamma = self.b / self.m
            self.quality_factor = self.natural_angular_frequency / self.gamma
        self.analytic_series_pos = []
        self.analytic_series_vel = []
        self.analytic_energy = []
        self.__coefficients = []
        self.solver()
        self.analytic_solution()
        self.data = []
        self.b_critical = 2 * np.sqrt(self.k * self.m)

        self.Euler_data = []
        self.B_Euler_data = []
        self.Verlet_data = []
        self.Euler_Cromer_data = []
        self.analytical_data = []
        self.time = np.array(range(0, self.no_steps, 1)) * self.h
        self.disturbed_Verlet_data = []

    def runSimulation(self):
        self.Euler_integrator()
        self.Better_Euler_integrator()
        self.Verlet_integrator()
        self.Euler_Cromer_integrator()
        print("Simulation has been executed")

    def getCoefficients(self):
        return self.__coefficients

    def Euler_integrator(self):
        position_series = [self.init_x]
        velocity_series = [self.init_v]
        for counter in range(1, self.no_steps, 1):
            v_n = velocity_series[len(velocity_series) - 1]
            x_n = position_series[len(position_series) - 1]
            a_n = (-self.b / self.m) * v_n + (-self.k / self.m) * x_n

            position_series.append(x_n + self.h * v_n)
            velocity_series.append(v_n + self.h * a_n)

        self.Euler_data = [position_series, velocity_series, self.energy_function(position_series, velocity_series)]

    def Better_Euler_integrator(self):
        position_series = [self.init_x]
        velocity_series = [self.init_v]
        for counter in range(1, self.no_steps, 1):
            v_n = velocity_series[len(velocity_series) - 1]
            x_n = position_series[len(position_series) - 1]
            a_0 = (-self.b / self.m) * v_n + (-self.k / self.m) * x_n

            position_series.append(x_n + self.h * v_n + 0.5 * self.h ** 2 * a_0)
            velocity_series.append(v_n + self.h * a_0)

        self.B_Euler_data = [position_series, velocity_series, self.energy_function(position_series, velocity_series)]

    def Euler_Cromer_integrator(self):
        position_series = [self.init_x]
        velocity_series = [self.init_v]
        for counter in range(1, self.no_steps, 1):
            v_n = velocity_series[len(velocity_series) - 1]
            x_n = position_series[len(position_series) - 1]

            temp = v_n - (self.k * self.h / self.m) * x_n  # v_n+1
            velocity_series.append(temp)
            position_series.append(x_n + self.h * temp)

        self.Euler_Cromer_data = [position_series, velocity_series,
                                  self.energy_function(position_series, velocity_series)]
        print(np.array(self.Euler_Cromer_data))

    def Verlet_integrator(self):
        position_series = [self.init_x]
        velocity_series = [self.init_v]

        D = 2 * self.m + self.b * self.h
        B = ((self.b * self.h) - (2 * self.m)) / D
        A = 2 * (2 * self.m - (self.k * self.h ** 2)) / D

        a_0 = (-self.b / self.m) * self.init_v + (-self.k / self.m) * self.init_x
        x_1 = self.init_x + self.init_v * self.h + 0.5 * a_0 * self.h ** 2  # obtained using a Taylor expansion of order 2
        position_series.append(x_1)

        for counter in range(1, self.no_steps, 1):
            position_series.append(A * position_series[counter] + B * position_series[counter - 1])

        # calculating velocities using an approximation of O(h^2)
        # the velocity is estimated using the mean value theorem
        for counter in range(1, self.no_steps, 1):
            velocity_series.append(
                (position_series[counter + 1] - position_series[counter - 1]) / (2 * self.h))  # +O(h^2)
        position_series = position_series[:len(position_series) - 1]

        self.Verlet_data = [position_series, velocity_series, self.energy_function(position_series, velocity_series)]

    def energy_function(self, position, velocity):
        temp_pos = np.array(position)
        temp_vel = np.array(velocity)
        return 0.5 * self.m * temp_vel ** 2 + 0.5 * self.k * temp_pos ** 2

    def convert_array(self, array):
        # operates on 1 dimensional arrays
        temp = []
        for entry in array:
            temp.append(entry)

        return temp

    def analytic_solution(self):
        # creates the analytic solution position series
        t_0 = 0
        for counter in range(0, self.no_steps, 1):
            self.analytic_series_pos.append(self.ana_position(t_0))
            self.analytic_series_vel.append(self.ana_velocity(t_0))
            t_0 += self.h
        print("solution found")
        self.analytic_energy = self.energy_function(self.analytic_series_pos, self.analytic_series_vel)

    def solver(self):
        #
        print((self.k / self.m))
        print((self.b ** 2 / (4 * self.m ** 2)))
        #
        A = 0
        B = 0
        temp = (self.b ** 2 / (4 * self.m ** 2))
        if self.b == 0:
            marker = 1
            print("The analytic solution is a simple harmonic motion")
            omega = np.sqrt(self.k / self.m)
            A = self.init_v / omega
            B = self.init_x

            # x = A*sin(omega*t) + B*cos(omega*t)
            # v = omega A cos(omega t) - omega*B*sin(omega t)

            self.__coefficients = [omega, 0, marker, A, B]
        elif (self.k / self.m) > temp:  # imaginary
            print("The solution is a lightly damped oscillation")
            marker = 3

            p = -1 * self.b / 2 * self.m
            q = np.sqrt((self.k / self.m) - self.b ** 2 / (4 * self.m ** 2))
            # initial conditions
            A = (- self.init_x * p + self.init_v) / q
            B = self.init_x
            self.__coefficients = [p, q, marker, A, B]
        elif (self.k / self.m) == temp:
            print("The solution is a critically damped oscillation")
            marker = 2  # repeated real solutions
            K = -1 * self.b / 2 * self.m
            # initial conditions
            A = self.init_x
            self.__coefficients = [K, 0, marker, A, B]
        elif (self.k / self.m) < temp:  # overdamped oscillation
            marker = 4
            print("The solution is an overdamped oscillation")
            p = -1 * self.b / 2 * self.m + np.sqrt(-(self.k / self.m) + self.b ** 2 / (4 * self.m ** 2))
            q = -1 * self.b / 2 * self.m - np.sqrt(-(self.k / self.m) + self.b ** 2 / (4 * self.m ** 2))

            # initial conditions
            A = (q * self.init_x - self.init_v) / (q - p)
            B = self.init_x - A
            self.__coefficients = [p, q, marker, A, B]
        print(marker)

    def ana_position(self, t):
        k_1 = self.__coefficients[0]
        k_2 = self.__coefficients[1]
        marker = self.__coefficients[2]
        A = self.__coefficients[3]
        B = self.__coefficients[4]
        if marker == 1:  # no damping solution
            return A * np.sin(self.natural_angular_frequency * t) + B * np.cos(self.natural_angular_frequency * t)
        elif marker == 2:  # regular damping (complex)
            return A * np.exp(k_1 * t)
        elif marker == 3:  # repeated root
            return np.exp(k_1 * t) * (A * np.sin(k_2 * t) + B * np.cos(k_2 * t))
        elif marker == 4:  # two real distinct solutions
            return A * np.exp(k_1 * t) + B * np.exp(k_2 * t)

    def ana_velocity(self, t):
        k_1 = self.__coefficients[0]
        k_2 = self.__coefficients[1]
        marker = self.__coefficients[2]
        A = self.__coefficients[3]
        B = self.__coefficients[4]
        if marker == 1:  # no damping solution
            return A * k_1 * np.cos(k_1 * t) - B * k_1 * np.sin(k_1 * t)
        elif marker == 2:  # regular damping (complex)
            return A * k_1 * np.exp(k_1 * t)
        elif marker == 3:  # repeated root
            return k_1 * np.exp(k_1 * t) * (A * np.sin(k_2 * t) + B * np.cos(k_2 * t)) + np.exp(k_1 * t) * (
                    A * k_2 * np.cos(k_2 * t) - B * k_2 * np.sin(k_2 * t))
        elif marker == 4:  # two real distinct solutions
            return A * k_1 * np.exp(k_1 * t) + B * k_2 * np.exp(k_2 * t)

    def plot_data(self):  # plots all on separate graphs
        # analytical solution
        figure3 = plt.figure()
        axes_5 = figure3.add_subplot(121)
        axes_5.plot(self.analytic_series_pos, self.analytic_series_vel, label="Analytical")
        axes_5.set_xlabel("position/ m")  # edit later if the functions don't exist
        axes_5.set_ylabel("velocity/ ms^-1")  # as above
        # energy plotting
        axes_6 = figure3.add_subplot(122)
        axes_6.plot(self.time, self.analytic_energy)
        axes_6.set_xlabel("time/ s")
        axes_6.set_ylabel("energy/ J")
        figure3.legend()

        # Euler method
        # plotting
        figure = plt.figure()
        axes_1 = figure.add_subplot(121)
        axes_1.plot(self.Euler_data[0], self.Euler_data[1], label="Euler")
        axes_1.set_xlabel("position/ m")  # edit later if the functions don't exist
        axes_1.set_ylabel("velocity/ ms^-1")  # as above
        # energy plotting

        axes_2 = figure.add_subplot(122)
        axes_2.plot(self.time, self.Euler_data[2])
        axes_2.set_xlabel("time/ s")
        axes_2.set_ylabel("energy/ J")
        figure.legend()
        # end plotting

        # Better Euler method
        # plotting
        figure2 = plt.figure()
        axes_3 = figure2.add_subplot(121)
        axes_3.plot(self.B_Euler_data[0], self.B_Euler_data[1], label="Better Euler")
        axes_3.set_xlabel("position/ m")  # edit later if the functions don't exist
        axes_3.set_ylabel("velocity/ ms^-1")  # as above
        # energy plotting
        axes_4 = figure2.add_subplot(122)
        axes_4.plot(self.time, self.B_Euler_data[2])
        axes_4.set_xlabel("time/ s")
        axes_4.set_ylabel("energy/ J")
        figure2.legend()
        # end plotting

        # Verlet method
        figure4 = plt.figure()
        axes_7 = figure4.add_subplot(121)
        axes_7.plot(self.Verlet_data[0], self.Verlet_data[1], label="Verlet")
        axes_7.set_xlabel("position/ m")
        axes_7.set_ylabel("velocity/ ms^-1")
        # energy plotting
        axes_8 = figure4.add_subplot(122)
        axes_8.plot(self.time, self.Verlet_data[2])
        axes_8.set_xlabel("time/ s")
        axes_8.set_ylabel("energy/ J")
        figure4.legend()

        # Euler Cromer method
        figure5 = plt.figure()
        axes_9 = figure5.add_subplot(121)
        axes_9.plot(self.Euler_Cromer_data[0], self.Euler_Cromer_data[1], label="Euler Cromer Method")
        axes_9.set_xlabel("position/ m")
        axes_9.set_ylabel("velocity/ ms^-1")
        # energy plotting
        axes_10 = figure5.add_subplot(122)
        axes_10.plot(self.time, self.Euler_Cromer_data[2])
        axes_10.set_xlabel("time/ s")
        axes_10.set_ylabel("energy/ J")
        figure5.legend()

    def plot_single(self):
        figure = plt.figure()
        axes_1 = figure.add_subplot(121)
        axes_1.set_xlabel("position/ m")
        axes_1.set_ylabel("velocity/ ms^-1")
        # energy_function
        energy = []
        axes_2 = figure.add_subplot(122)
        axes_2.set_xlabel("time/ s")
        axes_2.set_ylabel("Energy/ J")

        axes_1.plot(self.Euler_data[0], self.Euler_data[1], label="Euler")
        axes_2.plot(self.time, self.Euler_data[2])

        axes_1.plot(self.B_Euler_data[0], self.B_Euler_data[1], label="Improved Euler")
        axes_2.plot(self.time, self.B_Euler_data[2])

        axes_1.plot(self.Verlet_data[0], self.Verlet_data[1], label="Verlet")
        axes_2.plot(self.time, self.Verlet_data[2])

        axes_1.plot(self.Euler_Cromer_data[0], self.Euler_Cromer_data[1], label="Euler Cromer")
        axes_2.plot(self.time, self.Euler_Cromer_data[2])

        axes_1.plot(self.analytic_series_pos, self.analytic_series_vel, label="Analytic solution")
        axes_2.plot(self.time, self.analytic_energy)
        figure.legend()

    def save_data(self):
        # all files are saved as json dictionaries in the format "name of method": [position, velocity, energy]
        # the header of the file headers named appropriately contains the
        temp = [self.analytic_series_pos, self.analytic_series_vel, self.convert_array(self.analytic_energy)]

        data = {}
        data["Analytic"] = temp
        data["Euler"] = [self.Euler_data[0], self.Euler_data[1], self.convert_array(self.Euler_data[2])]
        data["Better Euler"] = [self.B_Euler_data[0], self.B_Euler_data[1], self.convert_array(self.B_Euler_data[2])]
        data["Verlet"] = [self.Verlet_data[0], self.Verlet_data[1], self.convert_array(self.Verlet_data[2])]
        data["Euler Cromer"] = [self.Euler_Cromer_data[0], self.Euler_Cromer_data[1],
                                self.convert_array(self.Euler_Cromer_data[2])]
        data["coefficients"] = [self.h, self.no_steps, self.b, self.m, self.k, self.init_x,
                                self.init_v]  # h, T, b, m, k, x, v

        with open(self.fileNameSave, 'w') as outfile:
            json.dump(data, outfile)
        outfile.close()

    def load_data(self):
        try:
            with open(self.fileNameLoad) as json_file:
                data = json.load(json_file)
                self.Euler_data = data["Euler"]
                self.B_Euler_data = data["Better Euler"]
                self.Verlet_data = data["Verlet"]
                self.analytic_series_pos = data["Analytic"][0]
                self.analytic_series_vel = data["Analytic"][1]
                self.analytic_energy = data["Analytic"][2]
                self.Euler_Cromer_data = data["Euler Cromer"]
                self.h, self.no_steps, self.b, self.m, self.k, self.init_x, self.init_v = data["coefficients"]
                print(self.h)
                print(self.no_steps)
                print(self.b)
                print(self.m)
                print(self.k)
                print(self.init_x)
                print(self.init_v)

            json_file.close()
            return True
        except:
            print("The file was not found")
            return False

    def find_accuracy(self):
        # finds the accuracy of the simulation by using the analytic energy as a baseline
        # this assigns a number of "fictitious energy" and also graphs the growth of the errors with time
        figure = plt.figure()
        axes_1 = figure.add_subplot(111)
        axes_1.set_ylabel("Energy error/ J")
        axes_1.set_xlabel("time/ s")
        fict_energy = []
        baseline = np.array(self.analytic_energy)

        # Euler's method
        temp = np.abs(np.array(self.Euler_data[2]) - baseline)
        axes_1.plot(self.time, temp, label="Euler")
        fict_energy.append(np.sum(temp))
        # Better Euler
        temp = np.abs(np.array(self.B_Euler_data[2]) - baseline)
        axes_1.plot(self.time, temp, label="Improved Euler")
        fict_energy.append(np.sum(temp))
        # Cromer
        temp = np.abs(np.array(self.Euler_Cromer_data[2]) - baseline)
        axes_1.plot(self.time, temp, label="Euler Cromer")
        fict_energy.append(np.sum(temp))
        # Verlet
        temp = np.abs(np.array(self.Verlet_data[2]) - baseline)
        axes_1.plot(self.time, temp, label="Verlet")
        fict_energy.append(np.sum(temp))
        temp = 0
        print("Euler: " + str(fict_energy[0]) + " J")
        print("Improved Euler: " + str(fict_energy[1]) + "J")
        print("Euler Cromer: " + str(fict_energy[2]) + "J")
        print("Verlet: " + str(fict_energy[3]) + "J")
        figure.legend()
        axes_1.set_title("h = " + str(self.h))

    def const_dist_Verlet_integrator(self, force, min, max):
        position_series = [self.init_x]
        velocity_series = [self.init_v]

        D = 2 * self.m + self.b * self.h
        B = (self.b * self.h - 2 * self.m) / D
        A = 2 * (2 * self.m - self.k * self.h ** 2) / D

        a_0 = (-self.b / self.m) * self.init_v + (-self.k / self.m) * self.init_x
        x_1 = self.init_x + self.init_v * self.h + 0.5 * a_0 * self.h ** 2  # obtained using a Taylor expansion of order 2
        position_series.append(x_1)

        for counter in range(1, self.no_steps, 1):
            if (counter * self.h > min) and (counter * self.h < max):
                position_series.append(
                    A * position_series[counter] + B * position_series[counter - 1] + (force / self.m * self.h ** 2))
            else:
                position_series.append(A * position_series[counter] + B * position_series[counter - 1])

        # calculating velocities using an approximation of O(h^2)
        # the velocity is estimated using the mean value theorem
        # the velocity is independent of the equation of motion. It just utilises the definition of velocity. If h is small
        # enough this approximation holds true
        for counter in range(1, self.no_steps, 1):
            velocity_series.append(
                (position_series[counter + 1] - position_series[counter - 1]) / (2 * self.h))  # +O(h^2)
        position_series = position_series[:len(position_series) - 1]

        self.disturbed_Verlet_data = [position_series, velocity_series,
                                      self.energy_function(position_series, velocity_series)]

    def funct_dist_Verlet_integrator(self, min, max, Amp, freq):
        position_series = [self.init_x]
        velocity_series = [self.init_v]

        D = 2 * self.m + self.b * self.h
        B = (self.b * self.h - 2 * self.m) / D
        A = 2 * (2 * self.m - self.k * self.h ** 2) / D

        a_0 = (-self.b / self.m) * self.init_v + (-self.k / self.m) * self.init_x
        x_1 = self.init_x + self.init_v * self.h + 0.5 * a_0 * self.h ** 2  # obtained using a Taylor expansion of order 2
        position_series.append(x_1)

        for counter in range(1, self.no_steps, 1):
            if (counter * self.h > min) and (counter * self.h < max):
                position_series.append(
                    A * position_series[counter] + B * position_series[counter - 1] + (
                            Amp * np.sin(freq * counter * self.h) / self.m * self.h ** 2))
            else:
                position_series.append(A * position_series[counter] + B * position_series[counter - 1])

        # calculating velocities using an approximation of O(h^2)
        # the velocity is estimated using the mean value theorem
        # the velocity is independent of the equation of motion. It just utilises the definition of velocity. If h is small
        # enough this approximation holds true
        for counter in range(1, self.no_steps, 1):
            velocity_series.append(
                (position_series[counter + 1] - position_series[counter - 1]) / (2 * self.h))  # +O(h^2)
        position_series = position_series[:len(position_series) - 1]

        self.disturbed_Verlet_data = [position_series, velocity_series,
                                      self.energy_function(position_series, velocity_series)]

    def push_testing(self, min, max, force, amp, freq):
        # constant force
        self.const_dist_Verlet_integrator(force, min, max)
        constant_data = self.disturbed_Verlet_data

        # sinusoidal force
        self.funct_dist_Verlet_integrator(min, max, amp, freq)
        function_data = self.disturbed_Verlet_data

        figure = plt.figure()
        axes_1 = figure.add_subplot(111)
        axes_1.plot(self.time, constant_data[0], label="constant force")
        # axes_1.plot(self.time, function_data[0], label="sinusoidal force")
        axes_1.plot(self.time, self.Verlet_data[0], label="undisturbed")
        # axes_1.plot(self.time, self.analytic_series_pos, label="analytical solution")
        figure.legend()

    def search(self, arr, x):
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

    def resonance_Plot(self):
        def test_func(x, a, b):
            return a * np.sin(b * x)

        temp = np.sqrt(self.k / self.m)
        freq_array = []
        for counter in range(0, 40, 1):
            freq_array.append(counter * temp * 0.1)

        amplitude = []
        for freq in freq_array:
            # get the data sets for a specific frequency
            self.funct_dist_Verlet_integrator(0, self.no_steps * self.h, self.init_x, freq)
            temp = self.disturbed_Verlet_data[0]
            # calculate amplitude
            params, params_covariance = optimize.curve_fit(test_func, self.time, temp, p0=[2, 2])
            amplitude.append((params[0]))
            # append to the arrays
        # plot the results
        figure = plt.figure()
        axes_1 = figure.add_subplot(111)
        axes_1.set_xlabel("Frequency/ Hz")
        axes_1.set_ylabel("Amplitude/ m")
        axes_1.plot(freq_array, amplitude, "+b")

    def Critical(self):
        temp = self.b
        b = [0.5 * self.b_critical, self.b_critical, 2 * self.b_critical]
        data = []
        for entry in b:
            self.b = entry
            self.Verlet_integrator()
            data.append(self.Verlet_data)

        # Verlet method
        # phase plots
        figure4 = plt.figure()
        axes_7 = figure4.add_subplot(121)
        axes_7.plot(data[0][0], data[0][1], label="0.5b")
        axes_7.plot(data[1][0], data[1][1], label="b")
        axes_7.plot(data[2][0], data[2][1], label="2b")
        axes_7.set_xlabel("position (m)")
        axes_7.set_ylabel("velocity (ms^-1)")
        # energy plotting
        axes_8 = figure4.add_subplot(122)
        axes_8.plot(self.time, data[0][2])
        axes_8.plot(self.time, data[1][2])
        axes_8.plot(self.time, data[2][2])
        axes_8.set_xlabel("time (s)")
        axes_8.set_ylabel("energy (J)")
        # position plot
        figure = plt.figure()
        axes_1 = figure.add_subplot(111)
        axes_1.plot(self.time, data[0][0], label="0.5b")
        axes_1.plot(self.time, data[1][0], label="b")
        axes_1.plot(self.time, data[2][0], label="2b")
        axes_1.set_ylabel("Position (m)")
        axes_1.set_xlabel("Time (s)")
        figure4.legend()
        figure.legend()

        self.b = temp

    def thing(self):
        plt.plot(self.time, self.Verlet_data[0], label="Verlet")
        plt.plot(self.time, self.analytic_series_pos, label="analytic")
        plt.plot(self.time, self.Euler_Cromer_data[0], label="Eruler Cromer")
        plt.plot(self.time, self.Euler_data[0], label="Euler")
        plt.plot(self.time, self.B_Euler_data[0], label="Better euler")

        plt.legend()


def getData():  # edit this function

    fileName = float(input("Enter the name of the simulation (don't forget the format)"))
    m = float(input("Enter the value for the mass of the particle: "))
    k = float(input("Enter the value of the spring constant: "))
    b = float(input("Enter the value of the damping constant "))
    T = float(input("Enter the time you want the simulation to run: "))
    h = float(input("Enter the time step in seconds: "))
    init_x = float(input("Enter the initial position"))
    init_v = float(input("Enter the initial velocity"))

    return m, k, b, T, h, init_x, init_v


def main():
    option = 0
    while option != "8":
        option = input("Select an option:")
        print("1. Run simulation")
        print("2. Load old simulation")
        if option == "1" or option == "2":
            if option == "1":
                m, k, b, T, h, init_x, init_v = getData()
                os = SHO(h, T, b, m, k, init_x, init_v)
            elif option == "2":
                print("Enter the name of the file or type 'none' if you want to use default")
                name = input()
                if name == "none":
                    os = SHO(0, 0, 0, 0, 0, 0, 0)
                    os.load_data()
                else:
                    os = SHO(0, 0, 0, 0, 0, 0, 0, name)
                    check = os.load_data()
                    if not check:
                        print("Goodbye!")
                        return 0
            option = input("Select an option:")
            print("3. Run critical damping simulation")
            print("4. Run simulation with the force appplied")
            print("5. Plot all of it")
            print("6. Save the simulation")
            print("7. Plot Resonance Curves")
            print("8. Leave")
            if option == "3":
                os.Critical()
            elif option == "4":
                min = input("Minimum time: ")
                max = input("Maximum time: ")
                force = input("Force magnitude: ")
                Amp = input("Sinusoidal force amplitude: ")
                freq = input("Sinusoidal force frequency: ")
                os.push_testing(min, max, force, Amp, freq)
            elif option == "5":
                os.plot_data()
                os.plot_single()
            elif option == "6":
                os.save_data()
                print("File was saved as data.txt")
            elif option == "7":
                os.resonance_Plot()
        else:
            print("Goodbye!")


def main2():
    m = 5.44
    k = 0.93
    b = 0.1
    T = 400  # max time
    h = 0.001  # time step

    os = SHO(h, T, b, m, k, 10, 0)

    os.runSimulation()
    os.thing()

    # os.push_testing(45.6, 100, 5, 5, 0.062832)
    # os.push_testing(57, 100, 5, 5, 6.2832)


main2()
plt.show()
