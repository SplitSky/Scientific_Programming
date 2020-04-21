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
    def __init__(self, time_step, max_time, b=0.01, m=1.0, k=1.0, init_x=0.0, init_v=0.0, fileNameSave="data.txt",
                 fileNameLoad="data.txt"):
        '''
        parameter name                  type                    description
        :param time_step:               float                   the time step used in calculations
        :param max_time:                float                   the time the simulation lasts
        :param b:                       float                   damping coefficient
        :param m:                       float                   mass of the oscillator
        :param k:                       float                   spring constant
        :param init_x:                  float                   initial position
        :param init_v:                  float                   initial velocity
        :param fileNameSave:            string                  filename used for saving data
        :param fileNameLoad:            string                  filename used for loading data
        self.no_steps                   float                   number of iterations
        self.natural_angular_frequency  float                   the natural frequency of the oscillation
        self.gamma                      float                   the damping constant
        self.quality_factor             float                   quality factor
        self.analytic_series_pos        array                   the position of the analytical solution
        self.analytic_series_vel        array                   the velocity of the analytical solution
        self.analytic_energy            array                   the array containing energy data of the analytical solution
        self.__coefficients             array                   coefficients used for the analytical solution
        self.b_britical                 float                   the critical damping constant
        self.Euler_data                 array                   the data from Euler's method
        self.B_Euler_data               array                   the data from improved Euler's method
        self.Verley                     array                   the data from Verlet's method
        self.Euler_Cromer_data          array                   the data from Euler-Cromer method
        self.analytical_data            array                   variable for storing analytical data
        self.time                       array                   the time array used in all simulations
        self.disturbed_Verlet_data      array                   the data of the Verlet method with force applied
        '''
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
        # data variables
        self.Euler_data = []
        self.B_Euler_data = []
        self.Verlet_data = []
        self.Euler_Cromer_data = []
        self.analytical_data = []
        self.time = np.array(range(0, self.no_steps, 1)) * self.h
        self.disturbed_Verlet_data = []

    def runSimulation(self):
        '''
        Runs the integrators as a single function
        '''
        self.Euler_integrator()
        self.Better_Euler_integrator()
        self.Verlet_integrator()
        self.Euler_Cromer_integrator()
        print("Simulation has been executed")

    def getCoefficients(self):
        # simple get function
        return self.__coefficients

    def Euler_integrator(self):
        '''
        parameter name                  type                    description
        position_series                 array                   stores the position temporarily
        velocity_series                 array                   stores the velocity temporarily
        v_n                             float                   stores the nth velocity term
        x_n                             float                   stores the nth position term
        a_n                             float                   stores the nth acceleration term
        '''
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
        '''
        parameter name                  type                    description
        position_series                 array                   stores the position temporarily
        velocity_series                 array                   stores the velocity temporarily
        v_n                             float                   stores the nth velocity term
        x_n                             float                   stores the nth position term
        a_n                             float                   stores the nth acceleration term
        '''
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
        '''
        parameter name                  type                    description
        position_series                 array                   stores the position temporarily
        velocity_series                 array                   stores the velocity temporarily
        v_n                             float                   stores the nth velocity term
        x_n                             float                   stores the nth position term
        a_n                             float                   stores the nth acceleration term
        temp                            float                   temporary variable - stores the v_n+1 term of the velocity
        '''
        position_series = [self.init_x]
        velocity_series = [self.init_v]
        for counter in range(1, self.no_steps, 1):
            v_n = velocity_series[len(velocity_series) - 1]
            x_n = position_series[len(position_series) - 1]

            a_0 = (-self.b / self.m) * v_n + (-self.k / self.m) * x_n

            temp = v_n + self.h * a_0  # v_n+1
            velocity_series.append(temp)
            position_series.append(x_n + self.h * temp)

        self.Euler_Cromer_data = [position_series, velocity_series,
                                  self.energy_function(position_series, velocity_series)]

    def Verlet_integrator(self):
        '''
        parameter name                  type                    description
        position_series                 array                   stores the position temporarily
        velocity_series                 array                   stores the velocity temporarily
        v_n                             float                   stores the nth velocity term
        x_n                             float                   stores the nth position term
        a_n                             float                   stores the nth acceleration term
        x_1                             float                   stores the second position of the oscillation
        D                               float                   temporary variable for ease of calculation
        B                               float                   temporary variable for ease of calculation
        A                               float                   temporary variable for ease of calculation
        '''
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
        '''
        parameter name                  type                    description
        temp_pos                        numpy array             stores the position array
        temp_vel                        numpy array             stores the velocity array
        :return: the array containing energy values
        '''
        temp_pos = np.array(position)
        temp_vel = np.array(velocity)
        return 0.5 * self.m * temp_vel ** 2 + 0.5 * self.k * temp_pos ** 2

    def convert_array(self, array):
        # operates on 1 dimensional arrays
        '''
        parameter name                  type                    description
        temp                            numpy array             the array holding the array being converting
        :param array:
        :return: converted array
        '''
        temp = []
        for entry in array:
            temp.append(entry)

        return temp

    def analytic_solution(self):
        # creates the analytic solution position series
        '''
        parameter name                  type                    description
        t_0                             float                   the time that the simulation is at
        '''
        t_0 = 0
        for counter in range(0, self.no_steps, 1):
            self.analytic_series_pos.append(self.ana_position(t_0))
            self.analytic_series_vel.append(self.ana_velocity(t_0))
            t_0 += self.h
        print("solution found")
        self.analytic_energy = self.energy_function(self.analytic_series_pos, self.analytic_series_vel)

    def solver(self):
        '''
        parameter name                  type                    description
        A                               float                   function coefficient
        B                               float                   function coefficient
        marker                          int                     the marker indicating the type of a solution
        p                               float                   function coefficient
        q                               float                   function coefficient
        K                               float                   function coefficient
        '''

        A = 0
        B = 0
        temp = (self.b ** 2 / (4 * self.m ** 2))
        if self.b == 0:
            marker = 1
            print("The analytic solution is a simple harmonic motion")
            omega = np.sqrt(self.k / self.m)
            A = self.init_v / omega
            B = self.init_x
            self.__coefficients = [omega, 0, marker, A, B]
        elif (self.k / self.m) > temp:  # imaginary
            print("The solution is a lightly damped oscillation")
            marker = 3
            p = -1 * self.b / (2 * self.m)
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

    def ana_position(self, t):
        '''
        parameter name                  type                    description
        k_1                             float                   function coefficient
        k_2                             float                   function coefficient
        marker                          int                     marker dictating the solution type
        A                               float                   function coefficient
        B                               float                   function coefficient
        :return returns the position of an analytic solution at time t
        '''
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
        '''
        parameter name                  type                    description
        k_1                             float                   function coefficient
        k_2                             float                   function coefficient
        marker                          int                     marker dictating the solution type
        A                               float                   function coefficient
        B                               float                   function coefficient
        :return returns the velocity of an analytic solution at time t
        '''
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
        '''
        parameter name                  type                    description
        axes_1                          object                  subplot object
        axes_2                          object                  subplot object
        axes_3                          object                  subplot object
        axes_4                          object                  subplot object
        axes_5                          object                  subplot object
        axes_6                          object                  subplot object
        axes_7                          object                  subplot object
        axes_8                          object                  subplot object
        axes_9                          object                  subplot object
        axes_10                         object                  subplot object
        figure                          object                  figure object
        figure2                         object                  figure object
        figure3                         object                  figure object
        figure4                         object                  figure object
        figure5                         object                  figure object
        '''
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
        # plots the single
        '''
        parameter name                  type                    description
        axes_1                          object                  subplot object
        axes_2                          object                  subplot object
        figure                          object                  figure object
        '''
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
        '''
        parameter name                  type                    description
        temp                            array                   stores the analytic data
        data                            dictionary              stores the data to be saved
        '''
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

        with open("data.txt", 'w') as outfile:
            json.dump(data, outfile)
        outfile.close()

    def load_data(self):
        '''
        parameter name                  type                    description
        data                            dictionary              stores the data to be saved
        json_file                       object                  json file object
        '''
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
                self.time = np.array(range(0, self.no_steps, 1)) * self.h

            json_file.close()
            return True
        except:
            print("The file was not found")
            return False

    def find_accuracy(self):
        # finds the accuracy of the simulation by using the analytic energy as a baseline
        # this assigns a number of "fictitious energy" and also graphs the growth of the errors with time
        '''
        parameter name                  type                    description
        axes_1                          object                  subplot object
        fict_energy                     numpy array             stores the error energy
        baseline                        numpy array             stores the analytic energy
        temp                            numpy array             stores the plotting value of the error energy
        figure                          object                  figure object
        '''
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
        print("The energy errors for b = " + str(self.b))
        print("Euler: " + str(fict_energy[0]) + " J")
        print("Improved Euler: " + str(fict_energy[1]) + "J")
        print("Euler Cromer: " + str(fict_energy[2]) + "J")
        print("Verlet: " + str(fict_energy[3]) + "J")
        figure.legend()
        axes_1.set_title("h = " + str(self.h))

    def const_dist_Verlet_integrator(self, force, min, max):
        '''
        parameter name                  type                    description
        position_series                 array                   stores the position temporarily
        velocity_series                 array                   stores the velocity temporarily
        v_n                             float                   stores the nth velocity term
        x_n                             float                   stores the nth position term
        a_n                             float                   stores the nth acceleration term
        x_1                             float                   stores the second position of the oscillation
        D                               float                   temporary variable for ease of calculation
        B                               float                   temporary variable for ease of calculation
        A                               float                   temporary variable for ease of calculation
        '''
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
        '''
        parameter name                  type                    description
        position_series                 array                   stores the position temporarily
        velocity_series                 array                   stores the velocity temporarily
        v_n                             float                   stores the nth velocity term
        x_n                             float                   stores the nth position term
        a_n                             float                   stores the nth acceleration term
        x_1                             float                   stores the second position of the oscillation
        D                               float                   temporary variable for ease of calculation
        B                               float                   temporary variable for ease of calculation
        A                               float                   temporary variable for ease of calculation
        '''
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
        '''
        parameter name                  type                    description
        constant_data                   array                   array of the data modified by a constant force
        function_data                   array                   array of the data modified by a sinusoidal force
        axes_1                          object                  subplot object
        figure                          object                  figure object
        '''
        # constant force
        self.const_dist_Verlet_integrator(force, min, max)
        constant_data = self.disturbed_Verlet_data

        # sinusoidal force
        self.funct_dist_Verlet_integrator(min, max, amp, freq)
        function_data = self.disturbed_Verlet_data

        figure = plt.figure()
        axes_1 = figure.add_subplot(111)
        axes_1.set_xlabel("Time (s)")
        axes_1.set_ylabel("Position (m)")
        axes_1.plot(self.time, constant_data[0], label="constant force")
        axes_1.plot(self.time, function_data[0], label="sinusoidal force")
        axes_1.plot(self.time, self.Verlet_data[0], label="undisturbed")
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
        '''
        parameter name                  type                    description
        temp                            float/array             temporary variable
        freq_array                      array                   array of angular frequencies
        amplitude                       array                   sum of the amplotudes for a resonance plot
        figure                          object                  figure object
        axes_1                          object                  subplot object
        '''

        temp = np.sqrt(self.k / self.m)
        freq_array = []
        for counter in range(0, 200, 1):
            freq_array.append(counter * temp * 0.01)

        amplitude = []
        for freq in freq_array:
            # get the data sets for a specific frequency
            self.funct_dist_Verlet_integrator(0, self.no_steps * self.h, self.init_x * 0.5, freq)
            temp = self.disturbed_Verlet_data[0]
            # calculate amplitude
            temp = np.abs(np.array(temp))
            amplitude.append(np.mean(temp))
            # append to the arrays
        # plot the results
        figure = plt.figure()
        axes_1 = figure.add_subplot(111)
        axes_1.set_xlabel("Frequency/ Hz")
        axes_1.set_ylabel("Amplitude/ m")
        axes_1.plot(freq_array, amplitude)

    def Critical(self):
        '''
        parameter name                  type                    description
        temp                            float                   stores the b value to avoid change
        b                               float                   damping coefficient
        data                            array                   stores the data of all the integrations
        figure4                         object                  figure object
        figure                          object                  figure object
        axes_1                          object                  subplot object
        axes_7                          object                  subplot object
        axes_8                          object                  subplot object
        '''
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

    def complete_Resonance(self):
        '''
        parameter name                  type                    description
        temp 2                           float                   stores the b value to avoid change
        figure                          object                  figure object
        axes_1                          object                  subplot object
        freq_array                      array                   array of angular frequencies
        amplitude                       array                   sum of the amplotudes for a resonance plot
        b_prime                         array                   stores the damping coefficients
        '''
        figure = plt.figure()
        axes_1 = figure.add_subplot(111)
        axes_1.set_xlabel("Frequency/ Hz")
        axes_1.set_ylabel("Amplitude/ m")

        temp2 = self.b
        temp = np.sqrt(self.k / self.m)
        freq_array = []
        for counter in range(0, 220, 1):
            freq_array.append(counter * temp * 0.01)
        b_prime = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 2.2, 4.5, 9]
        amplitude = []
        for b in b_prime:
            self.b = b
            for freq in freq_array:
                # get the data sets for a specific frequency
                self.funct_dist_Verlet_integrator(0, self.no_steps * self.h, self.init_x * 0.5, freq)
                temp = self.disturbed_Verlet_data[0]
                # calculate amplitude
                temp = np.abs(np.array(temp))
                amplitude.append(np.mean(temp))
                # append to the arrays

            axes_1.plot(freq_array, amplitude, label="b = " + str(b) + " kgs^-1")
            amplitude = []
        # plot the results

        figure.legend()
        self.b = temp2

    def find_accuracy_complete(self):
        '''
        parameter name                  type                    description
        temp                            float                   stores the b variable to avoid change
        '''
        temp = self.b
        self.b = 0
        self.find_accuracy()
        self.b = 0.1
        self.find_accuracy()
        self.b = 1
        self.find_accuracy()
        self.b = 2
        self.find_accuracy()
        self.b = 3
        self.find_accuracy()
        self.b = temp


def getData():  # edit this function
    '''
    parameter name                  type                    description
    m                               float                   mass
    k                               float                   spring constant
    b                               float                   damping coefficient
    T                               float                   total time
    h                               float                   time step
    init_x                          float                   initial position
    init_v                          float                   initial velocity
    '''
    m = float(input("Enter the value for the mass of the particle: "))
    k = float(input("Enter the value of the spring constant: "))
    b = float(input("Enter the value of the damping constant "))
    T = float(input("Enter the time you want the simulation to run: "))
    h = float(input("Enter the time step in seconds: "))
    init_x = float(input("Enter the initial position"))
    init_v = float(input("Enter the initial velocity"))
    return m, k, b, T, h, init_x, init_v


def main():
    '''
    parameter name                  type                    description
    option                          int                     option chosen
    name                            string                  file name
    os                              object                  the class object
    check                           boolean                 check variable
    '''
    option = 0
    while option != "9":

        print("1. Run simulation")
        print("2. Load old simulation")
        option = input("Select an option:")
        if option == "1" or option == "2":
            if option == "1":
                m, k, b, T, h, init_x, init_v = getData()
                os = SHO(h, T, b, m, k, init_x, init_v)
            elif option == "2":
                print("Enter the name of the file or type 'none' if you want to use default")
                name = input()
                if name == "none":
                    os = SHO(0.01, 100)
                    os.load_data()
                else:
                    os = SHO(0.01, 100, fileNameLoad=name)
                    check = os.load_data()
                    if not check:
                        print("Goodbye!")
                        return 0

            print("3. Run critical damping simulation")
            print("4. Run simulation with the force appplied")
            print("5. Plot all of it")
            print("6. Save the simulation")
            print("7. Plot Resonance Curves")
            print("8. Run default")
            print("9. Leave")
            option = input("Select an option:")
            if option == "3":
                os.Critical()
            elif option == "4":
                min = float(input("Minimum time: "))
                max = float(input("Maximum time: "))
                force = float(input("Force magnitude: "))
                Amp = float(input("Sinusoidal force amplitude: "))
                freq = float(input("Sinusoidal force frequency: "))
                os.push_testing(min, max, force, Amp, freq)
            elif option == "5":
                os.plot_data()
                os.plot_single()
            elif option == "6":
                os.save_data()
                print("File was saved as data.txt")
            elif option == "7":
                os.complete_Resonance()
            elif option == "8":
                print("Running default simulation")
                os.runSimulation()
                os.push_testing(45.6, 100, 2, 2, 0.062832)  # at zero amplitude
                os.push_testing(57, 100, 2, 2, 6.2832)  # at 3/4 of a cycle
                os.push_testing(53.2, 100, 2, 2, 0.41)
                os.find_accuracy_complete()
                os.Critical()
                os.complete_Resonance()
                os.plot_single()
                os.plot_data()
            plt.show()

        else:
            print("Goodbye!")


main()
plt.show()
