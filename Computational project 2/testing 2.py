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

plt.rcParams.update({'font.size': 14})
plt.style.use('default')
figure = plt.figure()
plt.rcParams.update({'errorbar.capsize': 2})


class SHO(object):
    def __init__(self, time_step, max_time, b=0.0, m=0.0, k=0.0, init_x=0.0, init_v=0.0):
        self.b = b
        self.m = m
        self.k = k
        self.position_series = [init_x]
        self.velocity_series = [init_v]
        self.h = time_step
        self.init_v = init_v
        self.init_x = init_x
        self.no_steps = int(np.rint(max_time / time_step))
        self.energy_series = []
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




    def getCoefficients(self):
        return self.__coefficients

    def Euler_integrator(self):
        for counter in range(1, self.no_steps, 1):
            v_n = self.velocity_series[len(self.velocity_series) - 1]
            x_n = self.position_series[len(self.position_series) - 1]
            a_n = (-self.b / self.m) * v_n + (-self.k / self.m) * x_n

            self.position_series.append(x_n + self.h * v_n)
            self.velocity_series.append(v_n + self.h * a_n)

    def Better_Euler_integrator(self):
        for counter in range(1, self.no_steps, 1):
            v_n = self.velocity_series[len(self.velocity_series) - 1]
            x_n = self.position_series[len(self.position_series) - 1]
            a_0 = (-self.b / self.m) * v_n + (-self.k / self.m) * x_n

            self.position_series.append(x_n + self.h * v_n + 0.5 * self.h ** 2 * a_0)
            self.velocity_series.append(v_n + self.h * a_0)

    def Verlet_integrator(self):
        D = 2 * self.m + self.b * self.h
        B = (self.b * self.h - 2 * self.m) / D
        A = 2 * (2 * self.m - self.k * self.h ** 2) / D

        a_0 = (-self.b / self.m) * self.init_v + (-self.k / self.m) * self.init_x
        x_1 = self.init_x + self.init_v * self.h + 0.5 * a_0 * self.h ** 2  # obtained using a Taylor expansion of order 2

        self.position_series.append(x_1)
        for counter in range(1, self.no_steps, 1):
            self.position_series.append(A * self.position_series[counter] + B * self.position_series[counter - 1])

        # calculating velocities using an approximation of O(h^2)
        # the velocity is estimated using the mean value theorem
        for counter in range(1, self.no_steps, 1):
            self.velocity_series.append(
                (self.position_series[counter + 1] - self.position_series[counter - 1]) / (2 * self.h))  # +O(h^2)
        print(len(self.position_series))
        self.position_series = self.position_series[:len(self.position_series) - 1]
        print(len(self.position_series))

    def Energy(self):
        temp_pos = np.array(self.position_series)
        temp_vel = np.array(self.velocity_series)

        E = 0.5 * self.m * temp_vel ** 2 + 0.5 * self.k * temp_pos ** 2
        self.energy_series = E

    def energy_function(self, position, velocity):
        temp_pos = np.array(position)
        temp_vel = np.array(velocity)

        E = 0.5 * self.m * temp_vel ** 2 + 0.5 * self.k * temp_pos ** 2
        return E

    def clearSeries(self):
        self.position_series = [self.init_x]
        self.energy_series = []
        self.velocity_series = [self.init_v]

    def analytic_solution(self):
        # creates the analytic solution position series
        t_0 = 0
        for counter in range(0, self.no_steps, 1):
            self.analytic_series_pos.append(self.ana_position(t_0))
            self.analytic_series_vel.append(self.ana_velocity(t_0))
            t_0 += self.h
        print("solution found")
        temp_pos = np.array(self.analytic_series_pos)
        temp_vel = np.array(self.analytic_series_vel)

        self.analytic_energy = 0.5 * self.m * temp_vel ** 2 + 0.5 * self.k * temp_pos ** 2

    def solver(self):
        A = 0
        B = 0
        temp = (self.b ** 2 / (4 * self.m ** 2))
        if self.b == 0:
            marker = 1
            print("o")
            omega = np.sqrt(self.k / self.m)
            A = self.init_v / omega
            B = self.init_x

            # x = A*sin(omega*t) + B*cos(omega*t)
            # v = omega A cos(omega t) - omega*B*sin(omega t)

            self.__coefficients = [omega, 0, marker, A, B]

        elif (self.k / self.m) > temp:  # imaginary
            print("a")
            marker = 3
            p = -1 * self.b / 2 * self.m
            q = np.sqrt((self.k / self.m) - self.b ** 2 / (4 * self.m ** 2))
            # initial conditions
            A = self.init_v - self.init_x
            B = self.init_x

            self.__coefficients = [p, q, marker, A, B]

        elif (self.k / self.m) == temp:
            print("b")
            marker = 2  # repeated real solutions
            K = -1 * self.b / 2 * self.m
            # initial conditions
            A = self.init_x
            self.__coefficients = [K, 0, marker, A, B]

        elif (self.k / self.m) < temp:  # normal solution (two real distinct roots)
            marker = 4

            p = -1 * self.b / 2 * self.m + np.sqrt((self.k / self.m) - self.b ** 2 / (4 * self.m ** 2))
            q = -1 * self.b / 2 * self.m - np.sqrt((self.k / self.m) - self.b ** 2 / (4 * self.m ** 2))

            # initial conditions
            A = (q * self.init_x - self.init_v) / (q - p)
            B = self.init_x - A
            self.__coefficients = [p, q, marker, A, B]

    def ana_position(self, t):
        k_1 = self.__coefficients[0]
        k_2 = self.__coefficients[1]
        marker = self.__coefficients[2]
        A = self.__coefficients[3]
        B = self.__coefficients[4]
        if marker == 1:  # no damping solution
            print("d")
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
                    A * k_2 * np.cos(k_2 * t) - B * k_1 * np.sin(k_2 * t))
        elif marker == 4:  # two real distinct solutions
            return A * k_1 * np.exp(k_1 * t) + B * k_2 * np.exp(k_2 * t)

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
    '''

    def plot_data(self):  # plots all on separate graphs

        position = []
        velocity = []
        t = np.array(range(0, self.no_steps, 1)) * self.h
        # analytical solution

        position.append(self.analytic_series_pos)
        velocity.append(self.analytic_series_vel)
        figure3 = plt.figure()
        axes_5 = figure3.add_subplot(121)
        axes_5.plot(position[0], velocity[0], label="Analytical")
        axes_5.set_xlabel("position/ m")  # edit later if the functions don't exist
        axes_5.set_ylabel("velocity/ ms^-1")  # as above
        # energy plotting
        self.Energy()
        axes_6 = figure3.add_subplot(122)
        axes_6.plot(t, self.analytic_energy)
        axes_6.set_xlabel("time/ s")
        axes_6.set_ylabel("energy/ J")
        figure3.legend()

        # Euler method
        self.Euler_integrator()
        position.append(self.position_series)
        velocity.append(self.velocity_series)

        # plotting
        figure = plt.figure()
        axes_1 = figure.add_subplot(121)
        axes_1.plot(position[1], velocity[1], label="Euler")
        axes_1.set_xlabel("position/ m")  # edit later if the functions don't exist
        axes_1.set_ylabel("velocity/ ms^-1")  # as above
        # energy plotting
        self.Energy()
        axes_2 = figure.add_subplot(122)
        axes_2.plot(t, self.energy_series)
        axes_2.set_xlabel("time/ s")
        axes_2.set_ylabel("energy/ J")
        figure.legend()
        # end plotting

        self.clearSeries()  # resets the variables
        # Better Euler method
        self.Better_Euler_integrator()
        position.append(self.position_series)
        velocity.append(self.velocity_series)

        # plotting
        figure2 = plt.figure()
        axes_3 = figure2.add_subplot(121)
        axes_3.plot(position[2], velocity[2], label="Better Euler")
        axes_3.set_xlabel("position/ m")  # edit later if the functions don't exist
        axes_3.set_ylabel("velocity/ ms^-1")  # as above
        # energy plotting
        self.Energy()
        axes_4 = figure2.add_subplot(122)
        axes_4.plot(t, self.energy_series)
        axes_4.set_xlabel("time/ s")
        axes_4.set_ylabel("energy/ J")
        figure2.legend()
        # end plotting

        # Verlet method
        self.clearSeries()
        self.Verlet_integrator()
        position.append(self.position_series)
        velocity.append(self.velocity_series)
        figure4 = plt.figure()
        axes_7 = figure4.add_subplot(121)
        axes_7.plot(position[3], velocity[3], label="Verlet")
        axes_7.set_xlabel("position/ m")  # edit later if the functions don't exist
        axes_7.set_ylabel("velocity/ ms^-1")  # as above
        # energy plotting
        self.Energy()
        axes_8 = figure4.add_subplot(122)
        axes_8.plot(t, self.energy_series)
        axes_8.set_xlabel("time/ s")
        axes_8.set_ylabel("energy/ J")
        figure4.legend()

        self.data = [position, velocity]

    def plot_single(self):
        t = np.array(range(0, self.no_steps, 1)) * self.h
        figure = plt.figure()
        axes_1 = figure.add_subplot(121)
        axes_1.set_xlabel("position/ m")
        axes_1.set_ylabel("velocity/ ms^-1")
        # energy_function
        energy = []
        axes_2 = figure.add_subplot(122)
        axes_2.set_xlabel("time/ s")
        axes_2.set_ylabel("Energy/ J")
        for counter in range(0, len(self.data[0])):
            E = self.energy_function(self.data[0][counter], self.data[1][counter])
            axes_1.plot(self.data[0][counter], self.data[1][counter])
            energy.append(E)
            axes_2.plot(t, E)



def main():
    m = 1.2
    k = 0.8
    b = 0.01
    T = 50  # max time
    h = 0.0001  # time step

    a = SHO(h, T, b, m, k, 10, 0)
    a.plot_data()
    a.plot_single()


main()
plt.show()
