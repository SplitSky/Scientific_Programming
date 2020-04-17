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


class SHM(object):
    def __init__(self, time_step, max_time, b=0, m=0, k=0, init_x=0, init_v=0):
        self.b = b
        self.m = m
        self.k = k
        self.position_series = [init_x]
        self.velocity_series = [init_v]
        self.h = time_step
        self.init_v = init_v
        self.init_x = init_x
        self.no_steps = int(np.rint(max_time / time_step))
        self.

    def Euler_integrator(self):
        for counter in range(1, self.no_steps):
            v_n = self.velocity_series[len(self.velocity_series) - 1]
            x_n = self.position_series[len(self.position_series) - 1]
            a_n = (-self.b / self.m) * v_n + (-self.k / self.m) * x_n

            self.position_series.append(x_n + self.h * v_n)
            self.velocity_series.append(v_n + self.h * a_n)


    def Better_Euler_integrator(self):
        for counter in range(1, self.no_steps):
            v_n = self.velocity_series[len(self.velocity_series) - 1]
            x_n = self.position_series[len(self.position_series) - 1]
            a_0 = (-self.b / self.m) * v_n + (-self.k / self.m) * x_n

            self.position_series.append(x_n + self.h * v_n + 0.5*self.h**2 * a_n)
            self.velocity_series.append(v_n + self.h * a_n)

    def Verlet_integrator(self):
        D = 2*self.m + self.b*self.h
        B = (self.b*self.h - 2*self.m) / D
        A = 2*(2*m - self.k*self.h**2)/D

        a_0 = (-self.b / self.m) * self.init_v + (-self.k / self.m) * self.init_x
        x_1 =  self.init_x + self.init_v*self.h + 0.5*a_0*self.h**2 # obtained using a Taylor expansion of order 2


        for counter in range(2, self.no_steps):


    def Energy(self):


    def plot(self,x, y):
        figure = plt.figure()
        axes_1 = figure.add_subplot(111)
        axes_1.plot(x, y)
        plt.xlabel("x / m")
        plt.ylabel("v / ms^-1")


