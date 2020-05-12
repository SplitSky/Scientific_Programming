# -*- coding: utf-8 -*-
"""
Author: Tomasz Neska
Date: 06/03/2020
Description: It's capitalism cyka
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


class person():
    def __init__(self, money, trade_ammount, ID):
        self.__money = money
        self.trade_ammount = trade_ammount
        self.ID = ID

    def give_money(self):
        self.__money -= self.trade_ammount
        return self.trade_ammount

    def take_money(self, amount):
        self.__money += amount

    def getMoney(self):
        return self.__money

    def getTradeMoney(self):
        return self.trade_ammount


def bubbleSort(arr):
    n = len(arr)

    # Traverse through all array elements
    for i in range(n - 1):
        # range(n) also work but outer loop will repeat one time more than needed.

        # Last i elements are already in place
        for j in range(0, n - i - 1):

            # traverse the array from 0 to n-i-1
            # Swap if the element found is greater
            # than the next element
            if arr[j].getMoney() > arr[j + 1].getMoney():
                arr[j], arr[j + 1] = arr[j + 1], arr[j]
    return arr
    # Driver code to test above


def capitalism(N, capita, trade):
    a = []
    for i in range(0, N, 1):
        a.append(person(capita, trade, i))

    max_time = 100000
    for i in range(0, max_time, 1):
        print("Iteration!")
        for human in a:

            index = random.randint(0, len(a) - 1)
            if index == human.ID:
                if (index + 1) == len(a):
                    index = 0
                else:
                    index += 1

            # actual trade
            if not human.getMoney() == 0:
                human.give_money()
                a[index].take_money(human.getTradeMoney())
            # end of trade

    a = bubbleSort(a)
    temp = []
    for entry in a:
        print("Person " + str(entry.ID) + " has " + str(entry.getMoney()))
        temp.append(entry.getMoney())

    x = range(0, len(temp), 1)
    plt.plot(x, temp)


capitalism(1000, 100, 1)
plt.show()
