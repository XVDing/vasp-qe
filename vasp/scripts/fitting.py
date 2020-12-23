# -*- coding: utf-8 -*-
"""
Created on 10:14 13-11-2020

@author: XYong Ding
mail to: dxy_vasp@163.com
python3: elastic_fitting.py
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import leastsq
import sys, os
curPath_soft = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath_soft)


# ******************************** linear fitting ***********************************


# ******************************** parabola fitting ***********************************
def parabola_fitting(data_x, data_y, label):
    # 二次函数的标准形式
    def func(params, x):
        a, b, c = params
        return a * x * x + b * x + c

    # 误差函数，即拟合曲线所求的值与实际值的差
    def error(params, x, y):
        return func(params, x) - y

    # 对参数求解
    def slovePara():
        p0 = [10, 10, 10]
        Para = leastsq(error, p0, args=(data_x, data_y))
        return Para

    # 输出最后的结果
    def solution():
        Para = slovePara()
        a, b, c = Para[0]
        return a, b, c
    a, b, c = solution()
    return a, b, c

# plot the fitting line
def para_plot(data_x, data_y, label):
    a, b, c = parabola_fitting(data_x, data_y, label)
    plt.figure(figsize=(8, 6))
    plt.scatter(data_x, data_y, color="green", label="sample data", linewidth=2)
    mid_x = max(data_x) - min(data_x)
    mid_y = max(data_y) - min(data_y)
    text_x = min(data_x) + mid_x / 5
    text_y = min(data_y) + mid_y / 1.2
    text_context = "y = " + str(round(a, 2)) + " x$^2$+ " + str(round(b, 2)) + " x+ " + str(round(c, 3))
    plt.text(text_x, text_y, text_context, fontname='arial', color="black", fontsize=14)
    # 画拟合直线
    x = np.arange(min(data_x), max(data_x) * 1.02, ((max(data_x) - min(data_x)) / 40))  ##在0-15直接画100个连续点
    y = a * x * x + b * x + c  ##函数式
    plt.plot(x, y, color="red", label="parabola fitting", linewidth=2)
    plt.legend()  # 绘制图例
    plt.savefig(label + ".png", dpi=300)
    plt.show()

if __name__ == '__main__':
    data = np.loadtxt("Energy-straineV.dat", skiprows=1, dtype=np.float)
    para_plot(data[:, 0], data[:, 1], "xx")

