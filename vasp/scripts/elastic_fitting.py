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
import math
import sys, os
curPath_soft = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath_soft)
import fitting as fes

curPath = os.getcwd()

def get_poscar(filename):
    poscar = open(filename, 'r')
    poscar_lines = poscar.readlines()
    pos = []
    pos_coordinates = []
    for flag_lines in range(2, 5):
        for i in range(0, 3):
            pos.append(float(poscar_lines[flag_lines].split()[i]))
    for flag_line in range(5, len(poscar_lines)):
        pos_coordinates.append(poscar_lines[flag_line])
    return pos, pos_coordinates

# get the energy from OUTCAR file
def get_energy():
    outcar = open("OUTCAR", 'r')
    outcar_lines = outcar.readlines()
    energy = []
    for flag_lines in outcar_lines:
        if "energy  without entropy" in flag_lines:
            energy.append(flag_lines.split("=")[-1])
    return float(energy[-1])

def get_volume():
    outcar = open("OUTCAR", 'r')
    outcar_lines = outcar.readlines()
    volume = []
    for flag_lines in outcar_lines:
        if "volume of cell" in flag_lines:
            volume.append(flag_lines.split(":")[-1])
    return float(volume[-1])
# fitting to get the energy-strain curve.
def elastic_fitting():
    delta = 0.005
    file_deform = ["c11", "c22", "c12", "c66"]
    energy_c11 = []
    energy_c22 = []
    energy_c12 = []
    energy_c66 = []
    volume = []
    for j in range(0, 4):
        for i in range(-6, 7):
            if j == 0:
                os.chdir(curPath + "/" + str(file_deform[j]) + "/" + str('strain') + str(delta * i))
                energy_c11.append(get_energy())
                volume.append(get_volume())
            elif j == 1:
                os.chdir(curPath + "/" + str(file_deform[j]) + "/" + str('strain') + str(delta * i))
                energy_c22.append(get_energy())
            elif j == 2:
                os.chdir(curPath + "/" + str(file_deform[j]) + "/" + str('strain') + str(delta * i))
                energy_c12.append(get_energy())
            elif j == 3:
                os.chdir(curPath + "/" + str(file_deform[j]) + "/" + str('strain') + str(delta * i))
                energy_c66.append(get_energy())
            else:
                print("Something wrong happens !")
    return energy_c11, energy_c22, energy_c12, energy_c66, volume

# output the energy-strain data
def output_energy():
    ev_to_j = 1.60217
    energy_c11, energy_c22, energy_c12, energy_c66, volume = elastic_fitting()
    volume1 = volume[6]
    v0 = ev_to_j * 100 / volume1
    e0 = energy_c11[6]
    strain = [-0.03, -0.025, -0.02, -0.015, -0.01, -0.005, 0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03]
    with open(curPath+"/Energy-strain.dat", "w", encoding='utf-8') as ev:
        ev.write("strain".ljust(12, ' ')+"C11".ljust(12, ' ')+"C22".ljust(12, ' ') + "C12".ljust(12, ' ') + "C66".ljust(12, ' ') + '\n')
        for i in range(0, 13):
            ev.write(str(strain[i])[0:10].ljust(12, ' ') + str((energy_c11[i]-e0))[0:8].ljust(12, ' ') +
                      str((energy_c22[i]-e0))[0:8].ljust(10, ' ') + str((energy_c12[i]-e0))[0:8].ljust(10, ' ')
                      + str((energy_c66[i]-e0))[0:8].ljust(10, ' ')+'\n')
    ev.close()

# plot the energy-strain curve
def energy_strain_plot():
    data = np.loadtxt(curPath + "/Energy-strain.dat", skiprows=1, dtype=np.float)
    energy_c11 = data[:, 1]
    energy_c22 = data[:, 2]
    energy_c12 = data[:, 3]
    energy_c66 = data[:, 4]
    strain = data[:, 0]
    plt.plot(strain, energy_c11, color="red", label = "uniaxial x direction")
    plt.plot(strain, energy_c22, color="blue", label = "uniaxial y direction")
    plt.plot(strain, energy_c12, color="black", label="biaxial")
    plt.plot(strain, energy_c66, color="cyan", label="shear")
    plt.legend(loc="best")
    plt.ylabel("Strain Energy (eV)", fontsize=14, fontname='arial')
    plt.xlabel("Strain", fontsize=14, fontname='arial')
    plt.xticks(strain)
    plt.savefig(curPath + "/ES.png", dpi=300)
    plt.show()

def energy_strain_fitting_plot():
    ev_to_j = 1.60217
    energy_c11, energy_c22, energy_c12, energy_c66, volume = elastic_fitting()
    volume1 = volume[6]
    v0 = ev_to_j * 100 / volume1
    pos, pos_coordinates = get_poscar("POSCAR")
    thickness = pos[-1]/10.0
    data = np.loadtxt(curPath + "/Energy-strain.dat", skiprows=1, dtype=np.float)
    strain = data[:, 0]
    label = ["uniaxial x direction", "uniaxial y direction", "biaxial", "shear"]
    label_text = ["C11", "C22", "C12", "C66"]
    color_scatter = ["black", "black", "black", "black"]
    color_fit = ["blue", "cyan", "orange", "red"]
    plt.figure(figsize=(8, 6))
    y_text = [0.8, 0.7, 0.6, 0.5]
    a_nihe = []
    d, e, f = fes.parabola_fitting(strain / 2, data[:, 4] * v0, label[3])
    for i in range(1, 5):
        a, b, c = fes.parabola_fitting(strain, data[:, i] * v0, label[i-1])
        g, h, k = fes.parabola_fitting(strain, data[:, i], label[i-1])
        mid_x = max(strain) - min(strain)
        mid_y = max(data[:, 3]) - min(data[:, 3])
        text_x = min(strain) + mid_x / 5
        text_y = min(data[:, i]) + mid_y * y_text[i-1]
        text_context = " "
        if i == 1:
            a_nihe.append(2*a)
            text_context = label_text[i-1] + " = " + str(round(2*a * thickness, 2))
        elif i == 2:
            a_nihe.append(2 * a)
            if abs(2*a - a_nihe[0])<=1:
                y_text[i ] = 0.7
                y_text[i + 1] = 0.6
                continue
            else:
                text_context = label_text[i - 1] + " = " + str(round(2 * a * thickness, 2))
        elif i == 3:
            text_context = "C12" + " = " + str(round(((a-(sum(a_nihe)/2.0)) * thickness), 2))
        elif i == 4:
            text_context = "C66" + " = " + str(round((d / 2.0) * thickness, 2))
        print("求解的曲线是:")
        print("y = " + str(round(a, 2)) + " x*x + " + str(round(b, 2)) + " x+ " + str(c))
        plt.text(text_x, text_y, text_context, fontname='arial', color="black", fontsize=14)
        plt.scatter(strain, data[:, i], color=color_scatter[i-1], label=label[i-1]+" sample data", linewidth=2)
        x = np.arange(min(strain), max(strain) * 1.01, ((max(strain) - min(strain)) / 40))     # 在0-15直接画100个连续点
        y = g * x * x + h * x + k  # 函数式
        plt.plot(x, y, color=color_fit[i-1], label=label[i-1], linewidth=2)
    # 画拟合直线
    plt.xlabel("Strain", fontsize=14)
    plt.ylabel("Strain Energy (eV/unit cell)", fontsize=14)
    leg = plt.legend(loc = "upper right")  # 绘制图例
    leg.get_frame().set_facecolor('cyan')
    leg.get_frame().set_linewidth(0.0)
    leg.get_frame().set_alpha(0.4)
    #leg.get_frame().set_edgecolor('b')
    plt.savefig(curPath +"/ES-fitting.png", dpi=300)
    plt.show()
if __name__ == '__main__':
    output_energy()
    energy_strain_fitting_plot()