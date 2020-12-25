# -*- coding: utf-8 -*-
"""
Created on 10:14 13-11-2020 

@author: X-Y Ding
mail to: dxy_vasp@163.com
python3: es.py
"""
import numpy as np
from matplotlib import pyplot as plt
import math
import sys, os
curPath_soft = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath_soft)

curPath = os.getcwd()

def get_poscar(filename = "POSCAR_orig"):
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

# deform in xx direction
def deform(delta=0.0, filename = "POSCAR_orig"):
    pos, pos_coordinates = get_poscar(filename)
    poscar = np.array(pos).reshape(3, 3)
    pos_matrix = np.matrix(poscar)
    unit_matrix = np.matrix([[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]])
    deform_matrix = np.matrix([[delta, 0, 0], [0, delta, 0], [0, 0, 0]])
    result = pos_matrix * (unit_matrix + deform_matrix)
    result_array = np.array(result)
    with open("POSCAR", "w", encoding='utf-8') as pos:
        pos.write("POSCAR: for 2D materials" + "\n")
        pos.write(str(1.0) + "\n")
        for i in range(0, 3):
            for j in range(0, 3):
                pos.write("   " + str(result_array[i][j])[:8].ljust(12, '0'))
            pos.write('\n')
        for h in range(0, len(pos_coordinates)):
            pos.write(pos_coordinates[h])
    pos.close()

def get_energy():
    outcar = open("OUTCAR", 'r')
    outcar_lines = outcar.readlines()
    energy = []
    for flag_lines in outcar_lines:
        if "energy  without entropy" in flag_lines:
            energy.append(flag_lines.split("=")[-1])
    return float(energy[-1])

def get_volume():
    pos, pos_coordinates = get_poscar("POSCAR")
    outcar = open("OUTCAR", 'r')
    outcar_lines = outcar.readlines()
    volume = []
    for flag_lines in outcar_lines:
        if "volume of cell" in flag_lines:
            volume.append(flag_lines.split(":")[-1])
    return float(volume[-1])/pos[-1]

def output_es():
    delta_x = [float(i * 0.01) for i in range(-10, 11)]
    energy = []
    volume = []
    for i in range(0, len(delta_x)):
        os.chdir(curPath + "/" + "strain" + str(delta_x[i]))
        energy.append(get_energy())
        volume.append(get_volume())
    os.chdir(curPath)
    with open("Energy_area.dat", "w", encoding='utf-8') as es:
        es.write("Energy area data" + "\n")
        for j in range(0, len(energy)):
            es.write("   " + str(volume[j])[:8].ljust(12, '0') + "  " + str(energy[j])[:8].ljust(12, '0') + "\n")
    es.close()
    return energy, volume

# manipulate deform
def manipulate_deform():
    os.system("cp POSCAR POSCAR_orig")
    delta_x = [float(i * 0.01) for i in range(-10, 11)]
    for i in range(0, len(delta_x)):
        deform(delta_x[i], "POSCAR_orig")
        os.system("mkdir " + "strain" + str(delta_x[i]))
        os.system("cp "+ curPath + "/POSCAR " + "./" + "strain" + str(delta_x[i]))
        os.system("cp "+ curPath + "/INCAR " + "./" +"strain" + str(delta_x[i]))
        os.system("cp "+ curPath + "/POTCAR " + "./" + "strain" + str(delta_x[i]))
        os.system("cp "+ curPath +"/KPOINTS " + "./" + "strain" + str(delta_x[i]))
    os.system("cp " + curPath_soft + "/submit_es.sh " + curPath + "/")
    os.system("cp POSCAR_orig POSCAR")

def plot_es():
    data = np.loadtxt("Energy_area.dat", skiprows=1, dtype=float)
    plt.plot(data[:,0], data[:,1], color="blue", ls="-.", label="Energy area")
    plt.legend()
    plt.savefig("Es.png", dpi=300)
    plt.show()
if __name__ == '__main__':
    manipulate_deform()




