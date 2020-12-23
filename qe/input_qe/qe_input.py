# -*- coding: utf-8 -*-
"""
Created on 19:54 04-11-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: input.py
"""
import numpy as np
from matplotlib import pyplot as plt
import sys, os
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)
import math

file_path = curPath

#------------------ FONT_setup ----------------------
font = {'family': 'arial',
        'color': 'black',
        'weight': 'normal',
        'size': 13.0,
        }
#------------------- Data manipulating ----------------------
#-------------------------  QE  -----------------------------
def get_rid_blank(filename):
    data = []
    file = open(filename, 'r', encoding='utf-8')
    file_data = file.readlines()
    for line in range(0, len(file_data)):
        if file_data[line] == '\n':
            continue
        else:
            file_lines = file_data[line]
            data.append(file_lines)
    return data
# obtain the poscar information
def poscar(filename="POSCAR"):
    poscar = open(filename, 'r', encoding='utf-8')
    poscar_lines = poscar.readlines()

    pos_vector = np.zeros(shape=(3, 3))
    for flag_lines in range(2, 5):
        for i in range(0, 3):
            pos_vector[flag_lines-2][i] = poscar_lines[flag_lines].split()[i]
    coordinates = np.zeros(shape=((len(poscar_lines)-8), 3))
    for flag_lines in range(8, len(poscar_lines)):
        for i in range(0, 3):
            coordinates[flag_lines-8][i] = poscar_lines[flag_lines].split()[i]
    return pos_vector, coordinates

# obtain the information from POSCAR
def get_element(filename = "POSCAR"):
    poscar = open(filename, 'r')
    pos_e = poscar.readlines()
    elements = pos_e[5].lstrip().split()
    element = []*len(elements)
    ele = []
    num_ele = []
    numbers = pos_e[6].lstrip().split()
    for flag_element in range(0, len(elements)):
        ele = elements[flag_element]
        element.append(ele)

    for flag_element_num in range(0, len(numbers)):
        ele = int(numbers[flag_element_num])
        num_ele.append(ele)
    return element, num_ele

def ChaCheng(a1, a2):
    c = []
    c1 = float(a1[1] * a2[2] - a1[2] * a2[1])
    c.append(c1)
    c2 = float(a1[2] * a2[0] - a1[0] * a2[2])
    c.append(c2)
    c3 = float(a1[0] * a2[1] - a1[1] * a2[0])
    c.append(c3)
    return c
def DianCheng(b1, b2):
    d1 = float(b1[0] * b2[0])
    d2 = float(b1[1] * b2[1])
    d3 = float(b1[2] * b2[2])
    d = d1 + d2 + d3
    return d

def get_poscar(filename = "POSCAR"):
    poscar = open(filename, 'r')
    poscar_lines = poscar.readlines()
    pos = []
    for flag_lines in range(2, 5):
        for i in range(0, 3):
            pos.append(float(poscar_lines[flag_lines].split()[i]))
    return pos

def get_reciprocal():
    pos_vector, pos_coordinate = poscar()
    rep = np.zeros(shape=(3, 3))
    a1 = [pos_vector[0][0], pos_vector[0][1], pos_vector[0][2]]
    a2 = [pos_vector[1][0], pos_vector[1][1], pos_vector[1][2]]
    a3 = [pos_vector[2][0], pos_vector[2][1], pos_vector[2][2]]
    volume = DianCheng(a1, ChaCheng(a2, a3))
    scalar = 2 * math.pi / volume
    for i in range(0, 3):
        rep[0][i] = scalar*ChaCheng(a2, a3)[i]
    for j in range(0, 3):
        rep[1][j] = scalar*ChaCheng(a3, a1)[j]
    for k in range(0, 3):
        rep[2][k] = scalar * ChaCheng(a1, a2)[k]
    reciprocal = np.array(rep)
    return reciprocal

# calculate the length of lattice vector
def kpoints(filename, scalar = 0.03):
    reprocal = get_reciprocal()/math.pi
    len_a = np.sqrt(reprocal[0][0]**2 + reprocal[0][1]**2 + reprocal[0][2]**2)
    len_b = np.sqrt(reprocal[1][0]**2 + reprocal[1][1]**2 + reprocal[1][2]**2)
    len_c = np.sqrt(reprocal[2][0]**2 + reprocal[2][1]**2 + reprocal[2][2]**2)
    k_a = math.ceil(len_a / (scalar+0.01))
    k_b = math.ceil(len_b / (scalar+0.01))
    k_c = math.ceil(len_c / (scalar+0.01))
    kpoints = open(filename, "a", encoding='utf-8')
    kpoints.write("\n")
    kpoints.write("K_POINTS {automatic}")
    kpoints.write("\n")
    kpoints.write(str(k_a) + "  " + str(k_b) + "  " + str(k_c))
    kpoints.write("\n")


# POSCAR to the position for QE
def pos_qe(filename, scalar = 0.03):
    element, num_ele = get_element()
    kpoints(filename, scalar = 0.03)
    pos_vector, pos_coordinates = poscar()
    pos = open(filename, "a", encoding='utf-8')
    pos.write("\n")
    pos.write("CELL_PARAMETERS {angstrom}" + "\n")
    for i in range(0, 3):
        for j in range(0, 3):
            pos.write("   " + str(pos_vector[i][j])[0:10].ljust(10, "0") + "   ")
        pos.write("\n")
    pos.write("\n")
    pos.write("ATOMIC_SPECIES" + "\n")
    for flag_type_ele in range(0, len(element)):
        pos.write(str(element[flag_type_ele]) + "\n")
    pos.write("\n")
    pos.write("ATOMIC_POSITIONS {crystal}" + "\n")
    index = 0
    for type_ele in range(0, len(element)):
        for flag_num in range(0, num_ele[type_ele]):
            pos.write(str(element[type_ele]) + "   ")
            for j in range(0, 3):
                pos.write("   " + str(pos_coordinates[index][j])[0:10].ljust(10, "0") + "   ")
            index = index + 1
            pos.write("\n")
    pos.write("\n")

def manipulate_qe():
    print(" **************************** What kind of calculations do you want ? ******************************* ")
    print("(1) vc-relax " + "\t" + "(2) relax ")
    print("(3) scf      " + "\t" + "(4) bands ")
    print("(5) dos/nscf " + "\t" + "(6) md    ")
    print("(7) vc-md ")
    type_cal = int(input("Input a Number: "))
    if type_cal == 1:
        os.system("cp "+file_path+"/vc-relax ./")
        pos_qe("vc-relax")
    elif type_cal == 2:
        os.system("cp "+file_path+"/relax ./")
        pos_qe("relax")
    elif type_cal == 3:
        os.system("cp "+file_path+"/scf ./")
        pos_qe("scf")
    elif type_cal == 4:
        os.system("cp "+file_path+"/band ./")
        pos_qe("band")
    elif type_cal == 5:
        os.system("cp "+file_path+"/nscf ./")
        pos_qe("nscf")
    elif type_cal == 6:
        print("Waiting")
    elif type_cal == 7:
        print("Waiting")
    else:
        print("You have input a wrong number !")

if __name__ == '__main__':
    manipulate_qe()

