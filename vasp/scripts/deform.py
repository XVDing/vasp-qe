# -*- coding: utf-8 -*-
"""
Created on 10:14 13-11-2020 

@author: Xian-Yong Ding
mail to: dxy_vasp@163.com
python3: get_data.py
"""
import numpy as np
from matplotlib import pyplot as plt
import math
import sys, os
curPath_soft = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath_soft)

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

# deform
def deform(c, filename = "POSCAR_orig"):
    pos, pos_coordinates = get_poscar(filename)
    poscar = np.array(pos).reshape(3, 3)
    pos_matrix = np.matrix(poscar)
    unit_matrix = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    ebtheta = np.matrix([[c[0], c[5]/2.0, c[4]/2.0], [c[5]/2.0, c[1], c[3]/2.0], [c[4]/2.0, c[3]/2.0, c[2]]])
    result = pos_matrix * (unit_matrix + ebtheta)
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

# manipulate deform
def manipulate_deform():
    os.system("rm -r c11")
    os.system("rm -r c22")
    os.system("rm -r c12")
    os.system("rm -r c66")
    c = [0, 0, 0, 0, 0, 0]
    deform(c,"POSCAR") # reform the POSCAR file
    os.system("cp POSCAR POSCAR_orig")
    delta = 0.005
    file_deform = ["c11", "c22", "c12", "c66"]
    for i in range(0, 4):
        os.system("mkdir " + str(file_deform[i]))
        for j in range(-6, 7):
            if i == 0:
                d = [delta * j, 0, 0, 0, 0, 0]
                deform(d)
            elif i == 1:
                d = [0, delta * j, 0, 0, 0, 0]
                deform(d)
            elif i == 2:
                d = [delta * j, delta * j, 0, 0, 0, 0]
                deform(d)
            elif i == 3:
                d = [0, 0, 0, 0, 0, delta * j]
                deform(d)
            else:
                print("Something wrong happens !")
            os.system("mkdir " + str(file_deform[i]) + "/" + str('strain') + str(delta * j))
            os.system("cp "+ curPath + "/POSCAR " + str(file_deform[i]) + "/" + str('strain') + str(delta * j))
            os.system("cp "+ curPath + "/INCAR " + str(file_deform[i]) + "/" + str('strain') + str(delta * j))
            os.system("cp "+ curPath + "/POTCAR " + str(file_deform[i]) + "/" + str('strain') + str(delta * j))
            os.system("cp "+ curPath +"/KPOINTS " + str(file_deform[i]) + "/" + str('strain') + str(delta * j))
    os.system("cp " + curPath_soft + "/submit_deform.sh " + curPath + "/")
    os.system("cp " + curPath + "/POSCAR_orig " + curPath + "/POSCAR")


if __name__ == '__main__':
    manipulate_deform()




