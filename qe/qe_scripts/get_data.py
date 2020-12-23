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

# get the spin information from INCAR file
def get_ispin(filename = "band"):
    incar = open(filename,'r')
    incar_lines = incar.readlines()
    flag=[]
    for flag_incar in range(0,len(incar_lines)):
        if 'ISPIN' in incar_lines[flag_incar]:
            flag = incar_lines[flag_incar].lstrip()
            spinline = incar_lines[flag_incar].split('=')[1]
            if flag[0]=='#':          # whether the ISPIN is commented or not
                    ispin = 1
            else:
                    ispin = int(spinline.split()[0])
            break
        else:
                ispin = 1
    return int(ispin)

# get the fermi level from OUTCAR
def get_fermiLevel(filename):
    outcar = open(filename,'r')
    outcar_lines = outcar.readlines()
    efermi = 0
    for flag_outcar in range(0,len(outcar_lines)):
        if 'Fermi' in outcar_lines[flag_outcar]:
            spinline = outcar_lines[flag_outcar].split()[4]
            efermi = float(spinline.split()[0])
            break
    return efermi
# get reciprocal unit cell
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

def get_poscar():
    poscar = open('POSCAR', 'r')
    poscar_lines = poscar.readlines()
    pos = []
    for flag_lines in range(2, 5):
        for i in range(0, 3):
            pos.append(float(poscar_lines[flag_lines].split()[i]))
    return pos

def get_reciprocal():
    pos = get_poscar()
    rep = []
    a1 = [pos[0], pos[1], pos[2]]
    a2 = [pos[3], pos[4], pos[5]]
    a3 = [pos[6], pos[7], pos[8]]
    volume = DianCheng(a1, ChaCheng(a2, a3))
    scalar = 2 * math.pi / volume
    for i in range(0, 3):
        b1 = scalar*ChaCheng(a2, a3)[i]
        rep.append(b1)
    for j in range(0, 3):
        b2 = scalar*ChaCheng(a3, a1)[j]
        rep.append(b2)
    for k in range(0, 3):
        b3 = scalar * ChaCheng(a1, a2)[k]
        rep.append(b3)
    return rep

# calculate k_length
def cal_k_length():
    k_coordinates, eigenval = get_band("band.dat")
    k_coordinates = k_coordinates.astype(float)
    eigenval = eigenval.astype(float)
    k_coordinate = np.zeros(shape=(len(k_coordinates), 3))
    for i in range(0, len(k_coordinates)):
        k_coordinate[i][0] = k_coordinates[i][0] * get_reciprocal()[0] + k_coordinates[i][1] * \
                                  get_reciprocal()[3] + k_coordinates[i][2]*get_reciprocal()[6]
        k_coordinate[i][1] = k_coordinates[i][0] * get_reciprocal()[1] + k_coordinates[i][1] * \
                                  get_reciprocal()[4] + k_coordinates[i][2]*get_reciprocal()[7]
        k_coordinate[i][2] = k_coordinates[i][0] * get_reciprocal()[2] + k_coordinates[i][1] * \
                                  get_reciprocal()[5] + k_coordinates[i][2]*get_reciprocal()[8]
    length = np.zeros(shape=(len(k_coordinate)))
    t, g, k = 0, 0, 0
    leng_pre = 0
    for flag_num in range(0, len(k_coordinate)):
        leng = np.sqrt((k_coordinate[flag_num][0] - t)**2 + (k_coordinate[flag_num][1] - g)**2
                       + (k_coordinate[flag_num][2]- k)**2)
        t = k_coordinate[flag_num][0]
        g = k_coordinate[flag_num][1]
        k = k_coordinate[flag_num][2]
        length[flag_num] = leng + leng_pre
        leng_pre = length[flag_num]
    return length

# transform the qe output file of relax to POSCAR
# get number of atoms
def get_natom(filename = "relax.out"):
    relax = open(filename, 'r')
    data = relax.readlines()
    natom = 0
    for i in range(0, len(data)):
        if "number of atoms/cell" in data[i]:
            natom = data[i].split('=')[1]
            break
    return int(natom)
# count atom type number
def grep_relax_type(filename = "relax"):
    relax = open(filename, 'r')
    data = relax.readlines()
    relax_type = "vc-relax"
    for i in range(0, len(data)):
        if "calculation" in data[i]:
             relax_type = data[i].split('=')[1]
    relax_type = relax_type.split("\"")[1]
    return relax_type

def get_pos_qe_vc(filename = "relax.out"):
    natom = get_natom(filename)
    relax = open(filename, 'r')
    pos = relax.readlines()
    cell = []
    pos_coordinate = []
    for i in range(0, len(pos)):
        if 'Begin final coordinates' in pos[i]:
            for j in range(0, 3):
                cell.append(pos[i + 5 + j].split())
            for k in range(0, natom):
                pos_coordinate.append(pos[i + 10 + k].split())
        elif 'End final coordinates' in pos[i]:
            break
    cell = np.array(cell).reshape(3, 3).astype(float)
    pos_coordinate = np.array(pos_coordinate).reshape(natom, 4)
    coordinate = np.zeros(shape=(natom, 3))
    for i in range(0, natom):
        for j in range(0, 3):
            coordinate[i][j] = pos_coordinate[i][j+1]
    atom_char = pos_coordinate[:, 0]
    atom_type = pos_coordinate[0, 0]
    cou_num = 0
    atom_type_num = []
    atom_type_char = []
    for i in range(0, natom):
        if atom_char[i] == atom_type:
            cou_num = cou_num + 1
        else:
            atom_type = atom_char[i]
            atom_type_num.append(cou_num)
            cou_num = 1
            continue
    cou_num = natom - int(sum(atom_type_num))
    atom_type_num.append(cou_num)
    type_num = 0
    for i in range(0, len(atom_type_num)):
        atom_type_char.append(atom_char[type_num])
        type_num = atom_type_num[i] + type_num
    return atom_type_num, atom_type_char, cell, coordinate
# get pos of qe relaxation
def get_pos_qe(filename = "relax.out"):
    relax_type = grep_relax_type()
    if relax_type == "vc-relax":
        atom_type_num, atom_type_char, cell, coordinate = get_pos_qe_vc(filename)
        return atom_type_num, atom_type_char, cell, coordinate
    elif relax_type == "relax":
        cell = []
        natom = get_natom(filename)
        relax_input = open("relax")
        input_r = relax_input.readlines()
        for i in range(0, len(input_r)):
            if "CELL_PARAMETERS" in input_r[i]:
                for k in range(0, 3):
                    cell.append(input_r[i+k+1].split())
        cell = np.array(cell).reshape(3, 3).astype(float)
        relax = open(filename, 'r')
        pos = relax.readlines()
        pos_coordinate = []
        for i in range(0, len(pos)):
            if 'Begin final coordinates' in pos[i]:
                for k in range(0, natom):
                    pos_coordinate.append(pos[i + 3 + k].split())
            elif 'End final coordinates' in pos[i]:
                break
        pos_coordinate = np.array(pos_coordinate).reshape(natom, 4)
        coordinate = np.zeros(shape=(natom, 3))
        for i in range(0, natom):
            for j in range(0, 3):
                coordinate[i][j] = pos_coordinate[i][j + 1]
        atom_char = pos_coordinate[:, 0]
        atom_type = pos_coordinate[0, 0]
        cou_num = 0
        atom_type_num = []
        atom_type_char = []
        for i in range(0, natom):
            if atom_char[i] == atom_type:
                cou_num = cou_num + 1
            else:
                atom_type = atom_char[i]
                atom_type_num.append(cou_num)
                cou_num = 1
                continue
        cou_num = natom - int(sum(atom_type_num))
        atom_type_num.append(cou_num)
        type_num = 0
        for i in range(0, len(atom_type_num)):
            atom_type_char.append(atom_char[type_num])
            type_num = atom_type_num[i] + type_num
        return atom_type_num, atom_type_char, cell, coordinate
# output poscar
def trans_qe_pos(filename = "relax.out"):
    atom_type_num, atom_type_char, cell, coordinate = get_pos_qe(filename)
    natom = sum(atom_type_num)
    num_type = len(atom_type_num)
    with open("CONTCAR", "w", encoding='utf-8') as pos:
        pos.write("CELL PARAMETERS" + '\n')
        pos.write("   " + str(1.0))
        pos.write('\n')
        for i in range(0, 3):
            pos.write("   ")
            for j in range(0, 3):
                pos.write(str(cell[i][j])[:10].ljust(10, '0') + "   ")
            pos.write('\n')
        for i in range(0, num_type):
            pos.write("  " + str(atom_type_char[i]))
        pos.write('\n')
        for j in range(0, num_type):
            pos.write("  " + str(atom_type_num[j]))
        pos.write('\n')
        pos.write("Direct")
        pos.write('\n')
        for h in range(0, natom):
            pos.write("   ")
            for k in range(0, 3):
                pos.write(str(coordinate[h][k])[:10].ljust(10, '0') + "   ")
            pos.write('\n')
    pos.close()

# transform the high-symmetry path of QE to KPATH
def trans_qe_kpath(filename = "band"):
    data = open(filename, 'r')
    band = data.readlines()
    num_high = 0
    high_k_coordinate = []
    for i in range(0, len(band)):
        if "K_POINTS" in band[i]:
            num_high = int(band[i+1])
            for j in range(0, num_high):
                high_k_coordinate.append(band[i+2+j])
    with open("KPATH", "w", encoding='utf-8') as kpath:
        kpath.write("K_POINTS {crystal_b}" + '\n')
        kpath.write(str(num_high) + '\n')
        for i in range(0, num_high):
            kpath.write(high_k_coordinate[i])
    kpath.close()

if __name__ == '__main__':
    trans_qe_pos()


