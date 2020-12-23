# -*- coding: utf-8 -*-
"""
Created on 19:54 04-11-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: pvasp.py
"""
import numpy as np
from matplotlib import pyplot as plt
import os, sys
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)
import math
import get_data as gd

#------------------- Data manipulating ----------------------
#-------------------------  QE  -----------------------------
# get rid of blank line
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

def get_kpath(filename = "KPATH"):
    high_kpt_point_char = []
    kpt = []
    kpath = []
    kpoints = get_rid_blank(filename)
    num_high_k = int(kpoints[1].split()[0])
    num_kline=int(kpoints[2].split()[3])
    high_kpt_points = np.zeros(shape=(num_high_k, 3))
    for flag_kpath in range(2, num_high_k+2):
        kpt_data1 = list(map(float, kpoints[flag_kpath].split()[0:3]))
        kpt_points_data = kpoints[flag_kpath].split('!')[1].strip()
        if kpt_points_data == "\Gamma":              # transform the gamma character
            kpt_points_data = u"Γ"
        high_kpt_point_char.append(kpt_points_data)      # get high symmetry points
        kpath.append(kpt_data1[0:3])
    for i in range(0, num_high_k):
        for j in range(0, 3):
            high_kpt_points[i][j] = kpoints[i+2].split('!')[0].split()[j]
    #print(num_high_k, num_kline, high_kpt_point, kpath)
    return num_high_k, num_kline, high_kpt_point_char, high_kpt_points

# get band data
# get nbnd
def get_nbnd(filename = "band.dat"):
    data = open(filename, 'r')
    data_band = data.readlines()
    j = 0
    band_lines = []
    length_lines = 0
    length_lines0 = 0
    for lines in range(2, len(data_band)):
        if len(data_band[lines].split()) == 3 and data_band[lines].split()[0][0] == '0':
           break
        else:
            band_lines.append(data_band[lines].split())
        length_lines = len(band_lines[lines-2]) + length_lines0
        length_lines0 = length_lines
    return length_lines

# get k_coordinate and band level
def get_band(filename = "band.dat"):
    nbnd = get_nbnd(filename)
    data = open(filename, 'r')
    data_band = data.readlines()
    nk = int(data_band[0].split('=')[2].split()[0])
    k_coordinate = []
    len_banddata = math.ceil(nbnd/10)
    band_index = []
    for i in range(0, nk):
        for j in range(0, len_banddata):
            for k in range(0, len(data_band[2 + i * (len_banddata+1) + j].split())):
                band_index.append(data_band[2 + i * (len_banddata+1) + j].split()[k])
        k_coordinate.append(data_band[1 + i * (len_banddata+1)].split())
    k_coordinate = np.array(k_coordinate).reshape(nk, 3)
    eigenval = np.array(band_index).reshape(nk, nbnd)
    return k_coordinate, eigenval

# calculate k_length
def cal_k_length(k_coordinates):
    k, eigenval = get_band("band.dat")
    k_coordinates = k_coordinates.astype(float)
    eigenval = eigenval.astype(float)
    k_coordinate = np.zeros(shape=(len(k_coordinates), 3))
    for i in range(0, len(k_coordinates)):
        k_coordinate[i][0] = k_coordinates[i][0] * gd.get_reciprocal()[0] + k_coordinates[i][1] * \
                                  gd.get_reciprocal()[3] + k_coordinates[i][2]*gd.get_reciprocal()[6]
        k_coordinate[i][1] = k_coordinates[i][0] * gd.get_reciprocal()[1] + k_coordinates[i][1] * \
                                  gd.get_reciprocal()[4] + k_coordinates[i][2]*gd.get_reciprocal()[7]
        k_coordinate[i][2] = k_coordinates[i][0] * gd.get_reciprocal()[2] + k_coordinates[i][1] * \
                                  gd.get_reciprocal()[5] + k_coordinates[i][2]*gd.get_reciprocal()[8]
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

# output klines.dat
def output_klines():
    num_high_k, num_kline, high_kpt_point_char, high_kpt_points = get_kpath()
    high_kpt_point = np.array(high_kpt_points).astype(float)
    length = cal_k_length(high_kpt_point)
    with open("klines.dat", "w", encoding='utf-8') as kline:
        kline.write("high_kpoints" + "\n")
        kline.write(str("0.").ljust(10, '0') + "  " + str(10) + "\n")
        kline.write(str("0.").ljust(10, '0') + "  " + str(-15) + "\n")
        for l_num in range(0, len(length)-1):
            kline.write(str(length[l_num])[:10].ljust(10, '0') + "  " + str(-15) + "\n")
            kline.write(str(length[l_num])[:10].ljust(10, '0') + "  " + str(10) + "\n")
            kline.write(str(length[l_num])[:10].ljust(10, '0') + "  " + str(-15) + "\n")
    kline.close()

# output data of band
def output_band(filename_band="band.dat", filename_scf="scf.out"):
    k_coordinate, eigenval = get_band(filename_band)
    k_length = cal_k_length(k_coordinate)
    eigenval = eigenval.astype(float)
    efermi = gd.get_fermiLevel(filename_scf)
    nbnd = len(eigenval[0])
    output_klines()
    with open("band_plot.dat", "w", encoding='utf-8') as band:
        band.write("band.dat" + "\n")
        for i in range(0, len(eigenval)):
            band.write(str(k_length[i])[:6].ljust(10, ' '))
            for j in range(0, nbnd):
                band.write(str(eigenval[i][j]-efermi)[:6].ljust(10, ' '))
            band.write("\n")
    band.close()

def bandplot(filename = ["relax.out", "band", "band.dat", "scf.out"], eng_r = [-10, 6, 2], color_noispin = "blue", color_ispin = "red"):
    gd.trans_qe_kpath(filename[1])
    gd.trans_qe_pos(filename[0])
    output_band(filename[2], filename[3])
    nbnd = get_nbnd(filename[2])
    num_high_k, num_kline, high_kpt_point_char, high_kpt_points = get_kpath("KPATH")
    length_k = np.zeros(shape=(num_high_k))
    band_data = np.loadtxt("band_plot.dat", skiprows=1, dtype=np.float64)
    klines_data = np.loadtxt("klines.dat", skiprows=1, dtype=np.float64)
    for i in range(0, num_high_k):
        length_k[i] = band_data[i*num_kline, 0]
    fig, ax = plt.subplots()
    for i in range(0, nbnd):
        plt.plot(band_data[:, 0], band_data[:, i+1], color=color_noispin)
    plt.xlim(min(length_k), max(length_k))
    plt.xticks(length_k)
    ax.set_xticklabels(high_kpt_point_char, rotation=0, fontsize=12, fontname='arial')
    plt.ylim(eng_r[0], eng_r[1])
    plt.ylabel("E - " + "E$_\mathrm{{F}}$" + " (eV)", fontsize=14, fontname='arial')
    ytick = np.arange(eng_r[0], eng_r[1] + 1, eng_r[2])
    a = int(len(ytick) / 2)
    plt.yticks(np.insert(ytick, a, 0), fontsize=12, )
    ax.axhline(y=0, xmin=0, xmax=1, linestyle='--', linewidth=0.5, color='0.5')
    for i in length_k[0:-1]:
        ax.axvline(x=i, ymin=0, ymax=1, linestyle='--', linewidth=0.5, color='0.5')
    plt.savefig("band.png", dpi=300)
    plt.show()
# manipulate band plot
def manipulate_bandplot():
    print("    ****************************************************************")
    print("    *This is a code used to plot kinds of band structure,written by*")
    print("    *                          XY Ding                             *")
    print("    ****************************************************************")
    print("\n")
    print("                       (^o^)GOOD LUCK!(^o^)                         ")
    print("\n")
    print("    ********************** Entering band plot ***********************")
    print("Whether the output file are relax.out, band, band.dat, scf.out or not !")
    print("The default output file for: relax.out, band for pwscf" 
          "band, scf and high_kpoint_path are band.dat and scf.out.")
    print("(1) yes")
    print("(2) no, let me re-define")
    initial_input = int(input("Input number"))
    if initial_input == 1:
        filename = ["relax.out", "band", "band.dat", "scf.out"]
        print("(1) energy range")
        print("(2) color")
        print("(3) Use default setting")
        inint = int(input("Input number:"))
        if inint == 1:
            energy = []
            min_energy = float(input("minimum energy:"))
            max_energy = float(input("maximum energy:"))
            scale_energy = float(input("energy scale:"))
            energy = [min_energy, max_energy, scale_energy]
            print("(1) continue to setting color")
            print("(2) end setting")
            con_end = int(input("Input number:"))
            if con_end == 1:
                color = str(input("which color do you want:"))
                bandplot(filename, energy, color)
            else:
                bandplot(filename, energy, "black")
        elif inint == 2:
            color = str(input("which color do you want:"))
            print("(1) continue to setting energy range")
            print("(2) end setting")
            con_end = int(input("Input number:"))
            if con_end == 1:
                min_energy = float(input("minimum energy:"))
                max_energy = float(input("maximum energy:"))
                scale_energy = float(input("energy scale:"))
                energy = [min_energy, max_energy, scale_energy]
                bandplot(filename, energy, color)
            else:
                bandplot(filename, [-10, 6, 2], color)
        else:
            bandplot(filename, [-10, 6, 2], "black")
    elif initial_input == 2:
        print("************************ Please input 4 filenames for reading ****************************")
        print("*************************** Split each filename with blank *******************************")
        print("******************* Default is: relax.out, band, band.dat, scf.out ***********************")
        filename = input("Input the 4 filenames: ").split()
        if len(filename) == 4:
            print("(1) energy range")
            print("(2) color")
            print("(3) Use default setting")
            inint = int(input("Input number:"))
            if inint == 1:
                energy = []
                min_energy = float(input("minimum energy:"))
                max_energy = float(input("maximum energy:"))
                scale_energy = float(input("energy scale:"))
                energy = [min_energy, max_energy, scale_energy]
                print("(1) continue to setting color")
                print("(2) end setting")
                con_end = int(input("Input number:"))
                if con_end == 1:
                    color = str(input("which color do you want:"))
                    bandplot(filename, energy, color)
                else:
                    bandplot(filename, energy, "black")
            elif inint == 2:
                color = str(input("which color do you want:"))
                print("(1) continue to setting energy range")
                print("(2) end setting")
                con_end = int(input("Input number:"))
                if con_end == 1:
                    min_energy = float(input("minimum energy:"))
                    max_energy = float(input("maximum energy:"))
                    scale_energy = float(input("energy scale:"))
                    energy = [min_energy, max_energy, scale_energy]
                    bandplot(filename, energy, color)
                else:
                    bandplot(filename, [-10, 6, 2], color)
            else:
                bandplot(filename, [-10, 6, 2], "black")
        else:
            print("You have input the wrong file number !")

if __name__ == '__main__':
    manipulate_bandplot()

