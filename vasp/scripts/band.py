# -*- coding: utf-8 -*-
"""
Created on 10:06 13-11-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: band.py
"""
import os, sys
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)
import numpy as np
from matplotlib import pyplot as plt
import math
import initial

# ******************** band ************************
# get K-PATH form KPOINTS and calculating the kpt
def get_kpt(filename = "KPOINTS"):
    high_kpt_num = 0
    high_kpt_point = []
    kpt = []
    kpath = []
    kpoints = open(filename, 'r').read().splitlines()
    num_kpt_line = int(kpoints[1].split()[0])
    for flag_kpath in range(4, len(kpoints), 3):
        kpt_data1 = list(map(float, kpoints[flag_kpath].split()[0:3]))
        kpt_points_data = kpoints[flag_kpath].split('!')[1].strip()
        if kpt_points_data == "\Gamma":              # transform the gamma character
            kpt_points_data = u"Γ"
        high_kpt_num = high_kpt_num + 1              # get high symmetry point number
        high_kpt_point.append(kpt_points_data)       # get high symmetry points
        kpath.append(kpt_data1[0:3])
    kpt_points_data1 = kpoints[-1].split('!')[1].strip()
    if kpt_points_data1 == "\Gamma":             # transform the gamma character
        kpt_points_data1 = u"Γ"
    high_kpt_point.append(kpt_points_data1)      # get high symmetry points
    kpt_data2 = list(map(float, kpoints[-1].split()[0:3]))
    kpath.append(kpt_data2)                          # get high symmetry path
    high_kpt_num = int(high_kpt_num + 1)                 # get high symmetry point number

    # calculating all the kpoints for the k axis of band
    # 1. transform the high symmetry points to its real coordinates in k space
    k_real = []
    k_middle = []
    k = [0, 0, 0]
    kx = [0, 0, 0]
    ky = [0, 0, 0]
    kz = [0, 0, 0]
    for flag_k in range(0, len(kpath)):
        for flag_d in range(0, 3):
            kx[flag_d] = kpath[flag_k][0] * initial.get_reciprocal()[flag_d]
            ky[flag_d] = kpath[flag_k][1] * initial.get_reciprocal()[flag_d+3]
            kz[flag_d] = kpath[flag_k][2] * initial.get_reciprocal()[flag_d+6]
            k[flag_d] = kx[flag_d] + ky[flag_d] + kz[flag_d]
        k_middle = [k[0], k[1], k[2]]                   # real high symmetry point coordinations
        k_real.append(k_middle)

    # calculating the kpath_length
    length = [0] * len(k_real)
    for cal_k in range(0, len(k_real)):
        if cal_k < len(k_real)-1:
            for i in range(0, 3):
                length[cal_k] = (
                    (k_real[cal_k + 1][0] - k_real[cal_k][0]) ** 2 + (k_real[cal_k + 1][1] - k_real[cal_k][1]) ** 2 + (
                                k_real[cal_k + 1][2] - k_real[cal_k][2]) ** 2) ** (0.5)
        else:
            break

    # get all the kpoints
    num_kpt = (num_kpt_line) * (high_kpt_num-1)
    kpt = []
    kn = 0.0
    lent = 0
    for flag_kpt_hn in range(0, high_kpt_num-1):
        for flag_kpt_pn in range(0, num_kpt_line):
            kn = lent + (length[flag_kpt_hn] / (num_kpt_line -1)) * flag_kpt_pn
            round(float(kn), 6)
            kpt.append(kn)
        lent = lent + length[flag_kpt_hn]
    return int(num_kpt), length, high_kpt_point, kpt

# output klines.dat
def output_klines(filename = "KPOINTS"):
    num_kpt, length, high_kpt_point, kpt = get_kpt(filename)
    len_high_k = 0
    for j in range(0, len(length)):
        length[j] = length[j] + len_high_k
        len_high_k = length[j]
    with open("klines.dat", "w", encoding='utf-8') as kline:
        kline.write("high_kpoints" + "\n")
        kline.write(str("0.").ljust(10, '0') + "  " + str(10) + "\n")
        kline.write(str("0.").ljust(10, '0') + "  " + str(-15) + "\n")
        for l_num in range(0, len(length) - 1):
            kline.write(str(length[l_num])[:10].ljust(10, '0') + "  " + str(-15) + "\n")
            kline.write(str(length[l_num])[:10].ljust(10, '0') + "  " + str(10) + "\n")
            kline.write(str(length[l_num])[:10].ljust(10, '0') + "  " + str(-15) + "\n")
    kline.close()



# get Energy level form Eigenval file
def get_eigenval():
    num_kpt, length, high_kpt_point, kpt = get_kpt()
    fermi = initial.get_fermiLevel()
    eig = []
    eig_up = []
    eig_dw = []
    energy_kpt = []
    energy_kpt_up = []
    energy_kpt_dw = []
#    energy = []
    nbd = initial.get_NBANDS() + 2
    ispin = initial.get_ispin()
    eigenval_lines = open('EIGENVAL', 'r').read().splitlines()
    if ispin == 1:
        for flag_nbd_num in range(0, nbd-2):
            for flag_k_num in range(0, int(num_kpt)+1):
                if flag_k_num < int(num_kpt):
                    flag_orbital = 8 + flag_nbd_num + flag_k_num * nbd
                    eig = float(eigenval_lines[flag_orbital].split()[1]) - fermi
                    round(float(eig), 6)
                    energy_kpt.append(eig)
        return energy_kpt
    elif ispin == 2:
        for flag_nbd_num in range(0, nbd-2):
            for flag_k_num in range(0, int(num_kpt)+1):
                if flag_k_num < int(num_kpt):
                    flag_orbital = 8 + flag_nbd_num + flag_k_num * nbd
                    eig_up = float(eigenval_lines[flag_orbital].split()[1]) - fermi
                    eig_dw = float(eigenval_lines[flag_orbital].split()[2]) - fermi
                    round(float(eig_up), 6)
                    round(float(eig_dw), 6)
                    energy_kpt_up.append(eig_up)
                    energy_kpt_dw.append(eig_dw)
        return energy_kpt_up, energy_kpt_dw
    else:
        return "No Ispin Information Obtained !"


# output band.dat
def output_band():
    output_klines()
    ispin = initial.get_ispin()
    if ispin == 1:
        num_kpt, length, high_kpt_point, kpoints = get_kpt()
        nbd = initial.get_NBANDS()
        kpt = np.transpose(np.array(kpoints))
        # kpt = kpoints.astype(np.float)
        eng = np.array(get_eigenval())
        energy = eng.reshape(nbd, num_kpt)
        # reverse the sequence
        #    rev_kpt = kpt[::-1]
        with open("band.dat", "w", encoding='utf-8') as band:
            band.write("Band Data" + "\n")
            for k_num in range(0, num_kpt):
                band.write(str(kpt[k_num])[:10].ljust(10, ' ') + "  ")
                for flag_nbd in range(0, nbd):
                    band.write(str(energy[flag_nbd][k_num])[:10].ljust(10, ' ') + "  ")
                band.write("\n")
        band.close()
    elif ispin == 2:
        num_kpt, length, high_kpt_point, kpoints = get_kpt()
        nbd = initial.get_NBANDS()
        kpt = np.transpose(np.array(kpoints))
        eng_up, eng_dw = get_eigenval()
        eig_up = np.array(eng_up)
        eig_dw = np.array(eng_dw)
        energy_up = eig_up.reshape(nbd, num_kpt)
        energy_dw = eig_dw.reshape(nbd, num_kpt)
        with open("band_up.dat", "w", encoding='utf-8') as band_up:
            band_up.write("Band Data" + "\n")
            for k_num in range(0, num_kpt):
                band_up.write(str(kpt[k_num])[:10].ljust(10, ' ') + "    ")
                for flag_nbd in range(0, nbd):
                    band_up.write(str(energy_up[flag_nbd][k_num])[:10].ljust(10, ' ') + "    ")
                band_up.write("\n")
        band_up.close()
        with open("band_dw.dat", "w", encoding='utf-8') as band_dw:
            band_dw.write("Band Data" + "\n")
            for k_num in range(0, num_kpt):
                band_dw.write(str(kpt[k_num])[:10].ljust(10, ' ') + "    ")
                for flag_nbd in range(0, nbd):
                    band_dw.write(str(energy_dw[flag_nbd][k_num])[:10].ljust(10, ' ') + "    ")
                band_dw.write("\n")
        band_dw.close()
    else:
        print("No Ispin Information Obtained !")


# plot band
def bandplot(eng_r = [-10, 6, 2], color_noispin = "black", color_ispin = "red"):
    output_band()
    ispin = initial.get_ispin()
    if ispin == 1:
        num_kpt, length, high_kpt_point, kpoints = get_kpt()
        nbd = initial.get_NBANDS()
        banddata = np.loadtxt("band.dat", skiprows=1, dtype=np.float64)
        kpt = banddata[:, 0]
        len_high_k = 0
        for j in range(0, len(length) - 1):
            length[j] = length[j] + len_high_k
            len_high_k = length[j]
        fig, ax = plt.subplots()
        for i in range(1, nbd + 1):
            plt.plot(kpt, banddata[:, i], color=color_noispin)
        plt.xlim(0, max(length))
        x_group_label = []
        x_start = 0
        x_group_label.append(x_start)
        for k in range(0, len(length) - 1):
            x_group_label.append(length[k])
        plt.xticks(x_group_label)
        ax.set_xticklabels(high_kpt_point, rotation=0, fontsize=12, fontname='arial')
        # plt.yticks(eng_r[0], eng_r[1], eng_r[2])
        plt.ylim(eng_r[0], eng_r[1])
        plt.ylabel("E - " + "E$_\mathrm{{F}}$" + " (eV)", fontsize=14, fontname='arial')
        ytick = np.arange(eng_r[0], eng_r[1] + 1, eng_r[2])
        a = int(len(ytick) / 2)
        plt.yticks(np.insert(ytick, a, 0), fontsize=12, )

        ax.axhline(y=0, xmin=0, xmax=1, linestyle='--', linewidth=0.5, color='0.5')
        for i in length[0:-1]:
            ax.axvline(x=i, ymin=0, ymax=1, linestyle='--', linewidth=0.5, color='0.5')
        plt.savefig("band.png", dpi=300)
        plt.show()
    elif ispin == 2:
        num_kpt, length, high_kpt_point, kpoints = get_kpt()
        nbd = initial.get_NBANDS()
        banddata_up = np.loadtxt("band_up.dat", skiprows=1, dtype=np.float64)
        banddata_dw = np.loadtxt("band_dw.dat", skiprows=1, dtype=np.float64)
        kpt = banddata_up[:, 0]
        len_high_k = 0
        for j in range(0, len(length) - 1):
            length[j] = length[j] + len_high_k
            len_high_k = length[j]
        fig, ax = plt.subplots()
        for i in range(1, nbd + 1):
            plt.plot(kpt, banddata_up[:, i], color=color_noispin)
            plt.plot(kpt, banddata_dw[:, i], color=color_ispin)
        plt.xlim(0, max(length))
        ax.legend(("Spin_UP", "Spin_DW"), loc='best')
        x_group_label = []
        x_start = 0
        x_group_label.append(x_start)
        for k in range(0, len(length) - 1):
            x_group_label.append(length[k])
        plt.xticks(x_group_label)
        ax.set_xticklabels(high_kpt_point, rotation=0, fontsize=12, fontname='arial')
        # plt.yticks(eng_r[0], eng_r[1], eng_r[2])
        plt.ylim(eng_r[0], eng_r[1])
        plt.ylabel("E - " + "E$_\mathrm{{F}}$" + " (eV)", fontsize=14, fontname='arial')
        ytick = np.arange(eng_r[0], eng_r[1] + 1, eng_r[2])
        a = int(len(ytick) / 2)
        plt.yticks(np.insert(ytick, a, 0), fontsize=12, )

        ax.axhline(y=0, xmin=0, xmax=1, linestyle='--', linewidth=0.5, color='0.5')
        for i in length[0:-1]:
            ax.axvline(x=i, ymin=0, ymax=1, linestyle='--', linewidth=0.5, color='0.5')
        plt.savefig("band.png", dpi=300)
        plt.show()
    else:
        print("No Ispin Information Obtained !")


# manipulate the fig
def manipulate_bandplot():
    print("    ****************************************************************")
    print("    *This is a code used to plot kinds of band structure,written by*")
    print("    *                          XY Ding                             *")
    print("    ****************************************************************")
    print("\n")
    print("                       (^o^)GOOD LUCK!(^o^)                         ")
    print("\n")
    ispin = initial.get_ispin()
    if ispin == 1:
        print("    *************************** ISPIN is 1 ****************************")
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
                bandplot(energy, color)
            else:
                bandplot(energy, "black")
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
                bandplot(energy, color)
            else:
                bandplot([-10, 6, 2], color)
        else:
            bandplot([-10, 6, 2], "black")
    elif ispin == 2:
        print("    *************************** ISPIN is 2 ****************************")
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
                color_up = str(input("color for Spin_UP:"))
                color_dw = str(input("color for Spin_DW:"))
                bandplot(energy, color_up, color_dw)
            else:
                bandplot(energy, "black", "red")
        elif inint == 2:
            color_up = str(input("color for Spin_UP:"))
            color_dw = str(input("color for Spin_DW:"))
            print("(1) continue to setting energy range")
            print("(2) end setting")
            con_end = int(input("Input number:"))
            if con_end == 1:
                min_energy = float(input("minimum energy:"))
                max_energy = float(input("maximum energy:"))
                scale_energy = float(input("energy scale:"))
                energy = [min_energy, max_energy, scale_energy]
                bandplot(energy, color_up, color_dw)
            else:
                bandplot([-10, 6, 2], color_up, color_dw)
        else:
            bandplot([-10, 6, 2], "black", "red")
    else:
        print("    ******************** ISPIN is not 1 or 2!  Good bye ! ********************")

if __name__ == '__main__':
    manipulate_bandplot()
