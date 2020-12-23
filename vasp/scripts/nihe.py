# -*- coding: utf-8 -*-
"""
Created on 6:07 PM 11-12-2020 

@author: X-Y Ding
mail to: dxy_vasp@163.com
python3: nihe.py
"""
import os, sys
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)
import numpy as np
from matplotlib import pyplot as plt

# get the fermi level from OUTCAR
def get_fermiLevel():
    outcar = open('OUTCAR','r')
    outcar_lines = outcar.readlines()
    for flag_outcar in range(0,len(outcar_lines)):
        if 'E-fermi' in outcar_lines[flag_outcar]:
            spinline = outcar_lines[flag_outcar].split(':')[1]
            efermi = float(spinline.split()[0])
            break
    return efermi

# get NBANDS from OUTCAR
def get_NBANDS():
    wann = open('wannier90.win','r')
    wann_lines = wann.readlines()
    nbands = 0
    for flag_wann in range(0,len(wann_lines)):
        if 'num_wann' in wann_lines[flag_wann]:
            nbands = wann_lines[flag_wann].split('=')[1]
            break
    return int(nbands)

def plot_band():
    efermi = get_fermiLevel()
    nbands = get_NBANDS()
    wann = np.loadtxt("wannier90_band.dat", dtype=np.float)
    data_band = np.loadtxt("BAND.dat", skiprows=2, dtype=np.float)
    kline_data = np.loadtxt("KLINES.dat", skiprows=3, dtype=np.float)
    k_num = int(len(wann[:, 0])/(nbands))
    wann[:,1] = wann[:,1]-efermi
    wann_x = wann[0:k_num, 0]

    with open("band_wannier90.dat", "w", encoding='utf-8') as band_wn:
        band_wn.write("Wannier90 band data" + "\n")
        for k in range(0, k_num):
            band_wn.write(str(wann_x[k])[:8].ljust(12, '0') + "   ")
            for h in range(0, nbands):
                band_wn.write(str(wann[h * k_num+k, 1])[:8].ljust(12, '0') + "   ")
            band_wn.write("\n")
    band_wn.close()

    fig, ax = plt.subplots()
    plt.plot(data_band[:, 0], data_band[:, 1], color='red', lw=1, label="PBE")
    plt.plot(kline_data[:, 0], kline_data[:, 1], color='grey', ls='-.', lw=0.5)
    plt.plot(wann_x, wann[0:k_num, 1], color='blue', lw=1, label="Wannier90")
    for j in range(1, nbands):
        plt.plot(wann_x, wann[k_num * j:k_num * j + k_num, 1], color='blue', lw=1)
    plt.legend(loc="upper right")
    plt.xlim(min(data_band[:,0]), max(data_band[:,0]))
    plt.ylabel("E - " + "E$_\mathrm{{F}}$" + " (eV)", fontsize=14, fontname='arial')
    plt.ylim(-3, 3)
    ax.xaxis.set_ticks([])
    plt.savefig("nihe.png", dpi=300)
    plt.show()
if __name__ == '__main__':
    plot_band()