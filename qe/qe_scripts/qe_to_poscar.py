# -*- coding: utf-8 -*-
"""
Created on 2:06 AM 28-11-2020 

@author: X-Y Ding
mail to: dxy_vasp@163.com
python3: qe_to_poscar.py
"""
import numpy as np
from matplotlib import pyplot as plt
import math
import os, sys
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)
import get_data as gd

# transform relaxation of Quantum Espresso to POSCAR
def manipulate():
    print("*********************** QE output to POSCAR Formate ******************************")
    print("(1) Use the default setting: relax.out !")
    print("(2) User re-define the output filename for relaxation !")
    init = int(input("input a number: "))
    if init == 1:
        gd.trans_qe_pos()
    elif init == 2:
        filename = str(input("Input the output filename for relaxation: "))
        gd.trans_qe_pos(filename)
    else:
        print("Wrong input numbers !")

if __name__ == '__main__':
    manipulate()