# -*- coding: utf-8 -*-
"""
Created on 10:21 14-10-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: main.py
"""

# import packages
import os, sys
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)

file_path = curPath + "/INCAR_templates"

# main control
print(" ********************************************************************")
print(" *This is a code used to simplify the vasp calculation, written by **")
print(" *****************           XY Ding                *****************")
print(" ********************************************************************")
print(" *************          (^o^)GOOD LUCK!(^o^)           **************")
print("\n")
print(" *********  Pre and post-preparation for VASP&QE calculation  **********")
print(" **********  eg:(1) pre-processing: incar template    ***************")
print(" **********  eg:(2) post-processing: band plot    *******************")
print(" ********************************************************************")
print(" (1) VASP")
print(" (2) Quantum Espresso")
seltct = int(input("Input a number: "))
print("\n")
if seltct == 1:
    print("********  Entering VASP calculation !  *******")
    import vasp.vasp

elif seltct == 2:
    print("********  Entering VASP calculation !  *******")
    import qe.qe
else:
    print("You are input a wrong number !")







