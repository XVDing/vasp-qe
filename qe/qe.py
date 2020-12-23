# -*- coding: utf-8 -*-
"""
Created on 10:21 14-10-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: qe.py
"""

import os, sys
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)

import qe_scripts.band as bd
import qe_scripts.get_data
import  qe_scripts.qe_to_poscar as qp
import input_qe.qe_input as initial


file_path = curPath + "/input_qe"

# main control
print(" ********************************************************************")
print(" *This is a code used to simplify the QE calculation, written by **")
print(" *****************           XY Ding                *****************")
print(" ********************************************************************")
print(" *************          (^o^)GOOD LUCK!(^o^)           **************")
print("\n")
print(" *********  Pre and post-preparation for QE calculation  ***********")
print(" **********  eg:(1) pre-processing: input template    ***************")
print(" **********  eg:(2) post-processing: band plot    *******************")
print(" ********************************************************************")
print(" (1) Input template for Quantum Espresso calculation")
print(" (2) Data processing for Quantum Espresso: eg. qe to POSCAR !")
print(" (3) Data processing for Quantum Espresso  calculation:eg. band and dos")
seltct_qe = int(input("Input a number: "))
print("\n")
if seltct_qe == 1:
    print("********  Incar template for QE calculation !  *******")
    initial.manipulate_qe()
elif seltct_qe == 2:
    print("Notation: You must have a input file of relaxation called")
    print("********************** relax ********************* ")
    qp.manipulate()
elif seltct_qe == 3:
    print("********  Data processing for QE calculation: band and dos ********")
    print("**** Notation: You must have a input file of relaxation called ****")
    print("****************************** relax ****************************** ")
    print("(1) band plot       " + "\t" + "(2) HSE06 band plot ")
    print("(3) Density of state" + "\t" + "(4) Projected band ")
    print("(5) Optics properties")
    process_num = int(input("Input a Number: "))
    if process_num == 1:
        print(" ************************* Enter band ploting ****************************")
        bd.manipulate_bandplot()
    elif process_num == 2:
        print(" ************************* Enter HSE06 ploting ****************************")
        print("Waiting")
    elif process_num == 3:
        print("************************** Enter DOS ploting ******************************")
    elif process_num == 4:
        print("*********************** Enter projected band ploting ******************************")
    elif process_num == 5:
        print(" ******************** Enter optics properties ploting **********************")
    else:
        print("please input a right number !")
else:
    print("You have input a wrong number !")

