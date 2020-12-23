# -*- coding: utf-8 -*-
"""
Created on 10:21 14-10-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: vasp.py
"""

import os, sys
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)

import scripts.initial

file_path = curPath + "/INCAR_templates"

# main control
print(" ********************************************************************")
print(" *This is a code used to simplify the vasp calculation, written by **")
print(" *****************           XY Ding                *****************")
print(" ********************************************************************")
print(" *************          (^o^)GOOD LUCK!(^o^)           **************")
print("\n")
print(" *********  Pre and post-preparation for vasp calculation  **********")
print(" **********  eg:(1) pre-processing: incar template    ***************")
print(" **********  eg:(2) post-processing: band plot    *******************")
print(" ********************************************************************")
print(" (1) Incar template for vasp calculation")
print(" (2) Data processing for vasp calculation:eg. band and dos")
print(" (3) Elastic constants calculation by using Energy-strain method ")
print(" (4) Wannier90 ")
seltct_vasp = int(input("Input a number: "))
print("\n")
if seltct_vasp == 1:
    print("********  Incar template for vasp calculation !  *******")
    print("(1) optimization" + "\t" + "(2) self-consistent")
    print("(3) band calculation" + "\t" + "(4) dos calculation")
    print("(5) HSE06 calculation" + "\t" + "(6) phonon")
    print("(7) Elastic calculation" + "\t" + "(8) md")
    incar = int(input("Input a Number: "))
    if incar == 1:
        os.system("cp " + file_path +"/INCAR_opt ./INCAR")
    elif incar == 2:
        os.system("cp " + file_path +"/INCAR_scf ./INCAR")
    elif incar == 3:
        os.system("cp " + file_path +"/INCAR_band ./INCAR")
    elif incar == 4:
        os.system("cp " + file_path +"/INCAR_dos ./INCAR")
    elif incar == 5:
        os.system("cp " + file_path +"/INCAR_hse ./INCAR")
    elif incar == 6:
        os.system("cp " + file_path +"/INCAR_phonon ./INCAR")
    elif incar == 7:
        os.system("cp " + file_path +"/INCAR_elastic ./INCAR")
    elif incar == 8:
        os.system("cp " + file_path +"/INCAR_md ./INCAR")
    else:
        print("you have print a wrong number !")
elif seltct_vasp == 2:
    print("********  Data processing for vasp calculation: band and dos ********")
    print("(1) band plot       " + "\t" + "(2) HSE06 band plot ")
    print("(3) Projected band  " + "\t" + "(4) Density of state ")
    print("(5) Optics properties")
    process_num = int(input("Input a Number: "))
    if process_num == 1:
        import scripts.band as bd
        print(" ************************** Enter band ploting ****************************")
        bd.manipulate_bandplot()
    elif process_num == 2:
        import scripts.hse06 as hse
        print(" ************************** Enter HSE06 ploting ****************************")
        hse.manipulate_bandplot()
    elif process_num == 3:
        print(" ****************************** waitting ***********************************")
    elif process_num == 4:
        import scripts.dos as dos
        print(" ************************** Enter DOS ploting ******************************")
        dos._manipulate_dos()
    elif process_num == 5:
        import scripts.optics as op
        print(" ******************** Enter optics properties ploting **********************")
        os.system("cp " + curPath + "/scripts/optics.sh ./")
        os.system("bash optics.sh")
        op.manipulate_plot()
    else:
        print("please input a right number !")
elif seltct_vasp == 3:
    print("*****************************  Mechanical properties  ***************************")
    print("(1) 2D Elastic calculation       ")
    print("(2) 2D Energy-area data      ")
    init = int(input("Input a number: "))
    if init == 1:
        print("********* Using Energy strain method to calculate the elastic constants ******")
        print("(1) prepare the calculating files")
        print("(2) process the results")
        in_p = int(input("Input a number: "))
        if in_p == 1:
            import scripts.deform as deform
            deform.manipulate_deform()
        elif in_p == 2:
            import scripts.elastic_fitting as ef
            ef.output_energy()
            ef.energy_strain_plot()
            ef.energy_strain_fitting_plot()
        else:
            print("waiting !")
    elif init == 2:
        print("********* obtain the es data *********")
        print("(1) prepare the calculating files")
        print("(2) process the results")
        in_s = int(input("Input a number: "))
        if in_s == 1:
            import scripts.es as es
            es.manipulate_deform()
        elif in_s == 2:
            import scripts.es as es
            es.output_es()
            es.plot_es()
        else:
            print("waiting !")
elif seltct_vasp == 4:
    print("********************************  Wannier90  *************************************")
    print("(1) Nihe Band structure      ")
    print("(2) Surface states      ")
    init = int(input("Input a number: "))
    if init == 1:
        import scripts.nihe as nihe
        nihe.plot_band()
    elif init == 2:
        print("waiting !")

else:
    print("You have input a wrong number !")




