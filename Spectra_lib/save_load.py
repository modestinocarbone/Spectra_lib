"""
 *  File:        save_load.py
 *  Author:      Modestino Carbone
 *  Created on:  02/10/2025
 *  Last update: 02/10/2025
 *
 *  Description:
 *      Here are collected only graphical functions.
 *
"""

import numpy as np
import os

# save file

def save(f, spec, file_name):
    
    if os.path.exists(file_name):
        print(f"the '{file_name}' already exixsts and it will be replaced")

    matrice = np.column_stack((f, spec))
    np.savetxt(file_name, matrice, delimiter='\t')

def savewstyle(f, spec, avg, dec, file_name):

    if os.path.exists(file_name):
        print(f"the '{file_name}' already exixsts and it will be replaced")

    matrice = np.column_stack((f, spec, avg, dec))
    np.savetxt(file_name, matrice, delimiter='\t')

# load file

def load(file_name):
    
    # check if present
    if os.path.exists(file_name):

        data = np.loadtxt(file_name, delimiter='\t')
        f = data[:, 0]
        spec = data[:, 1]

    else:
        raise Exception(f"File '{file_name}' not found")
    
    return f, spec

def loadwstyle(file_name):

    # check if present
    if os.path.exists(file_name):

        data = np.loadtxt(file_name, delimiter='\t')
        f    = data[:, 0]
        spec = data[:, 1]
        avg    = data[:, 0]
        dec = data[:, 1]

    else:
        raise Exception(f"File '{file_name}' not found")
    
    return f, spec, avg, dec


