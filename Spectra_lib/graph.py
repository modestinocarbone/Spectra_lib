# -*- coding: utf-8 -*-
"""
Created on Sat Mar 22 08:20:55 2025

@author: modes
"""

import matplotlib.pyplot as plt


def marker_plotter_MA(avgerages, decades, label):

    colors = [
    'lightblue', 'lightgreen', 'lightcoral',  'lightpink', 'lavender', 'lightcyan', 'peachpuff', 'mintcream', 'honeydew', 'thistle',
    'powderblue', 'moccasin', 'palegoldenrod', 'lemonchiffon', 'lightgoldenrodyellow', 'mistyrose'
    ]
    
    
    # Loop per tracciare le linee tratteggiate e colorare le porzioni
    for i  in range(len(decades)-1):
        
        # Aggiungere linee tratteggiate
        plt.axvline(x=decades[i-1], color='gray', linestyle='--', lw=1)
        plt.axvline(x=decades[i], color='gray', linestyle='--', lw=1)
        # Colorare la porzione
        
        if(label):
            if(avgerages[i]==1):
                plt.axvspan(decades[i],decades[i+1], color=colors[i], alpha=0.5, label=str(int(avgerages[i]))+' average')
            else:
                plt.axvspan(decades[i],decades[i+1], color=colors[i], alpha=0.5, label=str(int(avgerages[i]))+' averages')
        else:
                plt.axvspan(decades[i],decades[i+1], color=colors[i], alpha=0.5)
