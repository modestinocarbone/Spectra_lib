"""
 *  File:        graph.py
 *  Author:      Modestino Carbone
 *  Created on:  21/03/2025
 *  Last update: 01/10/2025
 *
 *  Description:
 *      Here are collected only graphical functions.
 *
"""

import matplotlib.pyplot as plt

def marker_plotter_MA(averages, decades, label):

    # Every decade has a fixed shade.
    # Number of maximum shades can be overcome by user.
    # A flexible implementation should be implemented.
    colors = [
        'lightblue', 'lightgreen', 'lightcoral',  'lightpink',
        'lavender', 'lightcyan', 'peachpuff', 'mintcream',
        'honeydew', 'thistle', 'powderblue', 'moccasin',
        'palegoldenrod', 'lemonchiffon', 'lightgoldenrodyellow',
        'mistyrose'
    ]
    
    # Draw loop
    for i  in range(len(decades)-1):
        
        # Add dashed lines
        plt.axvline(x=decades[i-1], color='gray', linestyle='--', lw=1)
        plt.axvline(x=decades[i], color='gray', linestyle='--', lw=1)
        
        # Color decade and insert legend
        if(label):
            if(averages[i]==1):
                plt.axvspan(decades[i], decades[i+1], color=colors[i], alpha=0.5, label=str(int(averages[i]))+' average')
            else:
                plt.axvspan(decades[i], decades[i+1], color=colors[i], alpha=0.5, label=str(int(averages[i]))+' averages')
        else:
                plt.axvspan(decades[i], decades[i+1], color=colors[i], alpha=0.5)
