# -*- coding: utf-8 -*-
"""
Created on Thu May 21 14:51:22 2020

@author: skidmore
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

SA = pd.read_csv('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production1/solvent-accessable-area-prod1-3.csv',
                        header= None)
SA = pd.read_csv('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production1/solvent-accessable-area-per-res-prod1-3.csv',
                        header= None)
sa = np.array(SA)

plt.figure(dpi=600)
plt.title('Solvent Accessable Surface All Productions', size = 15)
plt.ylabel('Area ($nm^2$)', size = 13)
plt.xlabel('Time (ps)', size = 13)
plt.plot(sa[:,0], sa[:,1], color='blue', label ='Prod1')
plt.plot(sa[:,0], sa[:,2], color='red', label ='Prod2')
plt.plot(sa[:,0], sa[:,3], color='green', label ='Prod3')
#plt.axvline(75, 0, 1, c ='black')
#plt.axvline(505, 0, 1, c='black')
#plt.plot(sa[,:0], sa[,:1], color='red', label ='Prod1')
plt.xlim([-10000, 1025000])
#plt.xlim([-10, 750])
plt.grid(True)
plt.legend(loc='best')





ROG_f = pd.read_csv('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production1/ROG-prod1-Nterm-TM-STAS-XYZ.csv',
                        header= None)
rog = np.array(ROG_f)

ROG_f2 = pd.read_csv('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production2/ROG-prod2-Nterm-TM-STAS-XYZ.csv',
                        header= None)
rog2 = np.array(ROG_f2)

ROG_f3 = pd.read_csv('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production3/ROG-prod3-Nterm-TM-STAS-XYZ.csv',
                        header= None)
rog3 = np.array(ROG_f3)

ROG_f4 = pd.read_csv('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/5da0-control/ROG-5da0-control-TM-STAS-XYZ.csv',
                        header= None)
rog4 = np.array(ROG_f4)

ROG_f5 = pd.read_csv('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/5da0-control2/ROG-5da0-control2-TM-STAS-XYZ.csv',
                        header= None)
rog5 = np.array(ROG_f5)

plt.rc('font', **font)
plt.figure(dpi=600)
plt.title('ROG production 1-3 and controls TM domain Z direction ', size = 15)
plt.ylabel('ROG (nm)', size = 13)
plt.xlabel('Time (ps)', size = 13)
plt.plot(rog[:,0],rog[:,9], color = 'blue', label = 'prod1')
plt.plot(rog2[:,0],rog2[:,9], color = 'red', label = 'prod2')
plt.plot(rog3[:,0],rog3[:,9], color = 'green', label = 'prod3')
plt.plot(rog4[:,0],rog4[:,4], color = 'purple', label = '5da0-control1')
plt.plot(rog5[:,0],rog5[:,4], color = 'orange', label = '5da0-control2')
plt.grid(True)
plt.xlim([-10000, 1025000])
plt.legend(loc='best')







