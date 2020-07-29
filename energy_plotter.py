# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 15:45:02 2020

@author: skidmore
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import savgol_filter

distances1 = pd.read_csv('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production1/production1_PRMSD.csv',
                        header= None)
distances1 = np.array(distances1)

distances2 = pd.read_csv('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production2/production2_PRMSD.csv',
                        header= None)
distances2 = np.array(distances2)

distances3 = pd.read_csv('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production3/production3_PRMSD.csv',
                        header= None)
distances3 = np.array(distances3)

distances = pd.read_csv('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/5da0-control/5da0_allatoms_PRMSD.csv',
                        header= None)
distances1 = np.array(distances)

tab = pd.read_csv('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/potential-energy.csv')
ys = []
ys.append(list(tab['Energy1']))
ys.append(list(tab['Energy2']))
ys.append(list(tab['Energy3']))

fig, ax1 = plt.subplots(dpi=600)
color = 'tab:red'
plt.title('N-term-trial Potential Energy')
ax1.set_ylabel('Potential Energy J/mol', color = color)
ax1.set_xlabel('Time (ps)')
ax1.plot(tab['Time'],tab['Energy3'], color = color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel('PRMSD ref to ref (nm) BB', color=color)
ax2.plot(df,di, color=color)
ax2.tick_params(axis='y', labelcolor=color)
fig.tight_layout()

#plt.plot(tab['Time'],tab['Energy3'])
#plt.plot(tab['Time'],tab['Energy1'])    

plt.plot(df,di)


df = []
#dd = []
d10 = []
d20 = []
#dc = []
d1l = []
#d2l = []
d10.append(0)
#d20.append(0)
for i in range(2001):
    df.append(i*.5)
    #dd.append(np.mean(distances[i]))
    #print(distances[1,4])
    d10.append(abs((distances[0,i+1])))
    #dc.append(abs((distances[972,i+1])))
    d1l.append(abs((distances[2000,i+1])))
    #d20.append(abs((distances2[0,i+1])))
    #d2c.append(abs((distances[972,i+1])))
    #d2l.append(abs((distances2[2000,i+1])))
#dc.insert(972, 0)
d1l.append(0)
#d2l.append(0)


d3f = []
#dd = []
d30 = []
#dc = []
d3l = []
d30.append(0)
for i in range(2028):
    d3f.append(i*.5)
    #dd.append(np.mean(distances[i]))
    #print(distances[1,4])
    d30.append(abs((distances3[0,i+1])))
    #dc.append(abs((distances[972,i+1])))
    d3l.append(abs((distances3[2027,i+1])))
#dc.insert(972, 0)
d3l.append(0)
#df.pop(-1)
#plt.scatter(tab['Time'], ys)




## compare energy to PRMSD ref by ref
fig, ax1 = plt.subplots(dpi=600)
color = 'tab:red'
plt.title('Production run 2 Potential Energy and PRMSD ')
ax1.set_ylabel('Potential Energy J/mol', color = color)
ax1.set_xlabel('Time (ps)')
ax1.plot(tab['Time'],tab['Energy1'], color = color)
ax1.tick_params(axis='y', labelcolor=color)
#yhat = savgol_filter(tab['Energy2'], 11, 3)
#ax1.plot(tab['Time'],yhat, color = 'tab:orange')

ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel('PRMSD ref to first frame (nm) ', color=color)
ax2.plot(df,d0, color=color)
ax2.tick_params(axis='y', labelcolor=color)
fig.tight_layout()



##framewise comparision, overlay
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}

plt.rc('font', **font)
plt.figure(dpi=600)
plt.title('N-term-trial STAS PRMSD first frame comparisons', size = 15)
plt.ylabel('PRMSD first frame to frame comparison', size = 13)
plt.xlabel('Time (ns)', size = 13)
plt.xticks([0,100,200,300,400,500,600,700,800,900,1000], rotation ='vertical')
plt.plot(df,d10, color = 'blue', label = 'Run 1')
#plt.plot(df,d20, color = 'red', label = 'Run 2')
#plt.plot(d3f,d30, color = 'green', label = 'Run 3')
plt.axvline(250, 0, 1, c ='black')
plt.axvline(300, 0, 1, c='black')
plt.axvline(650, 0, 1, c='purple')
plt.axvline(700, 0, 1, c='purple')
plt.axvline(800, 0, 1, c='orange')
plt.axvline(850, 0, 1, c='orange')
plt.xlim([-10, 1020])
plt.ylim([0, 0.7])
plt.grid(True)
#plt.legend(loc='best')

plt.figure(dpi=600)
plt.title('N-term-trial STAS PRMSD last frame comparisons', size = 15)
plt.ylabel('PRMSD last frame to frame comparison', size = 13)
plt.xlabel('Time (ns)', size = 13)
plt.xticks([0,100,200,300,400,500,600,700,800,900,1000], rotation ='vertical')
plt.plot(df,d1l, color = 'blue', label = 'Run 1')
#plt.plot(df,d2l, color = 'red', label = 'Run 2')
#plt.plot(d3f,d3l, color = 'green', label = 'Run 3')
plt.axvline(225, 0, 1, c='black')
plt.axvline(275, 0, 1, c='black')
plt.axvline(550, 0, 1, c ='purple')
plt.axvline(600, 0, 1, c='purple')
plt.axvline(800, .0, 1, c='orange')
plt.axvline(850, 0, 1, c='orange')
#plt.axvline(675, .4, .65, c='orange')
#plt.axvline(775, .4, .65, c='orange')
plt.xlim([-10, 1020])
plt.ylim([0, 0.7])
plt.grid(True)
#plt.legend(loc='best')

#second derivative to get centroid of flat region
fc = (np.diff(dc, 2)).tolist()
fc.append(0)
fc.append(0)

f10 = (np.diff(d10,1)).tolist()
f10.append(0)
f10.append(0)

fl = (np.diff(dl,2)).tolist()
fl.append(0)
fl.append(0)

plt.figure(dpi=600)
plt.title('Production run 3 PRMSD 2nd derivative frame comparisons')
plt.xlabel('Time (ps)')
plt.xticks([0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000], rotation ='vertical')
plt.grid(True)
plt.plot(df,f10, color = 'tab:blue', label = 'first frame')
#plt.plot(df,fl, color = 'tab:purple', label = 'last frame')
#plt.plot(df,fc, color = 'tab:green', label = 'centroid frame')
plt.legend(loc='best')

###now view as summed elementwise abs values and not
fff = [sum(x) for x in zip(f0, fl,fc)]
f0a = [abs(x) for x in f0]
fca = [abs(x) for x in fc]
fla = [abs(x) for x in fl]
fffa = [sum(x) for x in zip(f0a, fla,fca)]


plt.figure(dpi=600)
plt.title('Production run 3 PRMSD 2nd derivative summed')
plt.xlabel('Time (ps)')
plt.xticks([0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000], rotation ='vertical')
plt.grid(True)
plt.plot(df,fff)

plt.figure(dpi=600)
plt.title('Production run 3 PRMSD 2nd derivative summed absolute')
plt.xlabel('Time (ps)')
plt.xticks([0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000], rotation ='vertical')
plt.grid(True)
plt.plot(df,fffa)

##find region where we want the representative frame to come from then get centroid
lol=[]
lol2=[]
for x in range(0, 2000, 100):
    lol.append((x,x+100,np.mean(fffa[x:x+100])))
    lol2.append(np.mean(fffa[x:x+100]))
np.min(lol2)

beta = 1
m1 = 1500
m2 = 1600
index = np.exp(-beta*distances1[m1:m2] / distances1[m1:m2].std()).sum(axis=1).argmax()
print(index+m1)
index = np.exp(-beta*distances2[m1:m2] / distances2[m1:m2].std()).sum(axis=1).argmax()
print(index+m1)
index = np.exp(-beta*distances3[m1:m2] / distances3[m1:m2].std()).sum(axis=1).argmax()
print(index+m1)


ROG_f = pd.read_csv('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production1/ROG_prod_all.csv',
                        header= None)
rog = np.array(ROG_f)

ROG_f = pd.read_csv('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/5da0-control2/ROG-5da0-control2.csv',
                        header= None)
rog = np.array(ROG_f)

plt.rc('font', **font)
plt.figure(dpi=600)
plt.title('ROG 5da0 control 2 X, Y, Z', size = 15)
plt.ylabel('ROG (nm)', size = 13)
plt.xlabel('Time (ns)', size = 13)
plt.plot(df,rog[:,1], color = 'blue', label = 'X ROG')
plt.plot(df,rog[:,2], color = 'red', label = 'Y ROG')
plt.plot(df,rog[:,3], color = 'green', label = 'Z ROG')
plt.grid(True)
plt.xlim([-10, 1010])
plt.legend(loc='best')


'''
fig, ax1 = plt.subplots(dpi=600)
color = 'tab:red'
plt.title('5da0 control 2 ROG vs PRMSD to frame 1')
ax1.set_ylabel('ROG (nm)', color = color)
ax1.set_xlabel('Time (ns)')
ax1.plot(df,rog[:,1], color = 'red')
plt.grid(True)
#ax1.plot(d3f,rog[:,2], color = 'red')
#ax1.plot(d3f,rog[:,3], color = 'green')
ax1.tick_params(axis='y', labelcolor='red')

ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel('PRMSD to frame 1', color=color)
ax2.plot(df, d1l,color='blue')
#ax2.plot(df,d20, color='red')
#ax2.plot(d3f,d30, color='green')
ax2.tick_params(axis='y', labelcolor='blue')
#plt.plot(df,d2l, color = 'red', label = 'Run 2')
plt.xlim([-10, 1020])
plt.ylim([0, 0.7])
#plt.grid(True)
fig.tight_layout()'''






fig, ax1 = plt.subplots(dpi=600)
color = 'tab:red'
plt.title('5da0 control 2 ROG and PRMSD to first frame')
ax1.set_ylabel('ROG (nm)', color = color)
ax1.set_xlabel('Time (ns)')
ax1.plot(df,rog[:,1], color = 'red')
plt.axvline(475, 0, 1, c ='black')
plt.axvline(525, 0, 1, c='black')
plt.axvline(600, 0, 1, c='purple')
plt.axvline(650, 0, 1, c='purple')
plt.axvline(950, 0, 1, c='orange')
plt.axvline(1000, 0, 1, c='orange')
plt.grid(True)
ax1.tick_params(axis='y', labelcolor='red')

ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel('PRMSD to last frame', color=color)
ax2.plot(df, d10,color='blue')
ax2.tick_params(axis='y', labelcolor='blue')
plt.xlim([-10, 1020])
plt.ylim([0, 0.7])
fig.tight_layout()


print(distances1[760,1565])





ROG_3 = pd.read_csv('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production3/ROG_prod_3-sum-x-y-z.csv',
                        header= None)
rog3 = np.array(ROG_3)

ROG_2 = pd.read_csv('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production2/ROG_prod2-x-y-z.csv',
                        header= None)
rog2 = np.array(ROG_2)

ROG_1 = pd.read_csv('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production1/ROG_prod1-x-y-z.csv',
                        header= None)
rog1 = np.array(ROG_1)
                
plt.figure(dpi=600)
plt.title('ROG Zs All Productions Comparison', size = 15)
plt.ylabel('ROG (nm)', size = 13)
plt.xlabel('Time (ns)', size = 13)
#plt.xticks([0,100,200,300,400,500,600,700,800,900,1000], rotation ='vertical')
plt.plot(df,rog1[:,4], color = 'blue', label = 'Prod1')
plt.plot(df,rog2[:,4], color = 'red', label = 'Prod2')
plt.plot(d3f,rog3[:,4], color = 'green', label = 'Prod3')
#plt.plot(df,rog1[:,4], color = 'purple', label = 'Z')
plt.xlim([-10, 1020])
plt.grid(True)
plt.legend(loc='best')





basic = ['ARG','HIS','LYS']
acidic = ['ASP', 'GLU']
polar = ['SER', 'THR', 'ASN', 'GLN', 'GLY']
nonpolar = ['PRO', 'ALA','LEU', 'MET', 'PHE', 'TRP', 'TYR','VAL', 'ILE', 'CYS']
#plt.plot(tab['Time'],tab['Energy3'])
#plt.plot(tab['Time'],tab['Energy1'])    

#plt.plot(df,di)

fff= pd.read_csv('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/n_term-trial/n-term-prod3-compare.csv', header = None)
polarlst = []
basiccount = 0
acidiccount =0
polarcount = 0
nonpolarcount = 0
for x in fff[0]:
    numb = x[:3]
    y = x[3:]
    if y in basic:
        polarlst.append((numb, y, 'basic'))
        basiccount +=1
        polarcount +=1
    elif y in acidic:
        polarlst.append((numb, y, 'acidic'))
        acidiccount +=1
        polarcount +=1
    elif y in polar:
        polarlst.append((numb, y, 'polar'))
        polarcount +=1
    elif y in nonpolar:
        polarlst.append((numb, y, 'nonpolar'))
        nonpolarcount +=1








