# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 11:54:58 2020

@author: skidmore
"""

import os
import pandas as pd
#from __future__ import print_function
#%matplotlib inline
import numpy as np
import mdtraj as md

traj = md.load('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production1/every500thstep-prod1.trr', top='C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production1/before_completion/first-part/assembled.gro')

print(traj)
topology = traj.topology

print('How many atoms?    %s' % traj.n_atoms)
print('How many frames?    %s' % traj.n_frames)
print('How many residues? %s' % traj.n_residues)


frame_idx = 4 # zero indexed frame number
atom_idx = 9 # zero indexed atom index
print('Where is the fifth atom at the tenth frame?')
print('x: %s\ty: %s\tz: %s' % tuple(traj.xyz[frame_idx, atom_idx,:]))

print('All residues: %s' % [residue for residue in traj.topology.residues])

#for NA ions
## #of protein atoms X # of NA ion atoms, X # of frames
matrix = np.zeros((1602,415,int(traj.n_frames)))
matrix = np.zeros((1602,70182-69767,int(traj.n_frames)))

# distance for 2 points in 3d space
# d = ((x2 - x1)2 + (y2 - y1)2 + (z2 - z1)2)1/2 
frame_count = 0
for frame in traj.n_frames:
    frame_count += 1
    for atom1_index in range(1602):
        atom1 = np.array(traj.xyz[frame_count, atom1_index,:])
        # atom2 is ion/lipid, range from first to last of thet type. 
        for atom2_index in range(69767:70182):
            atom2 = np.array(traj.xyz[frame_count,atom2_index,:])
            distance = 
            matrix[,,frame_count] = distance

frame_count = 0
for frame in range(traj.n_frames):
    for atom1_index in range(1091,1602):
        atom1 = np.array(traj.xyz[frame, atom1_index,:])
        # atom2 is ion/lipid, range from first to last of thet type. 
        for atom2_index in range(330,415):
            atom2 = np.array(traj.xyz[frame,atom2_index,:])
            distance = (((atom2[0:1]-atom1[0:1])**2) + ((atom2[1:2]-atom1[1:2])**2)+((atom2[2:3]-atom1[2:3])**2))**.5
            matrix[atom1_index,atom2_index,frame] = distance
            if distance <= 5:
                i_list.append((atom1_index,atom2_index,frame, float(distance)))

frame_count = 0
i_list = []

for frame in range(traj.n_frames):
    for atom1_index in range(1602):
        atom1 = np.array(traj.xyz[frame, atom1_index,:])
        # atom2 is ion/lipid, range from first to last of thet type. 
        ff =0
        for atom2_index in range(69767,70182):
            
            atom2 = np.array(traj.xyz[frame,atom2_index,:])
            distance = (((atom2[0:1]-atom1[0:1])**2) + ((atom2[1:2]-atom1[1:2])**2)+((atom2[2:3]-atom1[2:3])**2))**.5
            matrix[atom1_index,ff,frame] = distance
            ff+=1
            if distance <= 5:
                i_list.append((atom1_index,atom2_index,frame, int(distance)))
            
            
            
print(matrix.shape)
print(matrix[109,414])
print(len(matrix[1601,414]))

##Find framewise interactions
for x in range(2002):
    frame_dict[x] = []
for n in range(1602):
    for m in range(415):
        frame=0
        for x in matrix[n,m]:
            #print(x)
            frame +=1
            #print(frame)
            if x <=5:
                frame_dict[frame].append([n,m])
                
for n in range(1602):
    for x in matrix[n,]:
        
                
                
                
    print(x[0,1])
    #print(len(x))
    for y in x:
        print(y)
        print(len(y))



count =0
for x in i_list:
    #select frame to look at. Right side of ==.
    if x[2] == 0:
        count +=1
        print(x)

i_dict = {}
for i in range(1091,1602):
    atom_ls = []
    for x in i_list:
        if x[0] == i:
            atom_ls.append((x))
        i_dict[i] = atom_ls

atom1 = np.array(traj.xyz[1, 1,:])
atom2 = np.array(traj.xyz[1, 2,:])
print(atom1[1:2])
print((atom1[0:1]) - (a2[0:1]))
##want to make a 3D matrix of columns = residues/protein atoms, rows = ions atoms/ lipid atoms, z = distances for column vs row atoms by frame











































resion = pd.read_csv("C:/Users/skidmore/Desktop/Prestin-MD/skidmore/1-24-2020/interactions/residue-to-ion-distances-500step.csv")
ionres = pd.read_csv("C:/Users/skidmore/Desktop/Prestin-MD/skidmore/1-24-2020/interactions/Ion-to-residue-distances-500step.csv")


resion.head()
print(resion['MET0'])
resiontime = {}
resdata ={}
for x in resion: 
    #print(x)
    y_count = 0
    rit = []
    timelst = []
    for y in resion[x]:
        #print(x, y, y_count)
        if float(y) <= 0.6:
            #print(y, y_count, resion.iloc[y_count]["Time"])
            #resiontime[x] = (y, resion.iloc[y_count]["Time"])
            
            if float(resion.iloc[y_count]["Time"]) in timelst:
                rit.append((float(y),float(resion.iloc[y_count]["Time"]+500)))
            else:
                rit.append((float(y),float(resion.iloc[y_count]["Time"])))
                timelst.append(float(resion.iloc[y_count]["Time"]))
        y_count += 1
    if rit != []:
        resdata[x] = rit
        ## THINGS IN resdata with large size, have long timescale interactions. 
del resdata['Time']

merged_dict ={}
md2 ={}
for j in ionres:
    k_count = 0
    irt = []
    timelst = []
    for k in ionres[j]:
        ion_dist_time = []
        if float(k) <= 0.6 :
            dtl = (float(k),float(ionres.iloc[k_count]["Time"]))
            #if dtl in resiontime.values():
            #    merged_dict[]
            #ion_dist_time.append((j,k,ionres.iloc[k_count]["Time"]))
            for ke, val in resdata.items():
                #print(ke, val, dtl)
                if dtl in val:
                    #print('YES')
                    irt.append((ke, j, dtl))
                   
                    #merged_dict[ke] = ion_dist_time
        k_count +=1
    ijlst = set()
    for i, j, k in irt:
        if (j+i) in ijlst:
            md2[str(j+i)].append(k)
        else:
            ijlst.add(j+i)
            #md2 shows specific residue, ion interactions. the size of the values reveal how long these interactions last.
            #NOte: other ions could be interacting at different times with this residue. SEE resdata
            md2[str(j+i)] = [k]
    ## merged_dict shows us which ions have the most interactions (see size of values)
    merged_dict[j] = irt
    
resid = {}    
for k, v in resdata.items():
    resid[int(k[3:])] = len(v)
    
    
import matplotlib.pyplot as plt

plt.figure(dpi=600)
plt.scatter(resid.keys(), resid.values())
plt.title('Interactions <= 0.6nm between Residues and Ions')
plt.xlabel("Residue")
plt.ylabel("Number of interactions")
plt.show()    
    
    
    
    
    
    
    
    
    