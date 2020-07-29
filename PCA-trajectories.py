# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 12:08:49 2020

@author: skidmore
"""

import numpy as np
import mdtraj as md
from sklearn import decomposition
traj = md.load('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production2/prod2-every500th.trr', top='C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production2/skid-production2.gro')

print(traj)

#do PCA w/ scikit learn
xyz = traj.xyz.reshape((traj.xyz.shape[0], -1))
pca = decomposition.PCA()
pca.fit(xyz)
xyz_transformed = pca.transform(xyz)
print(xyz_transformed)

#visualize
import matplotlib.pyplot as plt
plt.figure(dpi=600)
plt.scatter(xyz_transformed[:, 0], xyz_transformed[:, 1], c=range(len(xyz_transformed)))
plt.title('Primary component analysis 500ps frames')
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.show()

#Visualize with clusters highlighted and labled?
plt.figure(dpi=600)
plt.scatter(xyz_transformed[:, 0], xyz_transformed[:, 1], s=1)
plt.scatter(xyz_transformed[327, 0],xyz_transformed[327, 1],  c='r', s=10, )
plt.text(xyz_transformed[327, 0]+3,xyz_transformed[327, 1]+3, 'Clust1', fontsize=9, color ='r')
plt.scatter(xyz_transformed[699, 0],xyz_transformed[699, 1],  c='orange', s=10)
plt.text(xyz_transformed[699, 0]+3,xyz_transformed[699, 1]+3, 'Clust2', fontsize=9, color ='orange')
plt.scatter(xyz_transformed[1634, 0],xyz_transformed[1634, 1],  c='y', s=10)
plt.text(xyz_transformed[1634, 0]+3,xyz_transformed[1634, 1]+3, 'Clust3', fontsize=9, color ='y')
plt.scatter(xyz_transformed[84, 0],xyz_transformed[84, 1],  c='lime', s=10)
plt.text(xyz_transformed[84, 0]+3,xyz_transformed[84, 1]+3, 'Clust4', fontsize=9, color ='lime')
plt.scatter(xyz_transformed[229, 0],xyz_transformed[229, 1],  c='pink', s=10)
plt.text(xyz_transformed[229, 0]+3,xyz_transformed[229, 1]+3, 'Clust5', fontsize=9, color ='pink')
plt.title('Primary component analysis 500ps frames')
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.show()
