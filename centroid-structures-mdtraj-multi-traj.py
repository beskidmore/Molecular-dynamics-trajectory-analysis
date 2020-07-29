# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 09:25:21 2020

@author: skidmore
"""
## Centroid identification using MDtraj. This only gets the centroid structure of the frames provided. It doesn't do the clustering
from __future__ import print_function
#%matplotlib inline
import numpy as np
import mdtraj as md
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy
from scipy.cluster.vq import vq, kmeans, whiten
from scipy.spatial.distance import squareform
from sklearn.neighbors import KernelDensity
from sklearn import metrics
from scipy import stats
import gc
import pandas as pd
from joblib import Parallel, delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()


def fancy_dendrogram(*args, **kwargs):
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = scipy.cluster.hierarchy.dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('RMSD Average linkage hierarchical clustering')
        plt.xlabel('Sample index or (Cluster size)')
        plt.ylabel('Distance Score')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata


traj = md.load('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production1/every500thstep-prod1.trr', top='C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production1/before_completion/first-part/assembled.gro')
traj2 = md.load('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production2/500ps-frames/prod2-every500th.trr', top='C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production2/skid-production2.gro')
traj3 = md.load('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production3/prod3-every500.trr', top='C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production3/skid-production3.gro')
n_termtest = md.load('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/n_term-trial/n-term-production-500step.trr', top='C:/Users/skidmore/Desktop/Prestin-MD/skidmore/n_term-trial/n-term-production-last-frame.gro')
control = md.load('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/5da0-control/5da0-prod-500step.trr', top='C:/Users/skidmore/Desktop/Prestin-MD/skidmore/5da0-control/5da0-production.gro')
control2 = md.load('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/5da0-control2/5da0-control-2-500step.xtc', top='C:/Users/skidmore/Desktop/Prestin-MD/skidmore/5da0-control2/5da0-production2.gro')
control_merge = md.load('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/5da0-merged/control1-2-merge-500step.trr', top='C:/Users/skidmore/Desktop/Prestin-MD/skidmore/5da0-control2/5da0-production2.gro')
traj_merge = md.load('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/merged-traj-analysis/merged-traj-prod1and2and3.trr', top='C:/Users/skidmore/Desktop/Prestin-MD/skidmore/production2/skid-production2.gro')
assembly = md.load('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/assembly/prod1-3assembly-500step.trr', top='C:/Users/skidmore/Desktop/Prestin-MD/skidmore/assembly/assembled.gro')
tm_only = md.load('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/tm-domain-only/tm_prod_500step_center.xtc', top='C:/Users/skidmore/Desktop/Prestin-MD/skidmore/tm-domain-only/tm-self-assembly.gro')
print(traj)
traj = tm_only
topology = traj.topology

#control 5da0
atom_indices = topology.select('resid 0 to 466')
#control STAS
atom_indices = topology.select('resid 378 to 466')
#control TM
atom_indices = topology.select('resid 0 to 377')
atom_indices = topology.select('resid 0 to 430')
»
atom_indices = topology.select('name BB')
atom_indices = topology.select('resid 0 to 743')
#TM domain
atom_indices = topology.select('resid 74 to 504')
#STAS
atom_indices = topology.select('resid 505 to 743')

distances = np.empty((traj.n_frames, traj.n_frames)) 
for i in range(traj.n_frames):
    distances[i] = md.rmsd(traj, traj, i, atom_indices=atom_indices)
    
print('Max pairwise rmsd: %f nm' % np.max(distances))
print('Min pairwise rmsd: %f nm' % np.min(distances))
print('Mean pairwise rmsd: %f nm' % np.mean(distances))
print(np.where(distances == np.max(distances)))

##To compare between two different trajectories need to use gromacs trjcat



#md.r
#pool = multiprocessing.Pool(3)
#num_cores = multiprocessing.cpu_count()

## try to run in parallel
#def rmsder(traj_merge, traj_merge2, ranger, atom_indices):
#    for i in ranger: 
#        distances_merge[i] = md.rmsd(traj_merge, traj_merge2, i, atom_indices)

#results = Parallel(n_jobs=5)(delayed(rmsder)((traj_merge,traj_merge,range(traj_merge.n_frames),atom_indices_merge)))



#plot distribuition, need to NaN lower half of symmetric matrix

corr_dist = np.triu(distances, k=1)
corr_dist[np.tril_indices(corr_dist.shape[0], -1)] = np.nan
a = corr_dist[np.logical_not(np.isnan(corr_dist))]
#a= a.tolist()
a = np.array(a).reshape(-1,1)
bins = [0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,
                               0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,
                               0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1]
bins = [0.85,0.875,0.9,0.925,0.95,0.975,1]


plt.figure(dpi=600)
plt.title('PRMSD distribuiton TM only Production')
plt.ylabel('Count')
#plt.ylabel('Probability Density')
plt.xlabel('PRMSD')
#sns.distplot(corr_dist,bins = [0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,
#                               0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85],
#             kde = True, hist_kws={'edgecolor':'black'},kde_kws={'linewidth': 3})
sns.distplot(corr_dist,bins = [0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,
                               0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85],
             kde = True, hist = False)
sns.distplot(corr_dist,bins = [0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,
                               0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85],
             kde = False, hist = True, norm_hist=False, hist_kws={'edgecolor':'black'},kde_kws={'linewidth': 3})
#plt.plot([0.375, 0.375], [0,5], color = 'red')
#plt.plot([0.45, 0.45], [0,5], color = 'red')
#weights = np.ones_like(distances_merge)/float(len(distances_merge))
#plt.hist(corr_dist, , bins=50)
plt.show()

#count of bins, right =False matches behaviour of sns.distplot() and plt.hist()
(n,bins, patches) =plt.hist(a, bins = bins)
plt.show()
In [84]: print(n)


#test for nomalcy
k2, p = stats.normaltest(a)
alpha = 1e-3
print(p)
if p < alpha:  # null hypothesis: x comes from a normal distribution
    print("The null hypothesis can be rejected")
else:
    print("The null hypothesis cannot be rejected")

# try to get AUC aka. probability for bins

def plot_prob_density(a):
    plt.figure(figsize = (10, 7))

    unit = 1.5
    #x = np.linspace(a.min(), a.max(), 1000)[:, np.newaxis]
    x = np.linspace(0.0, 1.0, 1000)[:, np.newaxis]

    # Plot the data using a normalized histogram
    plt.hist(a, bins=50, density=True )
    #plt.hist(df_dinner, bins=10, density=True, label='Dinner Time', color='navy', alpha=0.2)

    # Do kernel density estimation
    kd_dist = KernelDensity(breadth_first=True, kernel='gaussian', bandwidth=0.005).fit(a)
    #kd_dinner = KernelDensity(kernel='gaussian', bandwidth=0.5).fit(df_dinner)

    # Plot the estimated densty
    kd_dist_vals = np.exp(kd_dist.score_samples(x))
    #kd_vals_dinner = np.exp(kd_dinner.score_samples(x))

    plt.plot(x, kd_dist_vals, color='orange')
    #plt.plot(x, kd_vals_dinner, color='navy')
    
    #plt.axvline(x=x_start,color='red',linestyle='dashed')
    #plt.axvline(x=x_end,color='red',linestyle='dashed')

    # Show the plots
    plt.xlabel('PRMSD', fontsize=15)
    plt.ylabel('Probability Density', fontsize=15)
    plt.legend(fontsize=15)
    plt.show()
    gc.collect()
    return kd_dist

def get_probability(start_value, end_value, eval_points, kd):
    
    # Number of evaluation points 
    N = eval_points                                      
    step = (end_value - start_value) / (N - 1)  # Step size

    x = np.linspace(start_value, end_value, N)[:, np.newaxis]  # Generate values in the range
    kd_vals = np.exp(kd.score_samples(x))  # Get PDF values for each x
    probability = np.sum(kd_vals * step)  # Approximate the integral of the PDF
    return probability.round(8)

#a = corr_dist.tolist()
a = corr_dist[np.logical_not(np.isnan(corr_dist))]
#a= a.tolist()
a = np.array(a).reshape(-1,1)

unit = 1.5
x = np.linspace(0.85, 0.875, 1000)[:, np.newaxis]
kd_dist = KernelDensity(kernel='gaussian', bandwidth=0.005).fit(a)
plot_prob_density(a)
#kd_dist_vals = np.exp(kd_dist.score_samples(x))
#plt.plot(x, kd_dist_vals, color='orange')


'''counter =0 
for x in a:
    if x >= 0.4 and x <=0.449999999999999:
        counter +=1'''
##Get the probability between two verticals. i.e. the AUC for given range        
'''print('Probability of PRMSD 0.4-0.45: {}'
      .format(get_probability(start_value = 0.0, 
              end_value = 0.025, 
              eval_points = 1000, 
              kd = kd_dist)))'''
bins = [0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,
                               0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85]
#bins = [0,0.025,0.05,0.075,0.1]
eval_points = 1000
kd = kd_dist
probabilities = []

 
#print(get_probability(0.95, 0.975, 1000, kd))
## get probability to fall in each bin
for x in range(len(bins)):
    start_value = bins[x]
    end_value = bins[x+1] 
    probabilities.append((start_value, end_value, get_probability(start_value, end_value, eval_points, kd)))


frequency, bins = np.histogram(distances_merge, bins = 25)
for b, f in zip(bins[1:], frequency):
    print(round(b, 1), ' '.join(np.repeat('*', f)))

#[0,0.05, 0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85]
################# Get distance between two cluster representatives, get centroid of whole trajectory ###############
maxindex = np.where(distances_merge == np.amax(distances_merge))
print(len(distances_merge))
print(distances[1414,200])

beta = 1
index = np.exp(-beta*distances / distances.std()).sum(axis=1).argmax()
index = np.exp(-beta*distances[clusts['3']] / distances[clusts['3']].std()).sum(axis=1).argmax()
print(index)
print(clusts['3'][index])
     
centroid = traj[index]
centroid_merge = traj_merge[index]
print(centroid)
print(centroid_merge)

################# Create clusters ########################


#not above. 
print(distances - distances.T)
#expect assert error. Don't know why.
#assert np.all(distances - distances.T < 1e-6)
reduced_distances = squareform(distances, checks=False)
#reduced_distances_merge = squareform(distances_merge, checks=False)
#print(len(reduced_distances))

#implement the linkage algorithm
linkage = scipy.cluster.hierarchy.linkage(reduced_distances, method='ward')
#ward_linkage = scipy.cluster.hierarchy.linkage(reduced_distances, method='ward')
#linkage_merge = scipy.cluster.hierarchy.linkage(reduced_distances_merge, method='average')

#plot resulting dendrogram
plt.figure(dpi=600)
plt.title('RMSD Ward linkage hierarchical clustering, 5da0 merged')
plt.ylabel('Distance Score')
#plt.xlabel('Sample index or (Cluster size)')
#change this for max-clusts/truncation of dendrograms
##### dendrogram of all frames ####
scipy.cluster.hierarchy.dendrogram(linkage,
                                       count_sort='descendent',
                                       leaf_rotation = 90,
                                       leaf_font_size=12.,
                                       show_contracted=True)
#plt.savefig('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/1-24-2020/scipy-average-clustering-minimal-26clusts.png')

###### dendrogram of p_count size #######
p_count =4
plt.figure(dpi = 600)
plt.ylabel('Distance Score')
fancy_dendrogram(
    linkage,
    truncate_mode='lastp',
    count_sort='descendent',
    p=p_count,
    leaf_rotation=90.,
    leaf_font_size=12.,
    show_contracted=True,
    annotate_above=.2,)
#plt.savefig('C:/Users/skidmore/Desktop/Prestin-MD/skidmore/merged-traj-analysis/scipy-average-clustering-8clusts.png', dpi =600)
plt.title('RMSD Ward linkage hierarchical clustering, 5da0 merged')
plt.show()
#print(linkage)
#print(len(linkage))




### K-means
whitened = whiten(distances)
#klink2, distortion = scipy.cluster.vq.kmeans(whitened, 2)
klink, distor = scipy.cluster.vq.kmeans(whitened, 3)
plt.figure(dpi =600)
plt.scatter(range(len(whitened[:, 1])), whitened[:, 1])
plt.scatter(range(len(klink[:, 1])), klink[:, 1], c='r')
plt.show()

#####
whitened = whiten(distances)
klink, distor = scipy.cluster.vq.kmeans(distances, 5)
plt.figure(dpi =600)
plt.title('K-means PRMSD by time, production run 3')
plt.ylabel('PRMSD')
plt.xlabel('Frame number')
plt.scatter(range(len(distances[:, 1])), distances[:, 1])
plt.scatter(range(len(klink[:, 1])), klink[:, 1], c='r')
plt.show()

#k-means try again
whitened = whiten(a)
book = np.array(whitened[0, whitened[2]])
scipy.cluster.vq.kmeans(whitened,book)


##k-means again
from pylab import plot,show
from numpy import vstack,array
from numpy.random import rand
from scipy.cluster.vq import kmeans,vq

centroids,_ = kmeans(distances,4)
# assign each sample to a cluster
idx,_ = vq(distances,centroids)

# some plotting using numpy's logical indexing
plot(distances[idx==0,0],distances[idx==0,1],'ob',
     distances[idx==1,0],distances[idx==1,1],'or')
plot(centroids[:,0],centroids[:,1],'sg',markersize=8)
show()


##Possible alternative, if we treat our data as 1D instead of 2D
#Kernal density estimate


# Make labels for X-axis
'''RRR = scipy.cluster.hierarchy.dendrogram(linkage, truncate_mode='lastp',p=p_count, count_sort='descendent', no_plot=True)
print(RRR['ivl'])
print(RRR)
labels = [str(i) for i in range(p_count)]
for ii in range(len(RRR["ivl"])):
    #print(ii)
    #print(RRR['ivl'][ii])
    #print(labels[ii])

    if str(clusts.keys()) == in str(RRR['ivl'][ii]):
        print(clusts.values(), RRR['ivl'][ii], RRR['ivl'])
# create a label dictionary
temp = {RRR["ivl"][ii]: clusts.keys(ii) for ii in range(len(RRR["ivl"]))}
def llf(xx):
    return "{} - custom label!".format(temp[xx])
'''



#fcluster retrival of clusters from dendogram
from scipy.cluster.hierarchy import fcluster
max_d = .5
clusters = fcluster(linkage, max_d, criterion='distance')
clusters


from scipy.cluster.hierarchy import fcluster
max_d = .4
clusters = fcluster(ward_linkage, max_d, criterion='distance')
clusters
##ALSO CAN DO THIS WAY FOR KNOWN p_count # OF CLUSTERS



##This is the default method, it is not ideal, visualization of dendrograms is likely better
####inco_clusts = fcluster(linkage, t=1)


##Here is the 'elbow' method. I have seen this used before and know some people like it, though it may not be reliable
last = linkage[-10:, 2]
last_rev = last[::-1]
idxs = np.arange(1, len(last) + 1) 
plt.figure(dpi = 600)
plt.title('Elbow method for cluster number selction, Ward method, 5da0 merged')
plt.ylabel('Variance explained (0 = 100%)')
plt.xlabel('Number of clusters')
plt.plot(idxs, last_rev, label= 'Highest order clusters by distance score')

acceleration = np.diff(last, 2)  # 2nd derivative of the distances
acceleration_rev = acceleration[::-1]
plt.plot(idxs[:-2] + 1, acceleration_rev, label = 'Second Derivative of distances')
plt.legend(loc='upper right')
plt.show()
elbows = acceleration_rev.argmax()   # if idx 0 is the max of this we want 5 clusters
print("clusters:", elbows)



##Using the Silhouette method THIS IS WITH scikit-learn

#distmatrix1 = scipy.spatial.distance.squareform(distances + distances.T, checks=False)
#ddgm = scipy.cluster.hierarchy.linkage(distmatrix1, method="average")
#nodes = scipy.cluster.hierarchy.fcluster(linkage, p_count, criterion="maxclust")
###change the p_count in max_clusters to view the silhouette score for each cluster size

max_clusts = fcluster(linkage, 4, criterion='maxclust')
metrics.silhouette_score(distances + distances.T , max_clusts, metric='euclidean', sample_size=len(max_clusts))

max_clusts_merge = fcluster(linkage_merge, p_count, criterion='maxclust')
metrics.silhouette_score(distances_merge + distances_merge.T , max_clusts_merge, metric='euclidean', sample_size=len(max_clusts_merge))

#plot sihouette scores against numb of clusters
sil = np.array([[2,0.6196653899100691],[3,0.5952386041856075], [4,0.5554474350497706],[5,0.5553622756190352],[6,0.4661076384925495],[7,0.49468614189854687],
                [8,0.39108567655740356],[9,0.3858865043660024],[10,0.37371800536840294],[20,0.28545191449953344],[25,0.24694104169154832]])
plt.figure(dpi =600)
plt.title('Silhouette scoring for Ward cluster count assertion, 5da0 Merged')
plt.ylabel('Silhouette score')
plt.xlabel('Number of clusters')
z = np.polyfit(sil[:,0], sil[:,1], 1)
p = np.poly1d(z)
plt.scatter(sil[:,0], sil[:,1])
plt.plot(sil[:,0], p(sil[:,0]))


##OK now lets try to identify the 'centroid' of each of the clusters from the scipy method

clusts ={}
numb = 0 
for x in max_clusts:
    #print(x)
    if str(x) in clusts:
        clusts[str(x)].append(numb)
    elif str(x) not in clusts:
        clusts[str(x)] = [numb]
    numb += 1
        #print(clusts)
#traj = traj_merge        
centroids = []        
for k in clusts:
    centroids.append(traj[clusts[k]])        
#print(len(clusts[str(3)]))       
#print(i for i in list(i for i in range(p_count)))
print('Cluster sizes:')
for i in range(1,p_count+1):
    print(str(i),len(clusts[str(i)]))
    
## find the centroid(not really centroid, but the structure with the lowest average rmsd to others of its group) of each grouping
avg_distance = {}
for k, v in clusts.items():
    minimums = {}
    for val in v:
        #print(k, val, v)
        #print(distances[val,v])
        #print(np.mean(distances[val,v]))
        minimums[val] = np.mean(distances[val,v])
    minuma = min(minimums, key=minimums.get)
    avg_distance[k] = minuma
print(avg_distance)

print(distances[2087,4711])



print(distances)
print(len(distances))
print(distances[2000,2001])
xs = []
f = []
## Time sorted distribution of PRMSD
for x in range(len(distances)):
    xs.append(x)
    f.append((distances[x, x+1]))
    #print(x, distances[x, x+1])
count =0 
flist = []
bunches =[]
count2 = 0
xss = []
for i in f:
    count2 +=1
    if count <4:
        flist.append(i)
        count +=1
    
    else:
        bunches.append(np.mean(flist))
        flist = []
        count = 0
        xss.append(xs[count2])
        
#f.pop(4000)        
xs.pop(-1)    

print(np.max(f))
'''plt.figure(dpi =600)
plt.title('PRMSD by time, production run 3, averaged every 5')
plt.ylabel('PRMSD')
plt.xlabel('Frame number')
z = np.polyfit()

plt.scatter(xss,bunches)
plt.plot()'''

plt.figure(dpi =600)
plt.title('PRMSD by time, TM domain only')
plt.ylabel('PRMSD')
plt.xlabel('Frame number')
#z = np.polyfit(xs,f, 10)
#z = np.polyfit(xs[0:1000], f[0:1000], 4)
z = np.polyfit(xs[0:603], f[0:603], 4)
p = np.poly1d(z)
#plt.scatter(xs[0:1000], f[0:1000])
plt.scatter(xs[0:603], f[0:603])
#plt.plot(xs[0:1000], p(xs[0:1000]), c='red')
plt.plot(xs[0:603], p(xs[400:800]), c='red')

print(np.min(f[450:650]))
fffd = f[450:650]
print(np.mean(f[450:650]))

np.std(f[450:650])


## production run 3 minimum 450 -> 650
##Minimum: 0.12164866179227829 , frame 538
## 
np.savetxt("5da0_allatoms_PRMSD.csv", distances, delimiter=",")
#d2 = distances.tolist()
import scipy.io as sio
adict = {}
adict['distances'] = distances
sio.savemat("5da0_allatoms_PRMSD.mat", adict)

plt.figure(dpi =600)
plt.title('PRMSD by time, production run 3, averaged every 5')
plt.ylabel('PRMSD')
plt.xlabel('Frame number')
z = np.polyfit(xss,bunches, 10)
p = np.poly1d(z)
plt.scatter(xss, bunches)
plt.plot(xss, p(xss), c='red')



###FFT of time sorted PRMSD
import matplotlib.pyplot as plotter
plt(f, sp)
f1 = np.array(f)
sp = np.fft.fft(f1)/(len(f1))
sp2 = sp[range(int(len(f1)/2))]

tpCount     = len(f1)
values      = np.arange(int(tpCount/2))
timePeriod  = tpCount/1
frequencies = values/timePeriod
sp3 = sp2[1:]
plt.plot(abs(sp3))


f1 = f1[1:]
sp2 = sp[1:]
freq = np.fft.fftfreq(f1.shape[-1])
plt.plot(freq, sp.real, freq, sp.imag)
plt.show()


t = np.arange(256)
ts = np.sin(t)
ss = np.fft.fft(np.sin(t))




## cannot figure out how to extract the centroid structures from the hierarchy.
#T = scipy.cluster.hierarchy.fcluster(linkage, .3, criterion = 'inconsistent')
#plt.scatter(linkage[:,0],linkage[:,1])















##Here is an attempt using scikitlearn. 

import time as time
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering

st = time.time()
n_clusters = 5  # number of regions
#ward = AgglomerativeClustering(n_clusters=n_clusters, linkage='ward')
#INITIATE WARD CLUSTERING on the linkage with n_clusters number of clusters
# group trajectoryies by frame into their clustered groups
fw = AgglomerativeClustering(n_clusters=n_clusters, linkage='average').fit_predict(linkage)

print(list(fw))
plt.title('Ward Agglomerative Clsutering by RMSD hierarchical clustering')
plt.xlabel('Trajectory frame (Taken every 500ps)')
plt.ylabel('Cluster assignment (26 options)')
plt.scatter(range(len(fw)),fw)
plt.show()

clusts ={}
numb = 0 
fw = list(fw)
for x in fw:
    print(x)
    numb += 1
    if str(x) in clusts:
        clusts[str(x)].append(numb)
    elif str(x) not in clusts:
        clusts[str(x)] = [numb]
        print(clusts)
        
centroids = []        
for k in clusts:
    centroids.append(traj[clusts[k]])        
print(len(clusts[str(2)]))       
print(i for i in list(i for i in range(5)))
for i in range(5):
    print(len(clusts[str(i)]))
    
## find the centroid(not really centroid, but the structure with the lowest average rmsd to others of its group) of each grouping
avg_distance = {}
for k, v in clusts.items():
    minimums = {}
    for val in v:
        #print(k, val, v)
        #print(distances[val,v])
        #print(np.mean(distances[val,v]))
        minimums[val] = np.mean(distances[val,v])
    minuma = min(minimums, key=minimums.get)
    avg_distance[k] = minuma
print(avg_distance)














'''def plot_prob_density(df_lunch, df_dinner, field, x_start, x_end):
    plt.figure(figsize = (10, 7))

    #unit = 1.5
    #x = np.linspace(df_lunch.min() - unit, df_lunch.max() + unit, 1000)[:, np.newaxis]

    # Plot the data using a normalized histogram
    plt.hist(df_lunch, bins=10, density=True, label='Lunch Time', color='orange', alpha=0.2)
    plt.hist(df_dinner, bins=10, density=True, label='Dinner Time', color='navy', alpha=0.2)

    # Do kernel density estimation
    kd_lunch = KernelDensity(kernel='gaussian', bandwidth=0.5).fit(df_lunch)
    kd_dinner = KernelDensity(kernel='gaussian', bandwidth=0.5).fit(df_dinner)

    # Plot the estimated densty
    kd_vals_lunch = np.exp(kd_lunch.score_samples(x))
    kd_vals_dinner = np.exp(kd_dinner.score_samples(x))

    plt.plot(x, kd_vals_lunch, color='orange')
    plt.plot(x, kd_vals_dinner, color='navy')
    
    plt.axvline(x=x_start,color='red',linestyle='dashed')
    plt.axvline(x=x_end,color='red',linestyle='dashed')

    # Show the plots
    plt.xlabel(field, fontsize=15)
    plt.ylabel('Probability Density', fontsize=15)
    plt.legend(fontsize=15)
    plt.show()
    gc.collect()
    return kd_lunch, kd_dinner

def get_probability(start_value, end_value, eval_points, kd):
    
    # Number of evaluation points 
    N = eval_points                                      
    step = (end_value - start_value) / (N - 1)  # Step size

    x = np.linspace(start_value, end_value, N)[:, np.newaxis]  # Generate values in the range
    kd_vals = np.exp(kd.score_samples(x))  # Get PDF values for each x
    probability = np.sum(kd_vals * step)  # Approximate the integral of the PDF
    return probability.round(4)'''
    






'''
distances = np.empty((traj.n_frames, traj.n_frames))
for k in clusts:    
    for i in clusts[k]: 
        #print(type(i))
        distances[int(k)] = md.rmsd(traj[i], traj[clusts[k]], int(i), atom_indices=atom_indices)
    print('Max pairwise rmsd: %f nm' % np.max(distances))
import statistics
beta = 1
print(beta*distances)
print(distances[0,1])
index = {}
vallst = []
for k, v in clusts.items():
    vallst = list()
    for val in v:
        print(k, val, v)
        print(distances[val,v])
        dv = distances[val,v]
        for d in dv:
            vallst.append(d)
        #index[k] = np.exp(-beta*distances[val,v] / distances[val,v].std()).sum().argmax()
        #print(vallst, type(vallst), statistics.stdev(vallst), type(statistics.stdev(vallst)))
    sdvl = statistics.stdev(vallst)
    index[k] = np.exp((x*-beta / sdvl for x in vallst)).sum(axis=1).argmax()

#index = np.exp(-beta*distances[clusts[k]] / distances[clusts[k]].std()).sum(axis=1).argmax()
print(index)

centroid = traj[index]
print(centroid)

print((-beta*distances[val,v] / distances[val,v].std()))
'''





'''       
print("Elapsed time: ", time.time() - st)
print("Number of clusters: ", np.unique(label).size)





##Try to use 
# Sum the vectors in each cluster
lens = {}      # will contain the lengths for each cluster
centroids = {} # will contain the centroids of each cluster
features = ()
for idx,clno in enumerate(T):
    centroids.setdefault(clno,np.zeros(2)) 
    centroids[clno] += features[idx,:]
    lens.setdefault(clno,0)
    lens[clno] += 1
# Divide by number of observations in each cluster to get the centroid
for clno in centroids:
    centroids[clno] /= float(lens[clno])


#maybe?? This doesn't really work because the x[0] values don't correspond to trajectory frames??
clusts ={}
for x in linkage:
    group = x[3]
    if group in clusts:
        clusts[group].append(x[0])
    elif group not in clusts:
        clusts[group] = [x[0]]

'''









