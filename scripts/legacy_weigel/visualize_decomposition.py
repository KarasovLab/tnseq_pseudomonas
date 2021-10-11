import numpy as np
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import pickle


def draw_vector(v0, v1, ax=None):
    ax = ax or plt.gca()
    arrowprops=dict(arrowstyle='->',
                    linewidth=2,
                    shrinkA=0, shrinkB=0)
    ax.annotate('', v1, v0, arrowprops=arrowprops)


def PCA(data, dims_rescaled_data=2):
    import numpy as NP
    from scipy import linalg as LA
    m, n = data.shape
    # mean center the data
    data -= data.mean(axis=0)
    # calculate the covariance matrix
    R = NP.cov(data, rowvar=False)
    # calculate eigenvectors & eigenvalues of the covariance matrix
    # use 'eigh' rather than 'eig' since R is symmetric, 
    # the performance gain is substantial
    evals, evecs = LA.eigh(R)
    # sort eigenvalue in decreasing order
    idx = NP.argsort(evals)[::-1]
    evecs = evecs[:,idx]
    # sort eigenvectors according to same index
    evals = evals[idx]
    # select the first n eigenvectors (n is desired dimension
    # of rescaled data array, or dims_rescaled_data)
    evecs = evecs[:, :dims_rescaled_data]
    # carry out the transformation on the data using eigenvectors
    # and return the re-scaled data, eigenvalues, and eigenvectors
    return NP.dot(evecs.T, data.T).T, evals, evecs


def plot_pca(data):
    from matplotlib import pyplot as MPL
    clr1 =  '#2026B2'
    fig = MPL.figure()
    ax1 = fig.add_subplot(111)
    data_resc, data_orig, what = PCA(data)
    ax1.plot(data_resc[:, 0], data_resc[:, 1], '.', mfc=clr1, mec=clr1)
    MPL.show()


def combine_host(pd_all):
    host_list=['Barley', 'NP29a', 'Canola', 'Col-0', 'Peas']
    pd_collapse=pd.DataFrame(index=host_list, columns=pd_all.columns)
    for host in host_list:
        relevant=[row for row in pd_all.index if host in row]
        sum_rows=pd_all.ix[relevant].sum(axis=0)
        pd_collapse.ix[host]=sum_rows


pd_all=pickle.load(open("/ebio/abt6_projects9/Pseudomonas_diversity/Tnseq/processed_reads/hiseq0081/pd_all.cpk", 'r'))
X=pd_all.ix[[rec for rec in pd_all.index if "317" not in rec and "NP29" not in rec]]
myData=np.array(X, dtype="float").T
from matplotlib.mlab import PCA
results = PCA(myData)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x = []
y = []
z = []
for item in results.Y:
 x.append(item[0])
 y.append(item[1])
 z.append(item[2])

plt.close('all') # close all latent plotting windows
fig1 = plt.figure() # Make a plotting figure
ax = Axes3D(fig1) # use the plotting figure to create a Axis3D object.
pltData = [x,y,z] 
ax.scatter(pltData[0], pltData[1], pltData[2], 'bo') # make a scatter plot of blue dots from the data
 
from sklearn.decomposition import PCA

X_perc=(X.T/X.T.sum()).T
pca = PCA(n_components = 2)
pca.fit(X_perc)
X_pca = pca.transform(X_perc)
projected = pca.fit_transform(X_perc)
print("original shape:   ", X.shape)
print("transformed shape:", X_pca.shape)

plt.scatter(projected[:, 0], projected[:, 1],
             edgecolor='none', alpha=0.5,
            cmap=plt.cm.get_cmap('spectral', 10))
plt.xlabel('component 1')
plt.ylabel('component 2')
plt.colorbar();
