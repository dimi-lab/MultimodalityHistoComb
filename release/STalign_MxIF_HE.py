# Refer to example code within STalign
# https://github.com/JEFworks-Lab/STalign/blob/main/docs/notebooks/xenium-heimage-alignment.ipynb
# Tested on Python=3.11.11
# Other dependencies are described in STalign original repo.
# 
# Please Note: This method requires to provide pairs of landmarks: pointsI and pointsJ

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import torch
import tifffile as tf
from STalign import STalign
from utils import get_cell_loc

HE_img_fn = "/Users/jjiang10/Data/OV_TMA/HE_H-13.tif"
HE_quant_fn = "/Users/jjiang10/Data/OV_TMA/AlignedCellQuant/H-13_1on1_HE_quant.tsv"
MxIF_quant_fn = "/Users/jjiang10/Data/OV_TMA/AlignedCellQuant/H-13_1on1_MxIF_quant.tsv" 

HE_centroids = get_cell_loc(HE_quant_fn)   # you may need to specify column names to get the cell coordinates
MxIF_centroids = get_cell_loc(MxIF_quant_fn)

xM = HE_centroids[:, 0]
yM = HE_centroids[:, 1]
xN = MxIF_centroids[:, 0]
yN = MxIF_centroids[:, 1]

he_img = tf.TiffFile(HE_img_fn).pages[0].asarray().astype(np.float)  # read the images

Inorm = STalign.normalize(he_img)
I = Inorm.transpose(2,0,1)
YI = np.array(range(I.shape[1]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
XI = np.array(range(I.shape[2]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
extentI = STalign.extent_from_x((YI,XI))


fig,ax = plt.subplots()
ax.scatter(xM,yM,s=1,alpha=0.2)
ax.scatter(xN,yN,s=1,alpha=0.1)
ax.set_title("Plot1")
ax.set_aspect('equal')


XJ,YJ,M,_ = STalign.rasterize(xM, yM, dx=30)
J = np.vstack((M, M, M)) # make into 3xNxM
# normalize
J = STalign.normalize(J)

# Randomly select 4 points from HE_centroids and the corresponding points from MxIF_centroids
indices = np.random.choice(len(xM), 4, replace=False)
selected_HE_points = HE_centroids[indices]
selected_MxIF_points = MxIF_centroids[indices]

temp = selected_HE_points[:,0].copy()
selected_HE_points[:,0]  = selected_HE_points[:,1]
selected_HE_points[:,1]  =  temp

temp = selected_MxIF_points[:,0].copy()
selected_MxIF_points[:,0]  = selected_MxIF_points[:,1]
selected_MxIF_points[:,1]  =  temp

pointsI = selected_HE_points
pointsJ = selected_MxIF_points

extentJ = STalign.extent_from_x((YJ,XJ))

L,T = STalign.L_T_from_points(pointsI,pointsJ)

# note points are as y,x
affine = np.dot(np.linalg.inv(L), [yM - T[0], xM - T[1]]) 
print(affine.shape)
xMaffine = affine[0,:] 
yMaffine = affine[1,:] 

fig,ax = plt.subplots()
ax.scatter(xMaffine,yMaffine,s=1,alpha=0.2)
ax.scatter(xN,yN,s=1,alpha=0.1)
ax.set_title("Plot2")
ax.set_aspect('equal')
# set device for building tensors
if torch.cuda.is_available():
    torch.set_default_device('cuda:0')
else:
    torch.set_default_device('cpu')

#  run LDDMM
# specify device (default device for STalign.LDDMM is cpu)
if torch.cuda.is_available():
    device = 'cuda:0'
else:
    device = 'cpu'

# keep all other parameters default
params = {'L':L,'T':T,
          'niter':2000,
          'pointsI':pointsI,
          'pointsJ':pointsJ,
          'device':device,
          'sigmaM':0.15, 
          'sigmaB':0.10,
          'sigmaA':0.11,
          'epV': 10,
          'muB': torch.tensor([0,0,0]), # black is background in target
          'muA': torch.tensor([1,1,1]) # use white as artifact 
          }

out = STalign.LDDMM([YI,XI],I,[YJ,XJ],J,**params)
A = out['A']
v = out['v']
xv = out['xv']
phi = STalign.build_transform(xv,v,A,XJ=[YI,XI],direction='f')
tpointsIphiiJ = STalign.transform_image_target_to_source(xv,v,A,[YJ,XJ],J,[YI,XI])
phiipointsJ = STalign.transform_points_target_to_source(xv,v,A,pointsJ)
tpointsJ = STalign.transform_points_target_to_source(xv,v,A,np.stack([yN, xN], -1))
tpointsI = STalign.transform_points_source_to_target(xv,v,A,np.stack([yM, xM], -1))
if tpointsI.is_cuda:
    tpointsI = tpointsI.cpu()

fig,ax = plt.subplots()
ax.scatter(tpointsI[:, 0],tpointsI[:, 1],s=2,alpha=0.7)
ax.scatter(yN,xN,marker='^', s=2,alpha=0.3)
ax.set_title("Transformed H&E and MxIF cell centroids")
ax.set_aspect('equal')
plt.show()

print("Done")

