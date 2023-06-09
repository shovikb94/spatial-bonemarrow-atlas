#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np
from skimage import io
from matplotlib import pyplot as plt
#from deepcell.datasets import multiplex_tissue
from deepcell.utils.plot_utils import create_rgb_image
from PIL import Image
import pandas as pd
#import tifffile
#import imageio
import glob
from deepcell.applications import Mesmer
from deepcell.utils.plot_utils import make_outline_overlay


# In[4]:


mcmicroDir = '/mnt/isilon/tan_lab_imaging/FUSION/NBM27_H10_CITRATE_REIMAGE/H10/Scan1'
inOmeTiff = glob.glob(mcmicroDir + '/mesmer/*_nuclear_panmembrane.ome.tif')[0]

outDir = os.path.join(mcmicroDir, 'mesmer')
sampleName = os.path.basename(mcmicroDir)
#outFile = os.path.join(outDir, sampleName + '_nuclear_panmembrane.ome.tif')

print('Input OME.TIFF file: ' + inOmeTiff)
#print('Output OME.TIFF file: ' + outFile)
outFileCell    = os.path.join(outDir, sampleName + '_segmentation_cell.tif')
outFileNuclear = os.path.join(outDir, sampleName + '_segmentation_nuclear.tif')


# In[5]:


# Load pre-built Mesmer segmentation model
from deepcell.applications import Mesmer
app = Mesmer()
print('Training Resolution:', app.model_mpp, 'microns per pixel')


# In[ ]:


# Read stacked nuclear and pan-membrane OME.TIF and plot the image to verify data has been loaded correctly

im = io.imread(inOmeTiff)
X_train = np.transpose(im, (1,2,0))
X_train = X_train.reshape(1, X_train.shape[0], X_train.shape[1], X_train.shape[2])
rgb_images = create_rgb_image(X_train, channel_colors=['green', 'blue'])

### Segment cells by Mesmer and save cellular segmentation mask
print('Starting Cell Segmentation')
segmentation_predictions = app.predict(X_train, image_mpp=0.5)
im = Image.fromarray(segmentation_predictions[0,:,:,0]).save(outFileCell)
overlay_data = make_outline_overlay(rgb_data=rgb_images, predictions=segmentation_predictions)
fig, ax = plt.subplots(1, 2, figsize=(15, 15))
ax[0].imshow(rgb_images[0, ...])
ax[1].imshow(overlay_data[0, ...])
ax[0].set_title('Raw data')
ax[1].set_title('Cellular Predictions')
fig.savefig(os.path.join(outDir + '_segmentation_cell_overlay.tif'), dpi=300)
np.save("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/SB67_NBM27_H10_CITRATE_REIMAGE/mesmer/segmentation_predictions.npy",segmentation_predictions)
### Segment nuclei by Mesmer and save nuclear segmentation mask
print('Starting Nuclear Segmentation')
segmentation_predictions_nuc = app.predict(X_train, image_mpp=0.5, compartment='nuclear')
im = Image.fromarray(segmentation_predictions_nuc[0,:,:,0]).save(outFileNuclear)
overlay_data_nuc = make_outline_overlay(rgb_data=rgb_images, predictions=segmentation_predictions)
fig, ax = plt.subplots(1, 2, figsize=(15, 15))
ax[0].imshow(rgb_images[0, ...])
ax[1].imshow(overlay_data[0, ...])
ax[0].set_title('Raw data')
ax[1].set_title('Nuclear Predictions')
fig.savefig(os.path.join(outDir + '_segmentation_nuclear_overlay.tif'), dpi=300)
np.save("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/SB67_NBM27_H10_CITRATE_REIMAGE/mesmer/segmentation_predictions_nuclear.npy",segmentation_predictions_nuc)
# plt.close('all')


# In[ ]:





# In[ ]:




