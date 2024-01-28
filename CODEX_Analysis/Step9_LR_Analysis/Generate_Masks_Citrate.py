import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tifffile

sample_names = ['NBM67_H38_Citrate','NBM67_H41_Citrate','NBM67_H14_Citrate',
'NBM70_H10_Citrate','NBM70_H26_Citrate','NBM70_H32_Citrate']
segmentations = ['NBM67/mesmer/NBM67_segmentation_cell.tif','NBM67/mesmer/NBM67_segmentation_cell.tif','NBM67/mesmer/NBM67_segmentation_cell.tif',
'NBM70/mesmer/NBM70_segmentation_cell.tif','NBM70/mesmer/NBM70_segmentation_cell.tif','NBM70/mesmer/NBM70_segmentation_cell.tif']

for name, mask in zip(sample_names, segmentations):
	print(name)
	print(mask)
	clusters = 'Reference_Map_Masks/'+name+'_CODEX_ReferenceMap.csv' 
	output = 'Reference_Map_Masks/'+name+'_ReferenceMap_Mask.npy'
	output2 = 'Reference_Map_Masks/'+name+'_ReferenceMap_Mask.tiff'
	print('Reading Clusters')
	clusters = pd.read_csv(clusters)
	print('Reading Mask')
	mask = tifffile.imread(mask)

	# Initialize numpy array from segmentation mask
	print('Initializing numpy array')
	masks_out = mask.copy()
	masks_out = np.array(masks_out).astype(int)
	print(masks_out.dtype)

	#Convert those cells that have been filtered out to zero
	print('Filtering cells')
	filter_mask = np.isin(masks_out,clusters['CellID'], invert=True)
	masks_out[filter_mask] = 0

	# Create dictionary mapping Cell ID keys to cluster values
	print('Creating mapping')
	mapping = dict(zip(clusters["CellID"], clusters["Annotation"]))

	# Define function to change the values in input array by key-value pairs from dict
	def mp(entry):
	    return mapping[entry] if entry in mapping.keys() else entry

	print('Vectorizing array')
	mp = np.vectorize(mp)

	print('Converting Values') 
	out = mp(masks_out)

	print('Saving mask')
	np.save(output, out)
	tifffile.imwrite(output2, out)
