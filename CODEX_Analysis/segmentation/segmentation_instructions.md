## Instructions for Segmentation ##
1. Use segmentation tuning script to appropriately rescale markers and create merged membrane pseudochannel from segmentation markers of interest
   (in our case Vimentin, CD45, and NaKATPase) and output two-channel TIF compatible with Mesmer.
2. Run mesmer on HPC using the Mesmer shell script calling the python file. Note for very large images you may need to tile your image and run piece-wise, or crop to ROI.
3. Run the quantification notebook to output the CSV from which to create the Seurat object in Step 1. 
