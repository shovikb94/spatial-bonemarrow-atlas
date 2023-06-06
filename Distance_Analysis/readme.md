## workflow

1. Use `image_to_contour.py` to convert the annotation mask to contour points.
2. Use `run_health_celltype.py` , `run_health_neighbour.py` ,`run_aml_celltype.py` ,`run_aml_neighbour.py` to get the output data
3. Use codes in `combine_run.ipynb` to convert multiply .txt to a .csv table
4. Use codes in `summarize.ipynb` to get the final all-in-one table

## Input data

1. segmentation mask: example: `/mask` folder
2. Cells table: contain coordinates and labels (cell type, neighborhood)

