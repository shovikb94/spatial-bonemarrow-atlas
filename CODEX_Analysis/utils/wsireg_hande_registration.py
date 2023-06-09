from wsireg.wsireg2d import WsiReg2D
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

samples = ['H41']
# samples = ['H26', 'H33', 'H35', 'H38', 'H39']

for sample in samples:
    he_path = "data/"+sample+"_HE.tif"
    codex_path = "data/"+sample+"_CODEX.tif"

    save_name = "BM"+sample
    project_folder = "BM"

    HE_resolution = 0.2485
    CODEX_resolution = 0.5069

    # HE_resolution = 0.37
    # CODEX_resolution = 0.37

    # initialize registration graph
    reg_graph = WsiReg2D(save_name, project_folder)

    # add registration images (modalities)
    reg_graph.add_modality(
        "CODEX",
        codex_path,
        image_res=CODEX_resolution,
        # channel_names=["DAPI", "CD44", "CD45"],
        # channel_colors=["blue", "green", "red"],
        preprocessing={
            "image_type": "FL",
            "ch_indices": [0],
            "as_uint8": True,
            "contrast_enhance": True,
        },
    )

    reg_graph.add_modality(
        "HE",
        he_path,
        image_res=HE_resolution,
        preprocessing={"image_type": "BF", "as_uint8": True, "invert_intensity": True},
    )

    reg_graph.add_reg_path(
        "HE",
        "CODEX",
        thru_modality=None,
        reg_params=["rigid", "affine"],
    )

    reg_graph.register_images()
    print("Registration done")
    reg_graph.save_transformations()
    print("Transformations saved")
    reg_graph.transform_images(file_writer="ome.tiff")

    print("Done")