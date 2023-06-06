library(spatstat)
library(png)
library(stringr)
library(foreach)
library(doParallel)
library(parallel)


immune.combined <- read.csv('neighborhood.csv')

filenames <- c('H10.csv','H14.csv','H26.csv', 'H27.csv', 'H32.csv', 'H33.csv', 'H35.csv', 
                'H36.csv', 'H37.csv', 'H38.csv', 'H39.csv', 'H41.csv')
# extract the sample name from each filename
sample_names <- c('H10', 'H14', 'H26', 'H27', 'H32', 'H33', 'H35', 'H36', 'H37', 'H38', 'H39', 'H41')
sample_ids <- c("SB67_NBM27_H10_CODEX_Mesmer", "SB67_NBM28_H14_CODEX_Mesmer", 'SB67_NBM36_H26_CODEX_Mesmer',
                'SB67_NBM41_H27_CODEX_Mesmer', 'SB67_NBM31_H32_CODEX_Mesmer', 'SB67_NBM38_H33_CODEX_Mesmer',
                'SB67_NBM37_H35_CODEX_Mesmer', 'SB67_NBM33_H36_CODEX_Mesmer', 'SB67_NBM32_H37_CODEX_Mesmer',
                'SB67_NBM34_H38_CODEX_Mesmer', 'SB67_NBM40_H39_CODEX_Mesmer', 'SB67_NBM39_H41_CODEX_Mesmer')
neighborhood_names_zeroindex <- c('0' = "HSC / Mature Myeloid",
                        "1" = "Erythroid/Myeloid",
                        "2" = "PC/Arteriolar",
                        "3" = "Erythroid",
                        "4" = "Arteriolar",
                        "5" = "Erythroid",
                        "6" = "Lymphoid",
                        "7" = "Erythroid/Myeloid/Lymphoid",
                        "8" = "Early Myeloid / Endosteal",
                        "9" = "Myeloid/Lymphoid",
                        "10" = "HSPC/Intermediate Myeloid",
                        "11" = "Erythroid/Myeloid/Lymphoid",
                        "12" = "Erythroid/Myeloid",
                        "13" = "Early Myeloid / Arteriolar",
                        "14" = "Peri-Arterolar Lymphoid")
structures <- c('adipocyte', 'art', 'sinusoid', 'bone', 'macrophage',  'stroma')

n_permutations <- 100

cores=detectCores()
cl <- makeCluster(3) #not to overload your computer
registerDoParallel(cl)
foreach (s = 1:6) %dopar% {
    library(spatstat)
    library(stringr)
    structure <- structures[s]
    print(structure)
    coefficients_list_1 <- list()
    coefficients_list_2 <- list()
    coefficients_list_3 <- list()
    sample_list <- list()
    celltype_list <- list()
    distance_median_list <- list()
    distance_sd_list <- list()
    pop_list <- list()
    pvalue_list <- list()
    iqr_list <- list()
    for (i in 1:length(sample_names)) {
        sample_name <- sample_names[i]
        print(sample_name)
        filename <- filenames[i]
        sample_id <- sample_ids[i]
        sample_data <- subset(immune.combined, subset = orig.ident == sample_id)
        unique_nbs <- unique(sample_data$neighborhood10)
        unique_cell_types <- unique(sample_data$cluster_anno_l2)
        
        cts <- read.csv(paste0('h_contour/',structure,'_csv/', filename))
        
        # if structure is bone, sample 10%, if structure is macrophage or stroma, sample 1%, otherwise sample 20%
        if (structure == 'bone') {
            cts <- cts[sample(nrow(cts), 0.1*nrow(cts)), ]
        } else if (structure == 'macrophage' | structure == 'stroma') {
            cts <- cts[sample(nrow(cts), 0.01*nrow(cts)), ]
        } else {
            cts <- cts[sample(nrow(cts), 0.2*nrow(cts)), ]
        }

        xmin <- min(cts$x0)
        xmax <- max(cts$x1)
        ymin <- min(cts$y0)
        ymax <- max(cts$y1)
        sample_data_crop <- subset(sample_data, subset = x.coord > xmin & x.coord < xmax & y.coord > ymin & y.coord < ymax)
        sample_data_crop_xy <- sample_data_crop[, c('x.coord', 'y.coord')]
        for (j in 1:length(unique_nbs)) {
            cell_type <- unique_nbs[j]
            # print(cell_type)
            coords_crop <- subset(sample_data_crop, subset = neighborhood10 == cell_type)
            # coords_crop <- subset(coords, subset = x.coord > xmin & x.coord < xmax & y.coord > ymin & y.coord < ymax)
            x <- coords_crop$x.coord
            y <- coords_crop$y.coord
            population <- length(x)
            if (length(x) < 10) {
                next
            }
            
            cells <- ppp(x, y, c(xmin-10, xmax+10), c(ymin-10, ymax+10))
            bone_cts <- psp(cts$x0, cts$y0, cts$x1, cts$y1, window = cells$window)
            
            bm <-list(cell=cells,bone=bone_cts)
            bm$dbone <- distfun(bm$bone)

            f <- ppm(cell~ dbone, data=bm)
            
            distance <- bm$dbone(x,y)
            # remove NA values
            distance <- distance[!is.na(distance)]
            # remove infinite values
            distance <- distance[!is.infinite(distance)]
            distance_median <- median(distance)
            # distance_sd <- sd(distance)
            # compute the std of the distance from 25th to 75th percentile
            distance_sd <- sd(distance[which(distance > quantile(distance, 0.25) & distance < quantile(distance, 0.75))])
            iqr <- quantile(distance, 0.75) - quantile(distance, 0.25)
            # perform permutation test
            p_distance_median <- list()
            for (k in 1:n_permutations) {
                # print(k)
                # sample cell coordinates
                permuted_xy <- sample_data_crop_xy[sample(nrow(sample_data_crop_xy), population), ]
                permuted_x <- permuted_xy[, 1]
                permuted_y <- permuted_xy[, 2]
                p_distance <- bm$dbone(permuted_x, permuted_y)
                p_distance <- p_distance[!is.na(p_distance)]
                p_distance <- p_distance[!is.infinite(p_distance)]
                p_distance_median <- append(p_distance_median, median(p_distance))
            }
            
            coefficients_list_1 <- append(coefficients_list_1, summary(f)$coef$Estimate[1])
            coefficients_list_2 <- append(coefficients_list_2, summary(f)$coef$Estimate[2])
            sample_list <- append(sample_list, sample_name)
            celltype_list <- append(celltype_list, neighborhood_names_zeroindex[as.character(cell_type)])
            pop_list <- append(pop_list, population)
            distance_median_list <- append(distance_median_list, distance_median)
            distance_sd_list <- append(distance_sd_list, distance_sd)
            pvalue_list <- append(pvalue_list, sum(p_distance_median < distance_median) / n_permutations)
            iqr_list <- append(iqr_list, iqr)
            }
        # }
    }
    # save to dsv
    sample_vector <- as.character(sample_list)
    writeLines(sample_vector,paste0( "output_f/nb_sample_list_", structure, '.txt'))
    celltype_vector <- as.character(celltype_list)
    writeLines(celltype_vector, paste0("output_f/nb_celltype_list_", structure, '.txt'))
    coefficients_list_1_vector <- as.character(coefficients_list_1)
    writeLines(coefficients_list_1_vector, paste0("output_f/nb_coefficients_list_1_", structure, '.txt'))
    coefficients_list_2_vector <- as.character(coefficients_list_2)
    writeLines(coefficients_list_2_vector, paste0("output_f/nb_coefficients_list_2_", structure, '.txt'))
    pop_vector <- as.character(pop_list)
    writeLines(pop_vector, paste0("output_f/nb_pop_list_", structure,  '.txt'))
    distance_median_vector <- as.character(distance_median_list)
    writeLines(distance_median_vector, paste0("output_f/nb_distance_mean_list_", structure, '.txt'))
    distance_sd_vector <- as.character(distance_sd_list)
    writeLines(distance_sd_vector, paste0("output_f/nb_distance_sd_list_", structure, '.txt'))
    pvalue_vector <- as.character(pvalue_list)
    writeLines(pvalue_vector, paste0("output_f/nb_pvalue_list_", structure, '.txt'))
    iqr_vector <- as.character(iqr_list)
    writeLines(iqr_vector, paste0("output_f/nb_iqr_list_", structure, '.txt'))
}
