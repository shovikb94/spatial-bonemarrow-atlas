library(spatstat)
library(png)
library(stringr)
library(foreach)
library(doParallel)
library(parallel)

immune.combined <- read.csv('AML_NSM_combined.neighborhood.csv')

sample_names <- c('183', '191', '382', '1329', '1443')
filenames <- c('183.csv', '191.csv', '382.csv', '1329.csv', '1443.csv')
sample_ids <- c('SB67_NBM51_AML1_183_CODEX_Mesmer', 'SB67_NBM44_AML2_191_CODEX_Mesmer', 'SB67_NBM46_AML1_382_CODEX_Mesmer', 'SB67_NBM52_AML3_1329_CODEX_Mesmer', 'SB67_NBM54_AML3_1443_CODEX_Mesmer')

structures <- c('bone','macrophage',  'stroma')

n_permutations <- 100
cl <- makeCluster(3) #not to overload your computer
registerDoParallel(cl)
foreach (s = 1:3) %dopar% {
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
        unique_cell_types <- unique(sample_data$classified_cluster_anno_l2)
        
        cts <- read.csv(paste0('aml_contour/',structure,'_csv/', filename))
        if (structure == 'bone') {
            cts <- cts[sample(nrow(cts), 0.5*nrow(cts)), ]
        } else if (structure == 'macrophage' | structure == 'stroma') {
            cts <- cts[sample(nrow(cts), 0.1*nrow(cts)), ]
        } 
        xmin <- min(cts$x0)
        xmax <- max(cts$x1)
        ymin <- min(cts$y0)
        ymax <- max(cts$y1)
        sample_data_crop <- subset(sample_data, subset = x.coord > xmin & x.coord < xmax & y.coord > ymin & y.coord < ymax)
        sample_data_crop_xy <- sample_data_crop[, c('x.coord', 'y.coord')]
        for (j in 1:length(unique_cell_types)) {
            cell_type <- unique_cell_types[j]
            # print(cell_type)
            coords_crop <- subset(sample_data_crop, subset = classified_cluster_anno_l2 == cell_type)
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

            # compute the interquartile range
            distance_iqr <- quantile(distance, 0.75) - quantile(distance, 0.25)

            # perform permutation test
            p_distance_median <- list()
            for (k in 1:n_permutations) {
                # print(k)
                # permute the cell coordinates
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
            celltype_list <- append(celltype_list, as.character(cell_type))
            pop_list <- append(pop_list, population)
            distance_median_list <- append(distance_median_list, distance_median)
            distance_sd_list <- append(distance_sd_list, distance_sd)
            pvalue_list <- append(pvalue_list, sum(p_distance_median < distance_median) / n_permutations)
            iqr_list <- append(iqr_list, distance_iqr)
            }
        # }
    }

    sample_vector <- as.character(sample_list)
    writeLines(sample_vector,paste0( "output_f/aml_sample_list_", structure, '.txt'))
    celltype_vector <- as.character(celltype_list)
    writeLines(celltype_vector, paste0("output_f/aml_celltype_list_", structure, '.txt'))
    coefficients_list_1_vector <- as.character(coefficients_list_1)
    writeLines(coefficients_list_1_vector, paste0("output_f/aml_coefficients_list_1_", structure, '.txt'))
    coefficients_list_2_vector <- as.character(coefficients_list_2)
    writeLines(coefficients_list_2_vector, paste0("output_f/aml_coefficients_list_2_", structure, '.txt'))
    pop_vector <- as.character(pop_list)
    writeLines(pop_vector, paste0("output_f/aml_pop_list_", structure,  '.txt'))
    distance_median_vector <- as.character(distance_median_list)
    writeLines(distance_median_vector, paste0("output_f/aml_distance_median_list_", structure, '.txt'))
    distance_sd_vector <- as.character(distance_sd_list)
    writeLines(distance_sd_vector, paste0("output_f/aml_distance_sd_list_", structure, '.txt'))
    pvalue_vector <- as.character(pvalue_list)
    writeLines(pvalue_vector, paste0("output_f/aml_pvalue_list_", structure, '.txt'))
    iqr_vector <- as.character(iqr_list)
    writeLines(iqr_vector, paste0("output_f/aml_iqr_list_", structure, '.txt'))
}
