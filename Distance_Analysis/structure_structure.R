library(spatstat)
library(png)
library(stringr)

coefficients_list_1 <- list()
coefficients_list_2 <- list()
sample_list <- list()
structure_list_1 <- list()
structure_list_2 <- list()
distance_median_list <- list()
distance_sd_list <- list()
pvalue_list <- list()
iqr_list <- list()

structures <- c('adipocyte', 'art', 'bone', 'macrophage', 'sinusoid', 'stroma')

filenames <- c('H10.csv','H14.csv','H26.csv', 'H27.csv', 'H32.csv', 'H33.csv', 'H35.csv', 
                'H36.csv', 'H37.csv', 'H38.csv', 'H39.csv', 'H41.csv')
# extract the sample name from each filename
sample_names <- c('H10', 'H14', 'H26', 'H27', 'H32', 'H33', 'H35', 'H36', 'H37', 'H38', 'H39', 'H41')
sample_ids <- c("SB67_NBM27_H10_CODEX_Mesmer", "SB67_NBM28_H14_CODEX_Mesmer", 'SB67_NBM36_H26_CODEX_Mesmer',
                'SB67_NBM41_H27_CODEX_Mesmer', 'SB67_NBM31_H32_CODEX_Mesmer', 'SB67_NBM38_H33_CODEX_Mesmer',
                'SB67_NBM37_H35_CODEX_Mesmer', 'SB67_NBM33_H36_CODEX_Mesmer', 'SB67_NBM32_H37_CODEX_Mesmer',
                'SB67_NBM34_H38_CODEX_Mesmer', 'SB67_NBM40_H39_CODEX_Mesmer', 'SB67_NBM39_H41_CODEX_Mesmer')
n_permutations <- 100
for (k in 1:length(sample_names)) {
    sample_name <- sample_names[k]
    print(sample_name)
    filename <- filenames[k]
    sample_data <- read.csv(paste0('h_contour/random/', filename))[, c('x0', 'y0','structure')]

    for (i in 1:length(structures)) {
        structure <- structures[i]
        print(structure)

        cts <- read.csv(paste0('h_contour/',structure,'_csv/', filename))
        
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
        sample_data_crop <- sample_data[sample_data$x0 > xmin & sample_data$x0 < xmax & sample_data$y0 > ymin & sample_data$y0 < ymax, ]
        sample_data_crop_xy <- sample_data_crop[, c('x0', 'y0')]
        for (j in 1:length(structures)) {
            structure_2 <- structures[j]

            coords_crop <- subset(sample_data_crop, structure == structure_2)

            x <- coords_crop$x0
            y <- coords_crop$y0
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
            structure_list_1 <- append(structure_list_1, structure)
            structure_list_2 <- append(structure_list_2, structure_2)
            distance_median_list <- append(distance_median_list, distance_median)
            distance_sd_list <- append(distance_sd_list, distance_sd)
            pvalue_list <- append(pvalue_list, sum(p_distance_median < distance_median) / n_permutations)
            iqr_list <- append(iqr_list, distance_iqr)
        }
    }
}

coefficients_list_1_vector <- as.character(coefficients_list_1)
coefficients_list_2_vector <- as.character(coefficients_list_2)
sample_list_vector <- as.character(sample_list)
structure_list_1_vector <- as.character(structure_list_1)
structure_list_2_vector <- as.character(structure_list_2)
distance_median_list_vector <- as.character(distance_median_list)
distance_sd_list_vector <- as.character(distance_sd_list)
pvalue_list_vector <- as.character(pvalue_list)
iqr_list_vector <- as.character(iqr_list)

writeLines(coefficients_list_1_vector,"output/st_coefficients_list_1.txt") 
writeLines(coefficients_list_2_vector,"output/st_coefficients_list_2.txt")
writeLines(sample_list_vector,"output/st_sample_list.txt")
writeLines(structure_list_1_vector,"output/st_structure_list_1.txt")
writeLines(structure_list_2_vector,"output/st_structure_list_2.txt")
writeLines(distance_median_list_vector,"output/st_distance_median_list.txt")
writeLines(distance_sd_list_vector,"output/st_distance_sd_list.txt")
writeLines(pvalue_list_vector,"output/st_pvalue_list.txt")
writeLines(iqr_list_vector,"output/st_iqr_list.txt")

