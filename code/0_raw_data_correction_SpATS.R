#!/usr/bin/env Rscript
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=0_raw_data_correction_SpATS
#SBATCH --output=../logs/%x-%j.out

system('module purge; module load GCC/9.3.0  OpenMPI/4.0.3  R/4.0.3')

# DESCRIPTION: Spatial correction of raw fitness data for single and double mutants
# June 3, 2024

# R v 4.0.3
# install.packages('SpATS')
library(SpATS)

# Load the data
df <- read.csv('../data/double_mutant_fitness_data_05312024.txt', sep='\t', header=T)

# Collect fitted values for TSC
res <- data.frame()

# Apply spatial analysis per flat (for the most part, one set per flat)
for (set in unique(df$Set)) {
    set_df <- df[df$Set == set,] # set data

    # Remove rows with missing values in TSC (total seed count)
    set_df <- set_df[!is.na(set_df$TSC),]

    # For sets grown on multiple flats, run the model per flat
    if (length(unique(set_df$Flat)) > 1) {
        for (flat in unique(set_df$Flat)){
            flat_df <- set_df[set_df$Flat == flat,] # flat data

            # Transform fixed and random variables/ set WT as reference level
            flat_df$R <- factor(flat_df$Row)
            flat_df$C <- factor(flat_df$Column)
            flat_df$geno <- factor(flat_df$Genotype, levels=c('WT', 'MA', 'MB', 'DM'))

            # Fit spatial model
            fit.TSC <- SpATS(response = 'TSC',
                spatial = ~ SAP (Column, Row, nseg = c(length(levels(flat_df$C)), length(levels(flat_df$R)))),
                genotype = 'geno',
                fixed = NULL,
                genotype.as.random = FALSE,
                random = ~ C + R,
                data = flat_df,
                control = list(tolerance = 1e-03))
            
            # Save the model
            saveRDS(fit.TSC, paste0('../output/0_raw_data_correction/Set_', set, '_flat_', flat, '_spatial_model.rds'))

            # Plot results
            pdf(paste0('../output/0_raw_data_correction/Set_', set, '_flat_', flat, '_spatial_model.pdf'))
            plot(fit.TSC)
            dev.off()

            # Collect fitted values
            out <- cbind(fit.TSC$data, fit.TSC$fitted)
            res <- rbind(res, out)
        }
    }
    if (length(unique(set_df$Flat)) == 1) {
        # Set WT as the reference level
        set_df$R <- factor(set_df$Row)
        set_df$C <- factor(set_df$Column)
        set_df$geno <- factor(set_df$Genotype, levels=c('WT', 'MA', 'MB', 'DM'))

        fit.TSC <- SpATS(response = 'TSC', 
            spatial = ~ SAP (Column, Row, nseg = c(length(levels(set_df$C)), length(levels(set_df$R)))),
            genotype = 'geno',
            fixed = NULL,
            genotype.as.random = FALSE,
            random = ~ C + R,
            data = set_df,
            control = list(tolerance = 1e-03))

        saveRDS(fit.TSC, paste0('../output/0_raw_data_correction/Set_', set, '_spatial_model.rds'))

        pdf(paste0('../output/0_raw_data_correction/Set_', set, '_spatial_model.pdf'))
        plot(fit.TSC)
        dev.off()

        out <- cbind(fit.TSC$data, fit.TSC$fitted)
        res <- rbind(res, out)
    }
}

# Save the corrected data
write.csv(res, '../data/double_mutant_fitness_data_05312024_TSC_corrected_SpATS.txt', row.names=F, quote=F)