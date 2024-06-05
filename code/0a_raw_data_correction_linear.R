#!/usr/bin/env Rscript
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=0_raw_data_correction_linear
#SBATCH --output=../logs/%x-%j.out

system('module purge; module load GCC/11.2.0  OpenMPI/4.1.1  R/4.3.1')

# DESCRIPTION: Correction of raw fitness data for single and double mutants
# using a linear model to estimate the marginalized means of total seed count
# June 3, 2024

# R v 4.3.1
# install.packages('emmeans')
# install.packages('pbkrtest')
# install.packages('lmerTest')
library(emmeans)

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

            # Set WT as the reference level
            flat_df$Genotype <- factor(flat_df$Genotype, levels=c('WT', 'MA', 'MB', 'DM'))

            # Fit a linear model
            model <- lme4::lmer(TSC ~ Genotype + (1|Column) + (1|Row), data = flat_df)
            
            # Calculate estimated means
            emm <- as.data.frame(emmeans(model, ~ Genotype))

            # Save the model
            saveRDS(model, paste0('../output/0_raw_data_correction/Set_', set, '_flat_', flat, '_linear_model.rds'))

            # Collect fitted values
            out <- dplyr::left_join(flat_df, emm, by='Genotype')
            res <- rbind(res, out)
        }
    } 
    if (length(unique(set_df$Flat)) == 1) {
        # Set WT as the reference level
        set_df$Genotype <- factor(set_df$Genotype, levels=c('WT', 'MA', 'MB', 'DM'))
        
        # Fit a linear model
        model <- lme4::lmer(TSC ~ Genotype + (1|Column) + (1|Row), data = set_df)
        
        # Calculate estimated means
        emm <- as.data.frame(emmeans(model, ~ Genotype))

        # Save the model
        saveRDS(model, paste0('../output/0_raw_data_correction/Set_', set, '_linear_model.rds'))

        # Collect fitted values
        out <- dplyr::left_join(set_df, emm, by='Genotype')
        res <- rbind(res, out)
    }
}

# Save the corrected data
write.csv(res, '../data/double_mutant_fitness_data_05312024_TSC_corrected_linear.txt', row.names=F, quote=F)