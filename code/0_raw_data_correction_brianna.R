#!/usr/bin/env Rscript
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=0_raw_data_correction_brianna
#SBATCH --output=../logs/%x-%j.out

system('module purge; module load GCC/11.2.0  OpenMPI/4.1.1  R/4.3.1')

# DESCRIPTION: Correction of raw fitness data for single and double mutants
# using Brianna's linear model to estimate the marginalized means of total seed count
# June 4, 2024

# R v 4.3.1
library(emmeans)

# Load the data
df <- read.csv('../data/double_mutant_fitness_data_05312024.txt', sep='\t', header=T)

# Collect fitted values for TSC
res <- data.frame()

# Apply spatial analysis per set
for (set in unique(df$Set)) {
    set_df <- df[df$Set == set,] # set data

    # Remove rows with missing values in TSC (total seed count)
    set_df <- set_df[!is.na(set_df$TSC),]

    # Set WT as the reference level
    set_df$Genotype <- factor(set_df$Genotype, levels=c('WT', 'MA', 'MB', 'DM'))

    # For sets grown on multiple flats, include flat as random effect
    if (length(unique(set_df$Flat)) > 1) {
        # Fit a linear model
        model <- lme4::lmer(TSC ~ Genotype + (1|Column) + (1|Row) + (1|Flat), data = set_df)

        # Save the model
        saveRDS(model, paste0('../output/0_raw_data_correction/Set_', set, '_brianna_model.rds'))

        # Calculate estimated means
        emm <- as.data.frame(emmeans(model, ~ Genotype))
        out <- dplyr::left_join(set_df, emm, by='Genotype')
        res <- rbind(res, out)
        }
    if (length(unique(set_df$Flat)) == 1) {
        # Fit a linear model
        model <- lme4::lmer(TSC ~ Genotype + (1|Column) + (1|Row), data = set_df)

        # Save the model
        saveRDS(model, paste0('../output/0_raw_data_correction/Set_', set, '_brianna_model.rds'))

        # Calculate estimated means
        emm <- as.data.frame(emmeans(model, ~ Genotype))
        out <- dplyr::left_join(set_df, emm, by='Genotype')
        res <- rbind(res, out)
    }
}

# Save the corrected data
write.csv(res, '../data/double_mutant_fitness_data_05312024_TSC_corrected_brianna.txt', row.names=F, quote=F)