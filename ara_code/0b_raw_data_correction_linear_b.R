#!/usr/bin/env Rscript
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=0_raw_data_correction_brianna
#SBATCH --output=../logs/%x-%j.out

system('module purge; module load GCC/11.2.0  OpenMPI/4.1.1  R/4.3.1')

# DESCRIPTION: Correction of raw fitness data for single and double mutants
# using Brianna's linear model to estimate the marginalized means of total seed count
# June 4, 2024 Note: double_mutant_fitness_data_05312024.txt was corrected with the model label ~ Genotype + (1|Row) + (1|Column) + (1|Flat)
# Oct 13, 2024 Note: fitness_data_for_Kenia_09232024.xlsx was corrected with the model label ~ Genotype + (1|Type) + (1|Flat)

# R v 4.3.1
library(emmeans)
library(dplyr)
library(readxl)

# Load the data
# df <- read.csv('../ara_data/double_mutant_fitness_data_05312024.txt', sep='\t', header=T)
path <- '/home/seguraab/ara-kinase-prediction/data/20240923_melissa_ara_data/fitness_data_for_Kenia_09232024.xlsx'
df <- as.data.frame(read_xlsx(path, sheet='with_border_cells'))
labels <- c('GN', 'PG', 'DTB', 'LN', 'DTF', 'SN', 'WO', 'FN', 'SPF', 'TSC', 'SH')

counter <- 1
for (label in labels){
    print(paste0('LABEL ', label, ' -----------------------------------------'))

    # Collect fitted values for each label
    if (counter == 1) res <- data.frame()
    if (counter != 1) res2 <- data.frame()

    # If the label column is not numeric
    if (!is.numeric(df[,label])) df[,label] <- as.numeric(df[[label]])

    # Run linear regression per set
    for (set in unique(df$Set)) {
        print(paste0('Set ', set))
        set_df <- df[df$Set == set,] # set data

        # Remove rows with missing values in label
        set_df <- set_df[!is.na(set_df[,label]),]

        # For sets grown on multiple flats, include flat as random effect
        if (length(unique(set_df$Flat)) > 1) {
            # Set WT as the reference level
            set_df$Genotype <- factor(set_df$Genotype, levels=c('WT', 'MA', 'MB', 'DM'))

            # Fit a linear model
            model <- lme4::lmer(paste0(label, ' ~ Genotype + (1|Type) + (1|Flat)'), data = set_df)

            # Save the model
            saveRDS(model, paste0('../output/0_raw_data_correction_ara/Set_', set,
                '_', label, '_brianna_model.rds'))

            # Calculate estimated means
            emm <- as.data.frame(emmeans(model, ~ Genotype))
            emm <- rename_with(emm, ~ paste0(label, "_", .))
            out <- left_join(set_df, emm, by=c('Genotype'=paste0(label, '_Genotype')))
            if (counter == 1) res <- rbind(res, out)
            if (counter != 1) res2 <- rbind(res2, out)
        }
        if (length(unique(set_df$Flat)) == 1) {
            # Set WT as the reference level
            set_df$Genotype <- factor(set_df$Genotype, levels=c('WT', 'MA', 'MB', 'DM'))
            tryCatch({
                # Fit a linear model
                model <- lme4::lmer(paste0(label, ' ~ Genotype + (1|Type)'), data = set_df)

                # Save the model
                saveRDS(model, paste0('../output/0_raw_data_correction_ara/Set_',
                    set, '_', label, '_brianna_model.rds'))

                # Calculate estimated means
                emm <- as.data.frame(emmeans(model, ~ Genotype))
                emm <- rename_with(emm, ~ paste0(label, "_", .))
                out <- left_join(set_df, emm, by=c('Genotype'=paste0(label, '_Genotype')))
                if (counter == 1) res <- rbind(res, out)
                if (counter != 1) res2 <- rbind(res2, out)
            }, error = function(e) {
                print(e)
            }, warning = function(w) {
                print(w)
            })
        }
    }

    if (counter != 1) {
        print(dim(res))
        print(dim(res2))
        res <- left_join(res, res2, by=c('Set', 'Flat', 'Column', 'Row', 
            'Number', 'Type', 'Genotype', 'Subline', 'MA', 'MB'), keep=F)
    }
    counter <- counter + 1
}

# Save the corrected data
res <- res[, !grepl("\\.x$|\\.y$", names(res))]
res <- left_join(df, res, by=c('Set', 'Flat', 'Column', 'Row', 'Number', 'Type',
    'Genotype', 'Subline', 'MA', 'MB'))
res <- rename_with(res, ~ sub("\\.x$", "", .), .cols = ends_with(".x"))
# write.table(res, paste0('../ara_data/double_mutant_fitness_data_05312024_all_corrected_brianna.txt'), row.names=F, quote=F, sep='\t')
write.table(res, paste0('../ara_data/fitness_data_for_Kenia_09232024_all_corrected_brianna.txt'), row.names=F, quote=F, sep='\t')