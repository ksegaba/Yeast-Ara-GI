--------------------------------------------------------------------------------
# Folders:
### data/
- Data downloaded June 3, 2021:
    - Raw fitness data: double_mutant_fitness_data_05312024.txt
    - Metadata: fitness_data_notes_updated_05032023.xlsx
- Corrected total seed count (TSC) data (three approaches):
    - brianna_comparemean_tolmer_df_withrelative.csv
    - double_mutant_fitness_data_05312024_corrected_emmeans.txt
    - double_mutant_fitness_data_05312024_corrected_SpATS.txt

### code/
- 0_raw_data_correction_brianna.R (Linear model)
    - For sets (WT, MA, MB, DM) grown on four flats, row/col/flat are random effects and subline ID is a fixed effect; one model per set.
    - For sets grown on one flat, row/col are random effects and subline ID is a fixed effect; one model per set.
- 0_raw_data_correction_linear.R (Linear model)
    - For sets (WT, MA, MB, DM) grown on four flats, row/col are random effects and subline ID is a fixed effect; one model per set per flat.
    - For sets grown on one flat, row/col are random effects and subline ID is a fixed effect; one model per set.
- 0_raw_data_correction_SpATS.R (Spatial analysis)
    - For sets (WT, MA, MB, DM) grown on four flats, row/col are random effects and subline ID is a fixed effect; one model per set per flat.
    - For sets grown on one flat, row/col are random effects and subline ID is a fixed effect; one model per set.
- 1_compare_raw_data_correction.ipynb (Results comparison/agreement)
--------------------------------------------------------------------------------