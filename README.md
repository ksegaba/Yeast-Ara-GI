--------------------------------------------------------------------------------
# Folders:
### data/
- Data downloaded June 3, 2021:
    - Raw fitness data: double_mutant_fitness_data_05312024.txt
    - Metadata: fitness_data_notes_updated_05032023.xlsx
- Corrected total seed count (TSC) data (three approaches):
    - brianna_comparemean_tolmer_df_withrelative.csv
    - double_mutant_fitness_data_05312024_corrected_linear.txt
    - double_mutant_fitness_data_05312024_corrected_linear_b.txt
    - double_mutant_fitness_data_05312024_corrected_SpATS.txt

### code/
#### Raw data batch correction
- 0a_raw_data_correction_linear.R (Linear model)
    - For sets (WT, MA, MB, DM) grown on four flats, row/col are random effects and subline ID is a fixed effect; one model per set per flat.
    - For sets grown on one flat, row/col are random effects and subline ID is a fixed effect; one model per set.
- 0b_raw_data_correction_linear_b.R (Linear model)
    - For sets (WT, MA, MB, DM) grown on four flats, row/col/flat are random effects and subline ID is a fixed effect; one model per set.
    - For sets grown on one flat, row/col are random effects and subline ID is a fixed effect; one model per set.
- 0c_raw_data_correction_SpATS.R (Spatial analysis)
    - For sets (WT, MA, MB, DM) grown on four flats, row/col are random effects and subline ID is a fixed effect; one model per set per flat.
    - For sets grown on one flat, row/col are random effects and subline ID is a fixed effect; one model per set.
- 0d_compare_raw_data_correction.ipynb (Results comparison/agreement)

#### XGBoost regression model
- 1a_make_regression_feature_tables.ipynb
- 1b_xgb_regression.py

--------------------------------------------------------------------------------