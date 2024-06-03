# Spatial correction of raw fitness data for single and double mutants
# June 3, 2024

install.packages('SpATS')
library(SpATS)

df <- read.csv('../data/double_mutant_fitness_data_05312024.txt', sep='\t', header=T)

# yield data is centered and scaled


[11:09] Izquierdo Romero, Paulo
MI18$col_f = factor(MI18$P)
MI18$row_f = factor(MI18$R)
 
fit.NE18 = SpATS(response = "yield_kg_ha", 
    spatial = ~ SAP (P, R, nseg = c(10,10)),
    genotype = "ID",
    fixed = NULL,
    genotype.as.random = False,
    random = ~row_f + col_f,
    data= NE18,
    control = list(tolerance = 1e-03))