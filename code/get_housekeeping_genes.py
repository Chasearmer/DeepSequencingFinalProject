import pandas as pd

# Read in the data
housekeepers = pd.read_csv('datasets/Housekeeping_GenesHuman.csv', sep=';')
all_genes = pd.read_csv('datasets/de_Mettl3_KO_vs_WT_m6A.csv')

# print the column names
print(housekeepers.columns)
print(all_genes.columns)

# print shape of the dataframes
print(housekeepers.shape)
print(all_genes.shape)

# Convert both columns to the same case (e.g., uppercase)
housekeepers['Gene.name'] = housekeepers['Gene.name'].str.upper()
# all_genes['gene_id'] = all_genes['gene_id'].str.split('|').str[-1]
all_genes['gene_id'] = all_genes['gene_id'].str.upper()

# split all genes into two dataframes: housekeepers and non-housekeepers
housekeeper_genes = housekeepers['Gene.name'].tolist()
housekeepers = all_genes[all_genes['gene_id'].isin(housekeeper_genes)]
non_housekeepers = all_genes[~all_genes['gene_id'].isin(housekeeper_genes)]

# print the shapes of the dataframes
print(housekeepers.shape)
print(non_housekeepers.shape)

# save the dataframes to csv files
housekeepers.to_csv('datasets/housekeepers.csv', index=False)
non_housekeepers.to_csv('datasets/non_housekeepers.csv', index=False)

