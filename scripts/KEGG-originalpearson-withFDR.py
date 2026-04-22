##############################################################################
#
#      KEGG pathways Pearson correlation with FDR   ##########################
#
##############################################################################

##############################################################################
# Author: Ioanna Pagani (ioanna.pagani@gmail.com)                            #
# GitHub: promicrobial[](https://github.com/promicrobial)                    #
# Date: 04-10-25                                                             #
# License: MIT                                                               #
# Version: 1.0                                                               #
#                                                                            #
# Description: This script analyzes correlations between KEGG metabolic      #
# pathways in microbiome data from Endpoint control and biofortified arms.   #
# It calculates Pearson correlation coefficients between pathways, performs  #
# Benjamini-Hochberg FDR correction for multiple testing, filters significant#
# correlations (FDR-adjusted p-value < 0.05), and generates two heatmaps:    #
# one for the Pearson correlation coefficients and one for the FDR-adjusted  #
# p-values. Results are saved as CSV files and high-resolution PNG images.   #
#                                                                            #
# Dependencies:                                                              #
# - pandas, numpy, seaborn, matplotlib, scipy, statsmodels                   #
#                                                                            #
# Input:                                                                     #
# - 'Endpoint-control-biofortified-filteredfor70.csv' (contains columns:     #
#   Samples, Pathways, Reads)                                                #
#                                                                            #
# Output:                                                                    #
# - KEGG_Endpoint-control-biofortified_originalpearson_withFDR.csv           #
# - KEGG_Endpoint-control-biofortified_originalpearson_withFDR_filteredfor005.csv #
# - KEGG-Endpoint-control-biofortified_originalpearson_Correlation_Heatmap.png #
# - KEGG-Endpoint-control-biofortified_originalpearson_FDRpvalues_Heatmap.png #
#                                                                            #
# Usage: python KEGG-originalpearson-withFDR.py                              #
#                                                                            #
# Last updated: 04-14-26


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests

# Read the CSV file
df = pd.read_csv('Endpoint-control-biofortified-filteredfor70.csv')

# Pivot the table to create a matrix of metabolic pathways
pivot_df = df.pivot(index='Samples', columns='Pathways', values='Reads')

# Function to calculate the correlation matrix with p-values
def calculate_corr_with_pvalues(df):
    corr_matrix = pd.DataFrame(index=df.columns, columns=df.columns, dtype=float)
    pvals = pd.DataFrame(index=df.columns, columns=df.columns, dtype=float)

    for col in df.columns:
        for row in df.columns:
            if col != row:
                x = df[col].dropna()
                y = df[row].dropna()
                x, y = x.align(y, join='inner')
                if len(x) > 1 and len(y) > 1 and x.var() > 0 and y.var() > 0:
                    corr, pval = pearsonr(x, y)
                    corr_matrix.at[row, col] = corr
                    pvals.at[row, col] = pval
                else:
                    corr_matrix.at[row, col] = np.nan
                    pvals.at[row, col] = np.nan
            else:
                corr_matrix.at[row, col] = 1.0  # Set self-correlation to 1.0
                pvals.at[row, col] = np.nan

    return corr_matrix, pvals

# Calculate the correlation matrix and p-values
corr_matrix, pvals = calculate_corr_with_pvalues(pivot_df)


# Create a DataFrame for correlation matrix with p-values
corr_results = []
for col in corr_matrix.columns:
    for row in corr_matrix.index:
        if col != row and not np.isnan(corr_matrix.at[row, col]):
            corr_results.append([row, col, corr_matrix.at[row, col], pvals.at[row, col]])

corr_df = pd.DataFrame(corr_results, columns=['source', 'target', 'pearson_coefficient', 'p_value'])

# Step to perform Benjamini-Hochberg FDR correction
pvals_flat = corr_df['p_value'].values
fdr_corrected_pvals = multipletests(pvals_flat, method='fdr_bh')[1]  # BH FDR correction

# Add the FDR corrected p-values to the DataFrame
corr_df['BH_FDR_adjusted_p_value'] = fdr_corrected_pvals

# Save the correlation matrix with p-values to a CSV file
corr_df.to_csv('KEGG_Endpoint-control-biofortified_originalpearson_withFDR.csv', index=False)

# Step to filter results where BH-FDR adjusted p-values < 0.05
filtered_corr_df = corr_df[corr_df['BH_FDR_adjusted_p_value'] < 0.05]


# Sort the DataFrame by Pearson correlation coefficient in descending order before pivoting
sorted_filtered_corr_df = filtered_corr_df.sort_values(by='pearson_coefficient', ascending=False)

# Create a pivot table for the sorted, filtered correlation coefficients
filtered_corr_pivot = sorted_filtered_corr_df.pivot(index='source', columns='target', values='pearson_coefficient')

# Save the filtered correlation matrix with p-values to a CSV file
filtered_corr_df.to_csv('KEGG_Endpoint-control-biofortified_originalpearson_withFDR_filteredfor005.csv', index=False)

# Plot Pearson correlation heatmap with corrected tick labels
plt.figure(figsize=(14, 12))
sns.heatmap(filtered_corr_pivot, annot=False, cmap='coolwarm', center=0, cbar=True,
            xticklabels=filtered_corr_pivot.columns, yticklabels=filtered_corr_pivot.index)
plt.title('Endpoint control and biofortified arms: KEGG Pathway Correlation Heatmap (original Pearson)')
plt.xticks(fontsize=6)
plt.yticks(fontsize=6)
plt.tight_layout()
plt.savefig('KEGG-Endpoint-control-biofortified_originalpearson_Correlation_Heatmap.png')  
plt.show()

# Reshape the FDR-corrected p-values into a matrix
pvals_matrix = pd.DataFrame(index=corr_matrix.index, columns=corr_matrix.columns, dtype=float)
for i, row in corr_df.iterrows():
    pvals_matrix.at[row['source'], row['target']] = row['BH_FDR_adjusted_p_value']

# Plot FDR adjusted p-value heatmap
plt.figure(figsize=(14, 12))
sns.heatmap(pvals_matrix.astype(float), annot=False, cmap='viridis',
            xticklabels=pvals_matrix.columns, yticklabels=pvals_matrix.index, vmin=0, vmax=0.05)
plt.title('Endpoint control and biofortified arms: KEGG original Pearson correlation FDR Adjusted P-values Heatmap')
plt.xticks(fontsize=6)
plt.yticks(fontsize=6)
plt.tight_layout()
plt.savefig('KEGG-Endpoint-control-biofortified_originalpearson_FDRpvalues_Heatmap.png')
plt.show()
