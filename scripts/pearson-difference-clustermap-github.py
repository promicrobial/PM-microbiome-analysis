##############################################################################
# Author: Ioanna Pagani (ioanna.pagani@gmail.com)                            #
# GitHub: promicrobial[](https://github.com/promicrobial)                     #
# Date: 04-10-25                                                             #                                                           #
#                                                                            #
# Description: This script compares Pearson correlation networks between     #
# Baseline and Endpoint arms for KEGG metabolic pathways. It calculates the  #
# difference in Pearson coefficients between the two conditions, generates   #
# a hierarchical clustermap visualizing these differences, and identifies    #
# the top 40 most frequently involved pathways in significant correlation    #
# changes (|difference| > 0.2), including separate counts for positive and   #
# negative differences.                                                      #
#                                                                            #
# Dependencies:                                                              #
# - pandas, seaborn, matplotlib, scipy, numpy, collections                   #
#                                                                            #
# Input:                                                                     #
# - KEGG_Baseline-control-biofortified_originalpearson_withFDR_filteredfor005.csv #
# - KEGG_Endpoint-control-biofortified_originalpearson_withFDR_filteredfor005.csv #
#                                                                            #
# Output:                                                                    #
# - pearson-difference-Baseline-Endpoint-for-both-biofortified-control.csv   #
# - pearson-difference-Clustermap-Baseline-Endpoint-for-both-biofortified-control.png #
# - Top-40-most-frequent-pathways-Baseline-Endpoint-for-both-biofortified-control.csv #
#                                                                            #
# Usage: python pearson-difference-clustermap.py                                      #
#                                                                            #
# Last updated: 04-14-26                                                     #
##############################################################################


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import numpy as np
from collections import Counter

def compare_correlations(file1, file2, output_file, clustermap_file, top_pathways_file):
    # Read both files
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    # Print column names to debug
    print("Columns in file1:", df1.columns.tolist())
    print("Columns in file2:", df2.columns.tolist())

    # Ensure correct column names are used
    df1 = df1.rename(columns={"pearson coefficient": "pearson_coefficient", "p-value": "p_value", "fdr adjusted p-value": "BH_FDR_adjusted_p_value"})
    df2 = df2.rename(columns={"pearson coefficient": "pearson_coefficient", "p-value": "p_value", "fdr adjusted p-value": "BH_FDR_adjusted_p_value"})

    # Merge data on source and target pathway names
    merged_df = df1.merge(df2, on=["source", "target"], suffixes=("_1", "_2"))

    # Compute pearson coefficient difference
    merged_df["pearson_coefficient_difference"] = merged_df["pearson_coefficient_2"] - merged_df["pearson_coefficient_1"]

    # Select required columns and rename them
    output_df = merged_df[["source", "target", "pearson_coefficient_difference", "p_value_2", "BH_FDR_adjusted_p_value_2"]]
    output_df.columns = ["source", "target", "pearson coefficient difference", "p-value", "fdr adjusted p-value"]

    # Save to CSV
    output_df.to_csv(output_file, index=False)
    print(f"Output saved to {output_file}")

    # Create pivot table for clustermap
    pivot_df = output_df.pivot(index="source", columns="target", values="pearson coefficient difference").fillna(0)

    # Cluster pathways based on positive/negative pearson differences
    linkage = sch.linkage(pivot_df, method='ward')
    dendro = sch.dendrogram(linkage, labels=pivot_df.index, no_plot=True)
    ordered_indices = dendro['leaves']
    ordered_labels = [pivot_df.index[i] for i in ordered_indices]
    pivot_df = pivot_df.loc[ordered_labels, ordered_labels]

    # Generate hierarchical clustermap
    plt.figure(figsize=(10, 8))
    sns.clustermap(pivot_df, row_cluster=True, col_cluster=True, method='ward', cmap="coolwarm", vmin=-1, vmax=1, figsize=(12, 10))
    plt.savefig(clustermap_file)
    plt.close()
    print(f"Clustermap saved to {clustermap_file}")

    # Identify top 20 pathways with highest frequency in pearson coefficient differences > 0.2 or < -0.2
    filtered_df = output_df[(output_df["pearson coefficient difference"] > 0.2) | (output_df["pearson coefficient difference"] < -0.2)]
    pathways = list(filtered_df["source"]) + list(filtered_df["target"])
    pathway_counts = Counter(pathways)

    # Compute separate counts for positive and negative differences
    pos_counts = Counter(list(filtered_df[filtered_df["pearson coefficient difference"] > 0.2]["source"]) +
                         list(filtered_df[filtered_df["pearson coefficient difference"] > 0.2]["target"]))
    neg_counts = Counter(list(filtered_df[filtered_df["pearson coefficient difference"] < -0.2]["source"]) +
                         list(filtered_df[filtered_df["pearson coefficient difference"] < -0.2]["target"]))

    top_20_pathways = pathway_counts.most_common(40)
    top_20_data = []
    for pathway, total_freq in top_20_pathways:
        pos_freq = pos_counts.get(pathway, 0)
        neg_freq = neg_counts.get(pathway, 0)
        top_20_data.append([pathway, total_freq, pos_freq, neg_freq])

    # Save top 20 pathways with breakdown to a file
    top_pathways_df = pd.DataFrame(top_20_data, columns=["Pathway", "Total Frequency", "Positive Difference Frequency", "Negative Difference Frequency"])
    top_pathways_df.to_csv(top_pathways_file, index=False)
    print(f"Top 40 pathways saved to {top_pathways_file}")

# Example usage
compare_correlations(
    "KEGG_Baseline-control-biofortified_originalpearson_withFDR_filteredfor005.csv",
    "KEGG_Endpoint-control-biofortified_originalpearson_withFDR_filteredfor005.csv",
    "pearson-difference-Baseline-Endpoint-for-both-biofortified-control.csv",
    "pearson-difference-Clustermap-Baseline-Endpoint-for-both-biofortified-control.png",
    "Top-40-most-frequent-pathways-Baseline-Endpoint-for-both-biofortified-control.csv"
)
