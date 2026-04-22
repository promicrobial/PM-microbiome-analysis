##############################################################################
# Author: Ioanna Pagani (ioanna.pagani@gmail.com)                            #
# GitHub: promicrobial[](https://github.com/promicrobial)                    #
# Date: 04-10-25                                                             #
# License: MIT                                                               #
# Version: 1.0                                                               #
#                                                                            #
# Description: This script performs differential network analysis between    #
# Baseline and Endpoint KEGG pathway correlation networks. It identifies     #
# unique edges, common edges with significant correlation changes (>0.2),    #
# and computes differential correlation coefficients. The script also        #
# calculates weighted degree centrality differences, generates a hierarchical#
# clustermap of differential correlations, and prepares nodes and edges      #
# files compatible with Gephi for network visualization.                     #
#                                                                            #
# Dependencies:                                                              #
# - pandas, numpy, seaborn, matplotlib, scipy, networkx, statsmodels         #
#                                                                            #
# Input:                                                                     #
# - KEGG_Baseline-control-biofortified_originalpearson_withFDR_filteredfor005.csv #
# - KEGG_Endpoint-control-biofortified_originalpearson_withFDR_filteredfor005.csv #
#                                                                            #
# Output:                                                                    #
# - differential_edges_endpoint-baseline.csv                                 #
# - network_stats_endpoint-baseline.csv                                      #
# - differential_clustermap_endpoint-baseline.png                            #
# - gephi_nodes_endpoint-baseline.csv                                        #
# - gephi_edges_endpoint-baseline.csv                                        #
#                                                                            #
# Usage: python pearson-differential-networks-github.py                      #
#                                                                            #
# Last updated: 04-14-26                                                     #
##############################################################################


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from collections import Counter
import networkx as nx
from statsmodels.stats.multitest import multipletests

def differential_network_analysis(file1, file2, output_diff_edges, output_network_stats, clustermap_file):
    # Read both files
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    # Rename columns to ensure consistency
    df1 = df1.rename(columns={"pearson": "pearson_coefficient", "p-value": "p_value", "BH_FDR_adjusted_p_value": "fdr_p_value"})
    df2 = df2.rename(columns={"pearson": "pearson_coefficient", "p-value": "p_value", "BH_FDR_adjusted_p_value": "fdr_p_value"})

    # Filter for significant correlations (FDR < 0.05)
    df1_sig = df1[df1["fdr_p_value"] < 0.05].copy()
    df2_sig = df2[df2["fdr_p_value"] < 0.05].copy()

    # Create network graphs
    G_Baseline = nx.Graph()
    G_Endpoint = nx.Graph()

    # Add edges to Baseline network
    for _, row in df1_sig.iterrows():
        G_Baseline.add_edge(row["source"], row["target"], weight=row["pearson_coefficient"])

    # Add edges to Endpoint network
    for _, row in df2_sig.iterrows():
        G_Endpoint.add_edge(row["source"], row["target"], weight=row["pearson_coefficient"])

    # Find differential edges
    Baseline_edges = set(G_Baseline.edges())
    Endpoint_edges = set(G_Endpoint.edges())
    unique_Baseline = Baseline_edges - Endpoint_edges
    unique_Endpoint = Endpoint_edges - Baseline_edges
    common_edges = Baseline_edges.intersection(Endpoint_edges)

    # Compute differences for common edges
    diff_edges = []
    for edge in common_edges:
        Baseline_weight = G_Baseline[edge[0]][edge[1]]['weight']
        Endpoint_weight = G_Endpoint[edge[0]][edge[1]]['weight']
        diff = Endpoint_weight - Baseline_weight
        if abs(diff) > 0.2:  # Threshold for significant difference
            diff_edges.append([edge[0], edge[1], Baseline_weight, Endpoint_weight, diff])

    # Add unique edges
    for edge in unique_Baseline:
        Baseline_weight = G_Baseline[edge[0]][edge[1]]['weight']
        diff_edges.append([edge[0], edge[1], Baseline_weight, 0.0, -Baseline_weight])
    for edge in unique_Endpoint:
        Endpoint_weight = G_Endpoint[edge[0]][edge[1]]['weight']
        diff_edges.append([edge[0], edge[1], 0.0, Endpoint_weight, Endpoint_weight])

    # Save differential edges
    diff_df = pd.DataFrame(diff_edges, columns=["source", "target", "Baseline_coefficient",
                                                "Endpoint_coefficient", "coefficient_difference"])
    diff_df.to_csv(output_diff_edges, index=False)
    print(f"Differential edges saved to {output_diff_edges}")

    # Compute network statistics (degree centrality)
    Baseline_degree = dict(G_Baseline.degree(weight='weight'))
    Endpoint_degree = dict(G_Endpoint.degree(weight='weight'))
    all_pathways = sorted(set(Baseline_degree.keys()).union(Endpoint_degree.keys()))

    # Prepare network stats DataFrame
    stats_data = []
    for pathway in all_pathways:
        Baseline_deg = Baseline_degree.get(pathway, 0)
        Endpoint_deg = Endpoint_degree.get(pathway, 0)
        stats_data.append([pathway, Baseline_deg, Endpoint_deg, Endpoint_deg - Baseline_deg])

    stats_df = pd.DataFrame(stats_data, columns=["Pathway", "Baseline Degree Centrality",
                                                 "Endpoint Degree Centrality", "Degree Difference"])
    stats_df = stats_df.sort_values(by="Degree Difference", key=abs, ascending=False)
    stats_df.to_csv(output_network_stats, index=False)
    print(f"Network statistics saved to {output_network_stats}")

    # Create pivot table for clustermap (differential correlations)
    pivot_df = pd.DataFrame(0.0, index=all_pathways, columns=all_pathways)
    for _, row in diff_df.iterrows():
        pivot_df.loc[row["source"], row["target"]] = row["coefficient_difference"]
        pivot_df.loc[row["target"], row["source"]] = row["coefficient_difference"]

    # Generate hierarchical clustermap
    plt.figure(figsize=(18, 16))
    sns.clustermap(pivot_df, row_cluster=True, col_cluster=True, method='ward', cmap="coolwarm",
                   vmin=-1, vmax=1, figsize=(18, 16), xticklabels=True, yticklabels=True)
    plt.title("Differential pearson Correlations (Endpoint vs. Baseline)")
    plt.savefig(clustermap_file, bbox_inches='tight')
    plt.close()
    print(f"Clustermap saved to {clustermap_file}")


# Example usage
differential_network_analysis(
    "KEGG_Baseline-control-biofortified_originalpearson_withFDR_filteredfor005.csv",
    "KEGG_Endpoint-control-biofortified_originalpearson_withFDR_filteredfor005.csv",
    "differential_edges_endpoint-baseline.csv",
    "network_stats_endpoint-baseline.csv",
    "differential_clustermap_endpoint-baseline.png"
)


def prepare_gephi_files(edges_file, stats_file, nodes_output, edges_output):
    # Read the input files
    edges_df = pd.read_csv(edges_file)
    stats_df = pd.read_csv(stats_file)

    # --- Prepare Nodes File ---
    # Create a nodes DataFrame with unique pathways
    pathways = sorted(set(edges_df["source"]).union(edges_df["target"]))
    nodes_df = pd.DataFrame({"Id": pathways, "Label": pathways})

    # Merge with stats_df to add node attributes
    nodes_df = nodes_df.merge(stats_df, left_on="Id", right_on="Pathway", how="left")
    nodes_df = nodes_df[["Id", "Label", "Baseline Degree Centrality", "Endpoint Degree Centrality", "Degree Difference"]]

    # Fill missing values (in case some pathways in edges are not in stats)
    nodes_df.fillna({"Baseline Degree Centrality": 0, "Endpoint Degree Centrality": 0, "Degree Difference": 0}, inplace=True)

    # Save nodes file
    nodes_df.to_csv(nodes_output, index=False)
    print(f"Nodes file saved to {nodes_output}")

    # --- Prepare Edges File ---
    # Create edges DataFrame for Gephi
    edges_df_gephi = edges_df[["source", "target", "coefficient_difference", "Baseline_coefficient", "Endpoint_coefficient"]].copy()
    edges_df_gephi.columns = ["Source", "Target", "Weight", "Baseline_Coefficient", "Endpoint_Coefficient"]

    # Add edge type (undirected for correlation networks)
    edges_df_gephi["Type"] = "Undirected"

    # Add edge category based on coefficient_difference
    edges_df_gephi["Edge_Category"] = edges_df_gephi["Weight"].apply(
        lambda x: "Positive Difference" if x > 0.2 else "Negative Difference" if x < -0.2 else "Small Difference"
    )

    # Save edges file
    edges_df_gephi.to_csv(edges_output, index=False)
    print(f"Edges file saved to {edges_output}")

# Example usage
prepare_gephi_files(
    edges_file="differential_edges_endpoint-baseline.csv",
    stats_file="network_stats_endpoint-baseline.csv",
    nodes_output="gephi_nodes_endpoint-baseline.csv",
    edges_output="gephi_edges_endpoint-baseline.csv"
)
