
import subprocess
import os
import sys
import glob
from os import mkdir
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import tkinter as tk
from tkinter import filedialog as fd, filedialog
import openpyxl
from matplotlib.ticker import ScalarFormatter
from pathlib import Path
print("NU BÖRJAR DET COOLA WOHO!")
# ======================================================================================================================
#chose folder with created readcount, no other files in this specefic folder!!!!
# =============================================================================w=========================================
BASE = Path(__file__).resolve().parent
ROOT = BASE.parent   # folder where main.py lives
folder = BASE/"RESULTS"
#files = filedialog.askdirectory(title= "select folder which has the generated READCOUNTS")#
files = BASE/"READCOUNT"

# ======================================================================================================================
#Import all text-files in folder. Sort and select data from readcount file. Keeping File-ID for alla CpG-markers
# =============================================================================w=========================================
all_rows = []
print(f"{files}")
for file in files.iterdir():
    if file.name.startswith("readcount_"):
        with open(file, "r") as fh:
            for line in fh:
                line=line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t") #alla positioner som standard
                gene = parts[0]
                position = int(parts[1])
                ref_base = parts[2]
                total_depth = int(parts[3])
                stats_block = parts[4:]
                row = {"gene": gene,
                       "position": position,
                       "ref_base": ref_base,
                       "total_depth" : total_depth,
                       "source" : file.name}
                #extract read counts for A C T G N
                for block in stats_block:
                    if block.startswith(("A:","C:","G:","T:","N:")):
                        fields = block.split(":")
                        base = fields[0] #A C G T, definition.
                        read_count = int(fields[1]) #first integer after base
                        row[f"{base}_reads"] = read_count
                all_rows.append(row)
df = pd.DataFrame(all_rows)
print(df.columns.tolist()) # den printar ut 200 coloumns kan säkert fixa med parts.
print(f"exported data in readcounts to dataframe from {folder}")
# ======================================================================================================================
#calculation of Methylation percentage for each markers and sample. and specific id.
#=============================================================================w=========================================
df["%_methylation"] = 100*(df["C_reads"]/(df["C_reads"]+df["T_reads"]))
df["sample_gene_id"] = df["source"].astype(str) + "_" + df["gene"].astype(str)
cols_convert = ["%_methylation","position","total_depth"]
df[cols_convert] = df[cols_convert].apply(pd.to_numeric, errors = 'coerce')
# ======================================================================================================================
# TOTAL NUMBER OF READ PER SAMPLE
#=======================================================================================================================
summary_avg_read_depth = (
    df.groupby("source").agg(
    read_depth = ("total_depth", "mean"),
    Total_Coverage =("total_depth", "sum")))
summary_avg_read_depth = summary_avg_read_depth.rename(columns={"read_depth":"avg_read_depth"})
summary_avg_read_depth = summary_avg_read_depth.reset_index("source")
print(summary_avg_read_depth)
# ======================================================================================================================
# plot 1: total number of reads
#=======================================================================================================================
fig, ax = plt.subplots(figsize=(14, 6))
ax.barh(summary_avg_read_depth["source"], summary_avg_read_depth['Total_Coverage'], color='seagreen')
ax.set_xlabel('Total Coverage (Number of Reads)', fontsize=12)
ax.set_ylabel('sample', fontsize=12)
ax.set_title(f'Total Number of Reads per Sample - {files}-Panel', fontsize=14, fontweight='bold')
ax.invert_yaxis()
ax.grid(axis='x', alpha=0.3)
# Add value labels on bars
for i, v in enumerate(summary_avg_read_depth['Total_Coverage']):
    ax.text(v + 1000, i, str(int(v)), va='center', fontsize=10)
plt.tight_layout()
plt.savefig(BASE/"RESULTS/total_number_of_reads.png", dpi=300, bbox_inches="tight")
plt.close()
# ======================================================================================================================
df.to_excel(BASE/"Save_values_age_anlysis.xlsx",index = False) # fixing of coverage values
#=======================================================================================================================
#summary of coverage and removal of any under 1000
summary_avg_read_depth = summary_avg_read_depth.set_index("source")
df["coverage"] = df["source"].map(summary_avg_read_depth["Total_Coverage"])
removed_samples = df.loc[df["coverage"] <= 1000, "source"].unique()
df = df[df["coverage"] > 1000]
print(df["coverage"].describe())
#print(df.to_string()) #visar hela df.
marker_col = ["gene","source","total_depth"]
marker_coverage_df = df[marker_col].copy()
print(marker_coverage_df)
marker_median = (marker_coverage_df.groupby("gene")["total_depth"].median().sort_values())
marker_order = marker_median.index.tolist()
print(marker_order)
print(marker_median)
# ======================================================================================================================
#Plot 2. Hollow white circles are noticed as outliers!
#=======================================================================================================================
fig,ax = plt.subplots(figsize=(14,6))
sns.stripplot(x = 'gene', y = 'total_depth', data=marker_coverage_df,alpha=0.6,size=6,jitter= True, ax=ax,order= marker_order)
sns.boxplot(x ='gene',y='total_depth',data=marker_coverage_df,width = 0.5,showcaps=False, boxprops=dict(alpha=0.3),medianprops=dict(color='red',linewidth = 2),
            ax=ax, order= marker_order)
y_offset = 0.05*marker_coverage_df["total_depth"].max()
for i,gene in enumerate(marker_order):
    median_val = marker_median[gene]
    ax.text(i,median_val + y_offset,f'{median_val:.0f}',ha='center',va = 'bottom', fontsize = 9,fontweight = 'bold' )
plt.xticks(rotation=90)
ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=10))
ax.yaxis.set_major_formatter(ScalarFormatter())
ax.get_yaxis().get_major_formatter().set_scientific(False)
ax.get_yaxis().get_major_formatter().set_useOffset(False) # Ensure full figure coverage
ax.set_ylim(0, marker_coverage_df["total_depth"].max() * 1.1)
plt.title('Depth per Marker Across All Samples')
plt.xlabel('Marker')
plt.ylabel('Total Read Depth')
plt.tight_layout()
plt.savefig(BASE/"RESULTS/boxplot_of_each_gene_reads.png", dpi=300, bbox_inches="tight")
plt.close()
# ======================================================================================================================
#Calculate the age! THE COOL PART!
#=======================================================================================================================
print(df.head())
print(df.shape)
df["gene_pos"] = df["gene"].astype(str) + "_" + df["position"].astype(str)
methyl_matrix = df.pivot_table(
    index="source",
    columns="gene_pos",
    values="%_methylation")
print(methyl_matrix)
for source in methyl_matrix.index:
    missing_genes = methyl_matrix.columns[methyl_matrix.loc[source].isna()].tolist()
    if missing_genes:
        print(f"{source} is missing markers: {missing_genes}")
markers_coef = {

    #
    #för andra CpG sites i bed-filer, sperma:
    #

    #
    #Coeffecienter i 2021 Pisarek artiklen
    #
 
    #
    #Modiferad bed-file complex visage blood
    #

    #
    #Blood complex visage
    #

    #
    #Next model here
    #
}
AGE_calc = {}
for source in methyl_matrix.index:
    #unmark model you want to use. and # model you do not want to use
    #age = 22.72 #random test för semen visage, changed intercept.
    age = 32.7211426535856 #age semen visage model.
    #age = 28.0 # testar lite random bara!
    #age = 42.2668332790967 #age blood complex visage model
    for gene_pos in methyl_matrix.columns:
        methyl = methyl_matrix.loc[source, gene_pos]
        if pd.isna(methyl):
            continue
        if gene_pos not in markers_coef: #ändra till Gene_pos om indexinerings filen är korrekt!
         continue
        if gene_pos == "ELOVL2_435":
                print(gene_pos)
                methyl = methyl**2
        coef = markers_coef[gene_pos] #sätt gene_pos när indexering är korrekt.
        age= age +  methyl*coef
    AGE_calc[source] = age
methyl_matrix["calculated_age"] = pd.Series(AGE_calc)
print(AGE_calc)
methyl_matrix.to_excel(BASE/"methylation_for_each_sample:gene.xlsx",index = True) # fixing of coverage values
