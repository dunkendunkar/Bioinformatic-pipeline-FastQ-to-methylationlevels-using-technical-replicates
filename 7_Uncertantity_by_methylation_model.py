
from tkinter import filedialog as fd, filedialog
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#change files to if folder is changed!
# ======================================================================================================================
#Start by generation of model using existing data: SELECT FOLDER WITH METHYLATION LEVELS, ESTIMATED AGE, CHRONOLOGICAL AGE
# ======================================================================================================================
BASE = Path(__file__).resolve().parent
files =fd.askopenfilename(title= "SELECT FOLDER WITH METHYLATION LEVELS, ESTIMATED AGE, CHRONOLOGICAL AGE", initialdir= "/home/froste/PycharmProjects/PythonProject/")
df_tot = pd.read_excel(files)
print(df_tot.head(20))
age = df_tot["AGE"]
x_age = age
age_pred = df_tot["calculated_age"]
x_pred = age_pred
deviation = (age - age_pred).abs()
dev_not_abs = age - age_pred
df_tot["dev_not_abs"] = dev_not_abs
df_tot["deviation"] = deviation
exclude = ['AGE', 'calculated_age', 'source', 'deviation', "rep", "dev_not_abs"]
print(df_tot.to_string)
marker_cols = df_tot.select_dtypes(include="number").columns.tolist()
marker_cols = [c for c in marker_cols if c not in exclude]
print(marker_cols)
# ======================================================================================================================
#Generate Z-value and D-value for both the age group and the samples itself for fitting of linear model
# ======================================================================================================================
#                               -- Now it is time to create a specific methylation% confidence interval MAE generation ---
#                              -- It is based on the deviance of methylation from specific age brackets.               ---
# 1. Copy data and create age groups
data = df_tot.copy()
data["group"] = pd.cut(
    data["AGE"],
    bins=[10,45,55,np.inf],
    labels=["young","middle aged","above middle age"])
group_mean = {}
group_sd = {}
# 2. Compute group means, SDs, and deviations
for m in marker_cols:
    data[f"{m}_group_mean"] = data.groupby("group")[m].transform("mean")
    data[f"{m}_group_sd"]   = data.groupby("group")[m].transform("std")
    group_mean[m] = data.groupby("group")[m].mean()
    group_sd[m] = data.groupby("group")[m].std()
    # data[f"{m}_group_sd"].replace(0, np.nan, inplace=True)
# 3. Compute Z-scores
for m in marker_cols:
    data[f"{m}_z"] = (
        data[m] - data[f"{m}_group_mean"]
    ) / data[f"{m}_group_sd"]
# 4. Compute D-value (combined deviation score)
z_cols = [f"{m}_z" for m in marker_cols]
data["D_value"] = np.sqrt((data[z_cols]**2).sum(axis=1))

# ======================================================================================================================
#Fitting of linear model using available data Error(age - pred.age).abs = A*D+b
# ======================================================================================================================
from sklearn.linear_model import LinearRegression
data["D_value"] = data["D_value"] / np.sqrt(len(marker_cols)) #normalize error with /n, n = amount if markers
X = data[["D_value"]] #coefficent
y = data["deviation"]
reg = LinearRegression()
reg.fit(X, y)
plt.scatter(X,y,c = data["AGE"], cmap="viridis")
plt.xlabel("D_value (methylation deviation for age group)")
plt.ylabel("Estimated MAE for samples")
plt.title("Estimated MAE by methylation deviation from age group methylation valyes")
plt.colorbar(label="Age")
plt.show()
MAE_sample = reg.predict(X)
print(MAE_sample)
# ======================================================================================================================
#Predict using the actual samples:
# ======================================================================================================================'

df = pd.read_excel(BASE/"meth_calcAge_combinedReplicates.xlsx")
exclude = ["basename", "status_tot_reads", "status_reads_marker","Status_replicates","calculated_age",'AGE', 'calculated_age', 'source', 'deviation', "rep", "dev_not_abs"]
marker_cols = df_tot.select_dtypes(include="number").columns.tolist()
marker_cols = [c for c in marker_cols if c not in exclude]
# ======================================================================================================================
#Use already generated sd.mean for markers and mean.meth values for age bracket.
# ======================================================================================================================
# 1. Copy data and create age groups
data = df.copy()
data["group"] = pd.cut(
    data["calculated_age"],
    bins=[10,45,55,np.inf],
    labels=["young","middle aged","above middle age"])
# 2. Compute group means, SDs, and deviations
# 3. Compute Z-scores
# data = df.copy()  # your sample
sample_group = data["group"].iloc[0]
for m in marker_cols:
    data[f"{m}_z"] = (
        data[m] - group_mean[m][sample_group]
    ) / group_sd[m][sample_group]
z_cols = [f"{m}_z" for m in marker_cols]
data["D_value"] = np.sqrt((data[z_cols]**2).sum(axis=1))
data["estimated_MAE"] = reg.predict(data[["D_value"]])
# 4. Compute D-value (combined deviation score)
z_cols = [f"{m}_z" for m in marker_cols]
data["D_value"] = np.sqrt((data[z_cols]**2).sum(axis=1))
data["estimated_MAE"] = reg.predict(data[["D_value"]])
print(data["estimated_MAE"])
cols_to_drop = [c for c in data.columns if c.endswith("_z")]
data = data.drop(columns=cols_to_drop)
data.to_excel(BASE/"meth_calcAge_combinedReplicates.xlsx",index = True) # fixing of coverage values
# ======================================================================================================================
#Plot for seeing:
# ======================================================================================================================
data["group"] = data["group"].astype("category")
plt.figure(figsize=(10, 7))
# Scatterplot with color by group
sns.scatterplot(
    data=data,
    x="estimated_MAE",
    y="calculated_age",
    hue="group",
    palette={"young": "green", "middle aged": "yellow", "above middle age": "red"},
    s=120,
    edgecolor="black"
)
# Add D-value labels to each point
for _, row in data.iterrows():
    plt.text(
        row["estimated_MAE"],
        row["calculated_age"],
        f"{row['D_value']:.2f}",
        fontsize=8,
        ha="center",
        va="bottom"
    )
plt.xlabel("Estimated MAE")
plt.ylabel("Predicted Age")
plt.title("Estimated Age vs Estimated MAE From D-value established using age specific data")
plt.legend(title="Age Group")
plt.grid(True, linestyle="--", alpha=0.4)

plt.tight_layout()
plt.show()
# ======================================================================================================================
#
# ======================================================================================================================

# ======================================================================================================================
#
# ======================================================================================================================

# ======================================================================================================================
#
# ======================================================================================================================
