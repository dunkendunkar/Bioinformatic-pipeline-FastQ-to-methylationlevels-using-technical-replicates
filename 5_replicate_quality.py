
import pandas as pd
#change files to if folder is changed!
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
BASE = Path(__file__).resolve().parent
#files = fd.askopenfilename(title= "select folder which has the fastq-files for your sampels", initialdir= "/home/froste/PycharmProjects/PythonProject/Age/Data")
#df = pd.read_excel(file)
df = pd.read_excel(BASE/"methylation_for_each_sample:gene.xlsx")
# =======================================================================================================================
# Create base
# =======================================================================================================================
df["base"] = (
    df["source"]
    .str.replace(r"readcount_\d+-", "", regex=True)     # tar bort prefixet
    .str.replace(r"_S\d+", "", regex=True)              # tar bort _S4, _S5, _S6 ...
    .str.replace(r"-r[12]", "", regex=True)             # tar bort -r1 eller -r2
    .str.replace(r"_L001\.txt$", "", regex=True)        # tar bort _L001.txt
    .str.replace(r"^\d+-BL-", "", regex=True))
# =======================================================================================================================
# Ensure finding of replicates
# =======================================================================================================================
print(df[["base","source"]])
df["rep"] = df["source"].str.extract(r"(?:-r([12])|rep([12]))").bfill(axis=1).iloc[:,0]
df["rep"] = "r" + df["rep"]
print(df.to_string)
# =======================================================================================================================
# Pass and Fail calculation if pair based on base name differ more the 5% for the methylation
# =======================================================================================================================
results = []
exclude = ["sample", "calculated_age", "base", "rep", "source", "AGE"]
marker_cols = [c for c in df.columns if c not in exclude]
df[marker_cols] = df[marker_cols].apply(pd.to_numeric, errors="coerce")
for base, group in df.groupby("base"):
    r1_rows = group[group["rep"] == "r1"]
    r2_rows = group[group["rep"] == "r2"]
    if r1_rows.empty or r2_rows.empty:
        print(f"Skipping {base}: missing r1 or r2")
        continue
    r1 = r1_rows.iloc[0]
    r2 = r2_rows.iloc[0]
    control = {"base": base, "Status": "Pass"}
    for col in marker_cols:
        diff = abs(r1[col] - r2[col])
        control[col] = diff
        if diff > 7.5:
            control["Status"] = "Fail"
    results.append(control)
control = pd.DataFrame(results)
# ===========================================================================================================================================================
# Control dataframe now contains the correct. Now we remove failed samples to other excel-file. Passed samples are passed on to uncertainty model.
# ============================================================================================================================================================
Control_passed = control[control["Status"] == "Pass"].copy()
Control_failed = control[control["Status"] == "Fail"].copy()
# Showcase reads for failed samples
df_reads = pd.read_excel(BASE/"Save_values_age_anlysis.xlsx")
#extract the source names of all failed samples
Control_failed = Control_failed.merge(df[["base", "source"]],on="base",how="left")
cols = ["gene","C_reads", "T_reads", "total_depth","source"]
Control_failed = Control_failed.merge(df_reads[cols], on="source", how="left")
Control_failed = Control_failed.drop(columns=["source"])
# =================================================================================================================================================================================
# ===========================================================================================================================================================
# New excel file with mean meth values, passed or failed status + calculation of age. Save in Results.
# ============================================================================================================================================================
#calculate mean meth_values in df:
exclude = ["sample", "base", "rep", "source"]
marker_cols = [c for c in df.columns if c not in exclude]
df[marker_cols] = df[marker_cols].apply(pd.to_numeric, errors="coerce")
mean_rows = []
for base, group in df.groupby("base"):
    r1_rows = group[group["rep"] == "r1"]
    r2_rows = group[group["rep"] == "r2"]
    if r1_rows.empty or r2_rows.empty:
        print(f"Skipping {base}: missing r1 or r2")
        continue
    r1 = r1_rows.iloc[0]
    r2 = r2_rows.iloc[0]
    entry = {"base": base}
    for col in marker_cols:
        entry[col] = (r1[col] + r2[col]) / 2
    mean_rows.append(entry)
all = pd.DataFrame(mean_rows)
df2 = []
df2 = all[["base"] + marker_cols]

#
#
#recalculate age:
#
#
markers_coef = {

    #
    #för andra CpG sites i bed-filer, sperma:

    #Modiferad bed-file complex visage blood

    #
    #Blood complex visage
    #
 
    #
    #Next model here
    #
}
AGE_calc = []
for idx, row in df2.iterrows():
    #unmark model you want to use. and # model you do not want to use
    #age = 22.72 #random test för semen visage, changed intercept.
    age = 32.7211426535856 #age semen visage model.
    #age = 28.0 # testar lite random bara!
    #age = 42.2668332790967 #age blood complex visage model
    # extract the row for this sample
    for m in marker_cols:
        methyl = row[m]
        if pd.isna(methyl):
            continue
        if m not in markers_coef:
         continue
        if m == "ELOVL2_435":
                print(m)
                methyl = methyl**2
        coef = markers_coef[m] #sätt gene_pos när indexering är korrekt.
        age = age +  methyl*coef
    AGE_calc.append(age)
df2["calculated_age"] = pd.Series(AGE_calc)
Control_failed = Control_failed.rename(columns={"Status": "Status_failed"})
Control_passed = Control_passed.rename(columns={"Status": "Status_passed"})
df2 = df2.merge(Control_failed[["base", "Status_failed"]], on="base", how="left")
df2 = df2.merge(Control_passed[["base", "Status_passed"]], on="base", how="left")
df2["Status_replicates"] = df2["Status_passed"].fillna(df2["Status_failed"])
df2 = df2.drop(columns=["Status_passed", "Status_failed"])
df2 = df2.drop_duplicates(subset="base")
df2.to_excel(BASE/"meth_calcAge_combinedReplicates.xlsx",index = True) # fixing of coverage values
print(control.to_string())
df_plot = control
print("COLUMNS IN DF_PLOT:", df_plot.columns.tolist())
marker_cols = [c for c in df_plot.columns if c not in ["base", "Status"]]
df_plot = df_plot.melt(
    id_vars=["base", "Status"],
    value_vars=marker_cols,
    var_name="Marker",
    value_name="%_meth"
)
df_plot["Color"] = df_plot["Status"].map({
    "Pass": "green",
    "Fail": "red"
})
plt.figure(figsize=(14, 6))
sns.boxplot(
    data=df_plot,
    x="Marker",
    y="%_meth",
    color="lightgrey",
    fliersize=3,
    linewidth=1,width= 0.6

)
sns.stripplot(
    data=df_plot,
    x="Marker",
    y="%_meth",
    size=4,
    alpha=1,
    color="grey",  # <— THIS MAKES DOTS WHITE
    edgecolor="black"  # optional: outline so they don’t disappear
)

plt.xticks(rotation=90)
plt.title("Methylation deviation per Marker (Replicates)", fontsize=14)
plt.ylabel("Methylation (%)", fontsize=12)
plt.xlabel("Marker", fontsize=12)
plt.tight_layout()
plt.legend(title="")
ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.savefig(BASE/"RESULTS/REPLICATED_DEVATION_PER_MARKER")
#
