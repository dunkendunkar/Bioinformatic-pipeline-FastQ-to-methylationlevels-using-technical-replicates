import pandas as pd
from patsy.contrasts import Sum
from pathlib import Path
BASE = Path(__file__).resolve().parent
df= pd.read_excel(BASE/"meth_calcAge_combinedReplicates.xlsx") # fixing of coverage values
df_all = pd.read_excel(BASE/"Save_values_age_anlysis.xlsx") # for reads
df_all["total_depth"] = pd.to_numeric(df_all["total_depth"], errors="coerce")
print(df.to_string())
# =======================================================================================================================
# =======================================================================================================================
# Control the reads for each marker: 1000 reads/marker per sample -> Fail.
# =======================================================================================================================
# =======================================================================================================================
def clean_basename(s):
    return (
        s.str.replace(r"^readcount_\d+-", "", regex=True)
         .str.replace(r"^\d+-BL-", "", regex=True)
         .str.replace(r"_S\d+_L001\.txt$", "", regex=True)
         .str.replace(r"-r[12]$", "", regex=True)
    )
df_all["basename"] = clean_basename(df_all["source"])
df["basename"] = clean_basename(df["base"])
print(df["basename"])
print(df_all["basename"])
qc1 = (
    df_all.groupby(["basename", "gene"])["total_depth"]
    .min()
    .reset_index()
)
qc1["status_reads_marker"] = qc1["total_depth"].apply(
    lambda x: "pass" if x >= 1000 else "fail"
)

# collapse to one row per sample
qc1 = qc1.groupby("basename")["status_reads_marker"].apply(
    lambda x: "fail" if "fail" in x.values else "pass"
).reset_index()
# =======================================================================================================================
# =======================================================================================================================
# =======================================================================================================================
# =======================================================================================================================
# Control the total reads for each sample: less then total reads of 2000xmarkers in model -> fail
# =======================================================================================================================
# =======================================================================================================================
#
qc2 = (
    df_all.groupby("basename")["total_depth"]
    .sum()
    .reset_index()
)
qc2["status_reads_tot"] = qc2["total_depth"].apply(
    lambda x: "pass" if x > 25000 else "fail"
)
print(qc2)
# =======================================================================================================================
# =======================================================================================================================
# =======================================================================================================================
# =======================================================================================================================
# Remove old QC columns if they exist
df = df.drop(columns=[c for c in df.columns if c.startswith("status_reads")], errors="ignore")
df = df.merge(qc1, on="basename", how="left")
df = df.merge(qc2, on="basename", how="left")
df.drop_duplicates(subset="basename",inplace=True)
print(df.to_string())
BASE = Path(__file__).resolve().parent
df.to_excel(BASE/"meth_calcAge_combinedReplicates.xlsx",index = True) # fixing of coverage values
# =======================================================================================================================
# =======================================================================================================================
# PLOT SOME READS for reads/marker and total reads!
# =======================================================================================================================
# =======================================================================================================================
