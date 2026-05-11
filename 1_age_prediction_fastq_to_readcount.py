
import subprocess
import tkinter as tk
#from tkinter import filedialog as fd, filedialog
import os
import pandas as pd
from pathlib import Path
import sys
BASE = Path(__file__).resolve().parent
ROOT = BASE.parent   # folder where main.py lives
#======================================================================================================================
# Welcome to age_prediction pipeline script.
#======================================================================================================================
#root = tk.Tk() #För manuellt!
#root.withdraw()#För manuellt!
#======================================================================================================================
# Choose directory with samples, CAN ONLY CONTAIN SAMPLE WITH FAST.GZ file format!
#NOTERA!!!!!
#Ingen mapp kan ha mellanslag i namnet, då funkar ej Paths! Alltså ex. "BloodDK" = OK, "Blood DK" = INTE OK!
#
# Notera att scriptet behöver ändras om man ej använder huvud scriptet "Run_all_scripts_Age. får då välja filer manuellt.
#======================================================================================================================
working_folder = BASE
#folder_with_samples = filedialog.askdirectory(title="select folder with the samples")
files = sys.argv[1] # om manuellt stänga av denna.
folder_with_samples = files
print("Files found:", os.listdir(folder_with_samples))
samples = {}
print("Files found:", os.listdir(folder_with_samples))
for filename in os.listdir(folder_with_samples):
    if not filename.endswith(".fastq.gz"):
        continue
    path_to_samples = os.path.join(folder_with_samples, filename)
    if "_R1_" in filename:
        base = filename.replace("_R1_","_X_")
        samples.setdefault(base,{})["R1"] = path_to_samples
    elif "_R2_" in filename:
        base = filename.replace("_R2_", "_X_")
        samples.setdefault(base, {})["R2"] = path_to_samples
#bygger lista av par av samples:
paired_samples = []
sample_number = 1
for base, reads in samples.items():
    if "R1" in reads and "R2" in reads:
        paired_samples.append({
            "sample_number" :sample_number,
            "R1": reads["R1"],
            "R2": reads["R2"]
        })
        sample_number= sample_number + 1
print("Paired Samples table")
print(paired_samples)
table = os.path.join(working_folder,"paired_samples_table.csv")
with open(table, "w", encoding = "utf-8-sig", newline="") as f:
    f.write("sample_number;R1_path;R2_path\n")
    for entry in paired_samples:
        f.write(f"{entry['sample_number']};{entry['R1']};{entry['R2']}\n")
    print(f"Paired sample table saved to {table}")
print("Number of paired samples:", len(paired_samples))
df = pd.read_csv(table, sep = ";", header = 0)
format_table = os.path.join(working_folder, "format_table.csv")
df = df.iloc[:,1:]
df.to_csv(working_folder / "format_table.csv", index=True, header=False)
print (df)
print(df.columns)
Path(table).unlink()
#======================================================================================================================
#
#======================================================================================================================
bed_file = sys.argv[2] # om manuellt stänga av denna.
print(bed_file)
#=============================================================================w=========================================
#
#======================================================================================================================
reference_file = sys.argv[3] # om manuellt stänga av denna.
print(reference_file)
#=============================================================================w=========================================
#======================================================================================================================
#if already made reference genome: delete to remove c2t-created error. Is done only if reference is already created.
#=============================================================================w=========================================
#Ref_bwameth_folder = Path("")
#base_name = ("Age_ref_bloodDK.fasta.bwameth")
#for file in Ref_bwameth_folder.glob(f"{base_name}*"):
#   file.unlink()
#   print("deleted", file)
#======================================================================================================================
#generation of reference genome files t2c types
#=============================================================================w=========================================
print(reference_file)
subprocess.run(["bwameth.py","index",reference_file], check = True)
#======================================================================================================================
#generation of SAM files:
# =============================================================================w=========================================
df = pd.read_csv(BASE/"format_table.csv", header = None, index_col = 0, dtype = {1:str,2:str})
print(df)
sam_folder = Path(BASE/"SAM")
# ======================================================================================================================
# Loop for all samples using Index from formated csv file:
# =============================================================================w=========================================
for i in df.index:
    R1 = df.loc[i,1]
    R2 = df.loc[i,2]
    sample_name = os.path.basename(R1).split("_R1")[0]
    #if "_R1_ in R1: #detta gör så att den hittar tekniska replikat om filen heter något med -r1- och sen _R1_ eller så.
        #sample_name = sample_name + "_r1"
    #else
        #sample_name = sample_name + "_r2"
    cmd = [
        "bwameth.py",
     "-t", "2",
     "--reference",
     reference_file,
     R1,
     R2]
    output_file = sam_folder / f"{sample_name}.sam"
    with open(output_file, "w") as out:
        result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    print("STDERR:", result.stderr)
    print("SAM saved to:",output_file)
    # ======================================================================================================================
    # generation of BAM files:
    # =============================================================================w=========================================
    bam_folder = Path(BASE/"BAM")
    output_file_bam = bam_folder/f"{sample_name}.bam"
    bam_file_sorted = bam_folder/f"{sample_name}.bam"
    bai_file = bam_folder/f"{sample_name}.bam.bai"
    cmd_to_bam = [
        "samtools",
        "view",
        "-b",
        output_file,
        "-o",
        output_file_bam]
    with open(output_file_bam, "wb") as out:
        subprocess.run(cmd_to_bam, stdout=out, check=True)
        print("STDERR:", result.stderr)
        print("Saved to BAM:", output_file_bam)
    # ======================================================================================================================
    # generation of sorted BAM files:
    # =============================================================================w=========================================
    # Bam index
    cmd_sort_bam_make_bai = [
        "samtools",
        "sort",
        "-o",
        bam_file_sorted,
        output_file_bam
    ]
    subprocess.run(cmd_sort_bam_make_bai, check=True)
    # ======================================================================================================================
    # generation of BAI files from sorted BAM files:
    # =============================================================================w=========================================
    cmd_make_bai = [
        "samtools",
        "index",
        bam_file_sorted,
        bai_file
    ]
    subprocess.run(cmd_make_bai)
    print("BAM sorted and indexed successfully.")
    # ======================================================================================================================
    # generation of BAM-readcount
    # =============================================================================w=========================================
    BASE = Path(__file__).resolve().parent
    READCOUNT_DIR = BASE / "READCOUNT"
    READCOUNT_DIR.mkdir(exist_ok=True)
    readcount = READCOUNT_DIR / f"readcount_{sample_name}.txt"

    cmd_bamreadcount = [
        "bam-readcount", # testar ny med base 20, maping q = 10
        "--max-warnings", "0", #0 orginal!
        "--min-base-quality", "30", # 30 orginal!
        "--min-mapping-quality", "30", # 30 orginal!
        "--site-list", bed_file,
        "-f", reference_file,
        bam_file_sorted]
    with open(readcount, "w") as out:
        subprocess.run(cmd_bamreadcount, stdout=out, check=True)
    print(f"Saved readcount_{sample_name}, success")
# ======================================================================================================================

