
import subprocess
import tkinter as tk
from tkinter import filedialog as fd, filedialog
import os
import pandas as pd
from pathlib import Path
print("QUALITY TESTING WITH MULTIQC!")
import sys
root = tk.Tk()
root.withdraw()
BASE = Path(__file__).resolve().parent
# ======================================================================================================================
#folder creation. Note selction of folder can contain other files then bam files.
# =============================================================================w=========================================
folder = BASE/"RESULTS"
#files = filedialog.askdirectory(title= "select folder which has the BAM-files for your sampels")
files = BASE/"BAM"
fastqc_out = BASE/"RESULTS/MultiQC/FastQC"
output = BASE/"RESULTS/MultiQC"
# Create all directories if they don't exist
os.makedirs(folder, exist_ok=True)
os.makedirs(files, exist_ok=True)
os.makedirs(fastqc_out, exist_ok=True)
os.makedirs(output, exist_ok=True)
# ======================================================================================================================
# Loop for all files in directory with .bam in name.
# =============================================================================w=========================================
for file in os.listdir(files):
    if file.endswith(".bam"):
        path_file = os.path.join(files, file)
        print(f"running FastQC on {path_file}")
        subprocess.run(["fastqc","-o" ,fastqc_out,path_file])
# ======================================================================================================================
#summary for all fastqc html report with a multiqc report
# =============================================================================w=========================================
subprocess.run(["multiqc",fastqc_out,"-o",output ])
print(f"generated MultiQC report in {output}")
# ======================================================================================================================
#remove all uncesseracry fastq files.
# =============================================================================w=========================================
for filename in fastqc_out:
    file_path = os.path.join(fastqc_out, filename)
    if os.path.isfile(file_path):
        os.remove(file_path)
# ======================================================================================================================
#Complete MULTIQC FOR BAM FILES
# =============================================================================w=========================================
print("Completed MultiQC for Samples with BAM files!")