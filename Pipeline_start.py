# =======================================================================================================================
# Running all the different script from one script! START HERE:
# =======================================================================================================================
# =======================================================================================================================
# Packages
# =======================================================================================================================
# =======================================================================================================================
import tkinter as tk
import pandas as pd
from tkinter import filedialog as fd, filedialog
import os
import pandas as pd
from pathlib import Path
root = tk.Tk()
from tkinter import messagebox
root.withdraw()
import subprocess
# =======================================================================================================================
#Controll Necessary Directories exists: SAM, BAM, READCOUNT, RESULT: ADD NAMES HERE IF ANYTHING MORE IS NEEDED
# =======================================================================================================================
#necessary = ["Code/BAM","Code/SAM","Code/READCOUNT", "Code/RESULTS","Code/MultiQC_bam"]   # example folder name
# Absolute path to the directory where the script is located
#script_dir = os.path.dirname(os.path.abspath("Code"))
# Create each folder safely
#for name in necessary:
#    path = os.path.join(script_dir, name)
#    os.makedirs(path, exist_ok=True)
#excel_files = [ "meth_calcAge_combinedReplicates.xlsx",
#                 "methylation_for_each_sample:gene.xlsx",
#                 "Save_values_age_anlysis.xlsx",
#                 "RESULTS/Bisulfate_controll_readcounts.xlsx"]
#BASE = os.path.dirname(os.path.abspath(__file__))
#CODE_DIR = os.path.join(BASE, "Code")
#for filename in excel_files:
#    path = os.path.join(CODE_DIR, filename)
#    if not os.path.exists(path):
#        os.makedirs(os.path.dirname(path), exist_ok=True)
#        pd.DataFrame().to_excel(path, index=False)
#csv_files = ["paired_samples_table.csv", "format_table.csv"]
#for filename in csv_files:
#    path = os.path.join(CODE_DIR, filename)
#    os.makedirs(os.path.dirname(path), exist_ok=True)
#    if not os.path.exists(path):
#        pd.DataFrame().to_csv(path, index=False)
#
# =======================================================================================================================
necessary = ["BAM", "SAM", "READCOUNT", "RESULTS", "MultiQC_bam"]
# Absolute path to the directory where the script is located
BASE = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.join(BASE, "Code")
# Create each folder safely
for name in necessary:
    os.makedirs(os.path.join(CODE_DIR, name), exist_ok=True)
excel_files = [ "meth_calcAge_combinedReplicates.xlsx",
                 "methylation_for_each_sample:gene.xlsx",
                 "Save_values_age_anlysis.xlsx",
                 "RESULTS/Bisulfate_controll_readcounts.xlsx"]

for filename in excel_files:
    path = os.path.join(CODE_DIR, filename)
    folder = os.path.dirname(path)
    os.makedirs(folder, exist_ok=True)
    if not os.path.exists(path):
        pd.DataFrame().to_excel(path, index=False)
csv_files = ["paired_samples_table.csv", "format_table.csv"]
for filename in csv_files:
    path = os.path.join(CODE_DIR, filename)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if not os.path.exists(path):
        pd.DataFrame().to_csv(path, index=False)
#
# =======================================================================================================================
# =======================================================================================================================
# Selection of part of script
# =======================================================================================================================
# =======================================================================================================================
#
# select if you want to run multiQC for BAM and FastQ files!
answer = messagebox.askyesno("MultiQC", "analysis for FastQ and the generated BAM files? Select Yes or No")
# select if you want to run FastQ->READCOUNT
answer_2 = messagebox.askyesno("FastQ->Readcount", "FastQ->readcount or start in READCOUNT folder? Select Yes for starting with FastQ files or No for starting in READCOUNT FOLDER")
#
# =======================================================================================================================
# =======================================================================================================================
# Select FastQ, Fasta, Bed reference files
# =======================================================================================================================
# =======================================================================================================================
#Select FastQ-files
files = filedialog.askdirectory(title= "select folder which has the fastq-files for your sampels", initialdir= "/home/froste/PycharmProjects/PythonProject/Age/Data")
#select FASTA
ref_folder = os.path.expanduser("/home/froste/PycharmProjects/Age")
reference_file = fd.askopenfilename(initialdir= "/home/froste/PycharmProjects/PythonProject/Age/Data", title ="Chose right reference genome for model, MUST BE FASTA file type", )
#select BED file
model_folder = os.path.expanduser("/home/froste/PycharmProjects/")
bed_file = fd.askopenfilename(initialdir= "/home/froste/PycharmProjects/PythonProject/Age/Data", title="Chose right BED_file for model.")
# Folder where main.py lives
# Path to the Code/ folder
base = CODE_DIR
# =======================================================================================================================
# =======================================================================================================================
# Run fastQ -> Readcount, Selected Yes or skipping to start in Readcount script
# =======================================================================================================================
# =======================================================================================================================
#
script_1 = os.path.join(base,"1_age_prediction_fastq_to_readcount.py")
if answer_2 == True:
    subprocess.run(["python",script_1 , files, bed_file, reference_file])
else:
    print("Skipping to READCOUNT folder results")
# =======================================================================================================================
# =======================================================================================================================
# Readcount -> Excel -> Save_age_values -> Methylation_levels_excel with estimated age
# =======================================================================================================================
# =======================================================================================================================
#
script_2 = os.path.join(base,"2_age_analysis_using_readcount.py")
subprocess.run(["python",script_2 ])
# =======================================================================================================================
# =======================================================================================================================
# MultiQC -> Selected Yes for analysis, No for skipping of it.
# =======================================================================================================================
# =======================================================================================================================
#
script_3 = os.path.join(base,"3_MultiQC_for_bam_files.py")
if answer:
    subprocess.run(["python",script_3]) #MultiQ_bam
else:
    print("Skipping MultiQC analysis.")
#
# =======================================================================================================================
# =======================================================================================================================
# Bisulfate controll
# =======================================================================================================================
# =======================================================================================================================
#
script_4 = os.path.join(base,"4_Bisulfat_Control.py")
subprocess.run(["python",script_4 , reference_file])
#
# =======================================================================================================================
# =======================================================================================================================
# Replication control for samples: Returns passed samples in methylation_levels, failed samples in Failed_samples_table.
# =======================================================================================================================
# =======================================================================================================================
#
script_5 = os.path.join(base,"5_replicate_quality.py")
subprocess.run(["python",script_5])
# =======================================================================================================================
# =======================================================================================================================
# Reads/Marker and total/reads_sample for alla markers in model: return Pass or fail. removes failed samples to Failed_samples Table
# =======================================================================================================================
# =======================================================================================================================
#
script_6 = os.path.join(base,"6_Control_reads.py")
subprocess.run(["python",script_6])
# =======================================================================================================================
# =======================================================================================================================
# =======================================================================================================================
#Predicted MAE for each sample
# =======================================================================================================================
# =======================================================================================================================
#
script_7 = os.path.join(base,"7_Uncertantity_by_methylation_model.py")
subprocess.run(["python",script_7 ])
# =======================================================================================================================
# =======================================================================================================================
#HTML-Report: Excel with predicted age + figures +
# =======================================================================================================================
# =======================================================================================================================
#
script_8 = os.path.join(base, "8_Folder_HTML_report.py")
subprocess.run(["python",script_8])
# =======================================================================================================================
# =======================================================================================================================
# Generated Result folder: Saves readcounts.txt, BAM + bai, Reads/markers,, predicted MAE in CSV + graph
# =======================================================================================================================
# =======================================================================================================================
print("Success! Analysis Complete")
# =======================================================================================================================
# =======================================================================================================================

