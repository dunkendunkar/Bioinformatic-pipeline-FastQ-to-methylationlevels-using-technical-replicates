import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#IFRÅN GAMLA
from pathlib import Path
from pathlib import Path
import shutil
import os
import shutil
import os
from pathlib import Path

# Base folder where Pipeline_start.py lives
BASE = Path(__file__).resolve().parent
CODE_DIR = BASE.parent
# Ask user for project folder name
folder_name = input("Enter folder name: ")
full_path = CODE_DIR / folder_name
full_path.mkdir(exist_ok=True)
# Subfolders to create inside the new project folder
subfolders = ["SAM", "BAM", "READCOUNT", "RESULTS"]
for sub in subfolders:
    (full_path / sub).mkdir(exist_ok=True)
print(f"Created project folder '{full_path}' with subfolders:")
for sub in subfolders:
    print(f" - {sub}")
# =====================================================================
# Function to copy files from source → target and delete originals
# =====================================================================
def move_files(src_folder, dst_folder, pattern=None):
    src_folder = Path(src_folder)
    dst_folder = Path(dst_folder)
    dst_folder.mkdir(exist_ok=True)
    for filename in os.listdir(src_folder):
        if pattern and not filename.startswith(pattern):
            continue
        src = src_folder / filename
        dst = dst_folder / filename
        if src.is_file():
            shutil.copy2(src, dst)
            src.unlink()
            print(f"Moved: {src} → {dst}")

# =====================================================================
# Copy + clean SAM, BAM, READCOUNT, RESULTS
# =====================================================================
move_files(BASE / "SAM",        full_path / "SAM")
move_files(BASE / "BAM",        full_path / "BAM")
move_files(BASE / "READCOUNT",  full_path / "READCOUNT", pattern="readcount")
move_files(BASE / "RESULTS",    full_path / "RESULTS")

# Copy MultiQC directory if it exists
multiqc_src = BASE / "RESULTS" / "MultiQC"
multiqc_dst = full_path / "RESULTS" / "MultiQC"
if multiqc_src.exists() and multiqc_src.is_dir():
    shutil.copytree(multiqc_src, multiqc_dst, dirs_exist_ok=True)
    shutil.rmtree(multiqc_src)
    print("Moved MultiQC directory.")
else:
    print("No MultiQC directory found.")

# =====================================================================
# Remove original folders completely (optional but recommended)
# =====================================================================
folders_to_remove = ["SAM", "BAM", "READCOUNT", "RESULTS", "MULTIQC"]
for folder in folders_to_remove:
    folder_path = BASE / folder
    if folder_path.exists():
        shutil.rmtree(folder_path)
        print(f"Removed folder: {folder_path}")

figure_files = [f for f in os.listdir(full_path / "RESULTS")
    if f.lower().endswith((".png", ".jpg", ".jpeg"))]

excel_html =[full_path /"Save_values_age_anlysis.xlsx",
full_path /"methylation_for_each_sample:gene.xlsx",
full_path /"meth_calcAge_combinedReplicates.xlsx"]

html_path= full_path / "RESULTS" / "age_report.html"
# --- 3. Build HTML ---
html = """
<html>
<head>
    <title>Analysis Report</title>
    <style>
        body { font-family: Arial; margin: 20px; }
        img { max-width: 900px; margin: 20px 0; border: 1px solid #ccc; display: block; }
        table { border-collapse: collapse; margin-top: 40px; width: 95%; }
        th, td { border: 1px solid #999; padding: 6px; font-size: 12px; }
        th { background: #eee; }
        h2 { margin-top: 50px; }
    </style>
</head>
<body>
<h1>Analysis Report</h1>
<p>This report shows all generated figures and Excel tables for the samples.</p>

<h2>Figures</h2>
"""
# --- 4. Insert images using RELATIVE paths ---
figure_files = sorted(figure_files)
for name in figure_files:
    rel_path = name  # HTML file is in the same folder as the images
    html += f'<img src="{rel_path}" alt="{name}">\n'
# --- 5. Insert Excel tables ---
html += "<h2>Excel Tables</h2>\n"
for path in excel_html:
    try:
        df_tmp = pd.read_excel(path)
        html += f"<h3>{path.name}</h3>\n"
        html += df_tmp.to_html(index=False)
    except Exception as e:
        html += f"<p>Error loading {path.name}: {e}</p>"
# --- 6. Close HTML ---
html += """
</body>
</html>
"""
# --- 7. Write HTML file ---
with open(html_path, "w", encoding="utf-8") as f:
    f.write(html)
print("HTML report created at:", html_path)


excel_files = [
    BASE / "Save_values_age_anlysis.xlsx",
    BASE / "methylation_for_each_sample:gene.xlsx",
    BASE / "meth_calcAge_combinedReplicates.xlsx"
]
for src in excel_files:
    dst = full_path /"RESULTS"/ src.name
    if src.exists():
        shutil.copy2(src, dst)
        print(f"Copied Excel file: {src.name}")
    else:
        print(f"WARNING: Excel file not found: {src}")


