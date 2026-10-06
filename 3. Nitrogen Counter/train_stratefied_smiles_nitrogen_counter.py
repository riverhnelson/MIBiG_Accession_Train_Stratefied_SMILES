import csv
from rdkit import Chem

# Input CSV file containing MIBiG accessions and SMILES strings and output CSV containing min/max nitrogen counts for each BGC
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()

input_file = filedialog.askopenfilename(
    title="Select input CSV",
    filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
)

if not input_file:
    exit()

output_file = filedialog.asksaveasfilename(
    title="Save output CSV",
    defaultextension=".csv",
    filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
)

if not output_file:
    exit()

with open(input_file, "r", newline="", encoding="utf-8") as infile:
    reader = csv.reader(infile)

    results = []

    for row in reader:

        # Skip empty rows
        if not row:
            continue

        mibig_accession = row[0]

        nitrogen_counts = []

        for smiles in row[1:]:

            # Skip empty SMILES strings
            if not smiles:
                continue

            # Convert SMILES to RDKit molecule object
            molecule = Chem.MolFromSmiles(smiles)

            # Print warning for SMILES that cannot be parsed
            if molecule is None:
                print(f"WARNING: could not parse SMILES for {mibig_accession}:")
                print(smiles)
                continue

            # Count nitrogen atoms in each smiles string
            nitrogen_count = sum(
                1 for atom in molecule.GetAtoms()
                if atom.GetSymbol() == "N"
            )

            nitrogen_counts.append(nitrogen_count)

        # If there are valid nitrogen counts, store min and max for each MIBiG accession
        if nitrogen_counts:
            min_nitrogen = min(nitrogen_counts)
            max_nitrogen = max(nitrogen_counts)

            results.append([
                mibig_accession,
                min_nitrogen,
                max_nitrogen
            ])

# Write results to output CSV file with column headers
with open(output_file, "w", newline="", encoding="utf-8") as outfile:
    writer = csv.writer(outfile)

    writer.writerow([
        "MIBiG accession",
        "Min chemical product nitrogen count",
        "Max chemical product nitrogen count"
    ])

    writer.writerows(results)


print(f"Nitrogen count output written to: {output_file}")