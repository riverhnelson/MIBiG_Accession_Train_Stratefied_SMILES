"""
This script reads the CSV file produced by mibig_npatlas_coconut_smiles.py and outputs a train stratefied SMILES CSV file.
The output file has the MIBiG accession in the first column, followed by one SMILES for each unique compound associated with that accession in columns 2+.
"""

import csv
import tkinter as tk
from tkinter import filedialog


# Priority order for selecting a single SMILES assigned to a single compound. In the Manual SMILES column, TRUE means that no structure was avaialble, so it is ignored.
SMILES_COLUMNS = [
    "NPAtlas SMILES",
    "MIBiG SMILES",
    "COCONUT SMILES",
    "Manual SMILES",
]


def choose_smiles(row):

    # NPAtlas
    smiles = row.get("NPAtlas SMILES", "").strip()

    if smiles:
        return smiles

    # MIBiG
    smiles = row.get("MIBiG SMILES", "").strip()

    if smiles:
        return smiles

    # COCONUT
    smiles = row.get("COCONUT SMILES", "").strip()

    if smiles:
        return smiles

    # Manual
    smiles = row.get("Manual SMILES", "").strip()

    if smiles and smiles.upper() != "TRUE":
        return smiles

    return ""


def main():

    root = tk.Tk()
    root.withdraw()

    # Select input CSV produced by mibig_npatlas_coconut_smiles.py
    input_file = filedialog.askopenfilename(
        title="Select input CSV",
        filetypes=[
            ("CSV files", "*.csv"),
            ("All files", "*.*"),
        ],
    )

    if not input_file:
        return

    # Name output CSV
    output_file = filedialog.asksaveasfilename(
        title="Save output CSV",
        defaultextension=".csv",
        filetypes=[
            ("CSV files", "*.csv"),
            ("All files", "*.*"),
        ],
    )

    if not output_file:
        return

    data = {}

    # Keep the original order of MIBiG accessions.
    accession_order = []

    with open(
        input_file,
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as f:

        reader = csv.DictReader(f)

        for row in reader:

            accession = row.get(
                "MIBiG accession",
                ""
            ).strip()

            compound = row.get(
                "Compound",
                ""
            ).strip()

            if not accession:
                continue

            # Create the accession even if it has no associated compound/SMILES.
            if accession not in data:

                data[accession] = {}

                accession_order.append(accession)

            if not compound:
                continue

            # Select the highest priority SMILES
            smiles = choose_smiles(row)

            if not smiles:
                continue

            # Only one SMILES per compound.
            if compound not in data[accession]:
                data[accession][compound] = smiles

    # Write output CSV with no column headers and the format of: MIBiG accession, SMILES of compound 1, SMILES of compound 2, etc.
    
    with open(
        output_file,
        "w",
        encoding="utf-8",
        newline=""
    ) as f:

        writer = csv.writer(f)

        for accession in accession_order:

            row = [accession]

            # Add one selected SMILES for each unique compound associated with the accession.
            for smiles in data[accession].values():
                row.append(smiles)

            writer.writerow(row)


if __name__ == "__main__":
    main()