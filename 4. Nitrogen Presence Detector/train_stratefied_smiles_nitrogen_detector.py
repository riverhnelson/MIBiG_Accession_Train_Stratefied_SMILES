import csv
import tkinter as tk
from tkinter import filedialog
from rdkit import Chem


def main():

    # Create file picker window
    root = tk.Tk()
    root.withdraw()

    # Select input CSV
    input_file = filedialog.askopenfilename(
        title="Select input CSV",
        filetypes=[
            ("CSV files", "*.csv"),
            ("All files", "*.*")
        ]
    )

    if not input_file:
        print("No input file selected.")
        return

    # Select output CSV location
    output_file = filedialog.asksaveasfilename(
        title="Save output CSV",
        defaultextension=".csv",
        filetypes=[
            ("CSV files", "*.csv"),
            ("All files", "*.*")
        ]
    )

    if not output_file:
        print("No output file selected.")
        return

    results = []

    # Read input CSV
    with open(
        input_file,
        "r",
        newline="",
        encoding="utf-8"
    ) as infile:

        reader = csv.reader(infile)

        for row in reader:

            # Skip empty rows
            if not row:
                continue

            mibig_accession = row[0].strip()

            # Track whether this accession has at least one SMILES
            has_smiles = False

            # Assume no nitrogen until one is found
            contains_nitrogen = 0

            # Check every SMILES associated with the accession
            for smiles in row[1:]:

                smiles = smiles.strip()

                # Skip empty SMILES strings
                if not smiles:
                    continue

                # At least one SMILES exists
                has_smiles = True

                # Convert SMILES to RDKit molecule
                molecule = Chem.MolFromSmiles(smiles)

                # Warn if RDKit cannot parse the SMILES
                if molecule is None:
                    print(
                        f"WARNING: could not parse SMILES "
                        f"for {mibig_accession}:"
                    )
                    print(smiles)
                    continue

                # Check whether this compound contains nitrogen
                has_nitrogen = any(
                    atom.GetSymbol() == "N"
                    for atom in molecule.GetAtoms()
                )

                # If at least one compound contains nitrogen,
                # the accession is classified as 1.
                if has_nitrogen:
                    contains_nitrogen = 1
                    break

            # If there were no SMILES at all, skip this accession
            if not has_smiles:
                continue

            results.append([
                mibig_accession,
                contains_nitrogen
            ])

    # Write output CSV
    with open(
        output_file,
        "w",
        newline="",
        encoding="utf-8"
    ) as outfile:

        writer = csv.writer(outfile)

        writer.writerow([
            "MIBiG accession",
            "Chemical products contain nitrogen (1 = TRUE and 0 = FALSE)"
        ])

        writer.writerows(results)

    print()
    print("Nitrogen presence analysis complete.")
    print("Output written to:")
    print(output_file)


if __name__ == "__main__":
    main()