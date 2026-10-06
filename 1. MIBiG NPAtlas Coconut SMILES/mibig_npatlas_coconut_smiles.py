"""
Add COCONUT database CSV and MIBiG 4.0 database JSON SMILES to train_stratefied_smiles.csv.
Compounds that do not have SMILES will be marked as "Manual SMILES" = TRUE
DOI and PMID will be extracted from MIBiG JSON and added to the output. Links to DOI and PMID are also added to the output.
Compound names will be filled in from NPAtlas_as.tsv if the NPAtlas SMILES matches exactly.

INPUTS:
1. train_stratefied_smiles.csv: column 1 = MIBiG accession; columns 2+ = NPAtlas SMILES
2. Folder containing MIBiG 4.0 JSON files (mibig_json_4.0): used to add MIBiG SMILES
3. NPAtlas_as.tsv: used to fill in missing compound names by text matching to NPAtlas SMILES
4. COCONUT database CSV (coconut_csv-10_2026.csv): used to add COCONUT smiles

OUTPUT CSV COLUMNS:
Compound
NPAtlas SMILES
MIBiG SMILES
COCONUT SMILES
Manual SMILES
DOI
PMID
Reference links
"""

import csv
import json
import re
import sys
from pathlib import Path
from rdkit import Chem
import tkinter as tk
from tkinter import filedialog, messagebox


# Required due to the large CSV fields of the COCONUT database.
csv.field_size_limit(sys.maxsize)


OUTPUT_COLUMNS = [
    "MIBiG accession",
    "Compound",
    "NPAtlas SMILES",
    "MIBiG SMILES",
    "COCONUT SMILES",
    "Manual SMILES",
    "DOI",
    "PMID",
    "Reference links",
]


DOI_RE = re.compile(
    r"10\.\d{4,9}/[-._;()/:\w]+",
    re.I,
)


PMID_RE = re.compile(
    r"(?:pubmed:|pmid:|pubmed\.ncbi\.nlm\.nih\.gov/)"
    r"(\d{6,9})",
    re.I,
)


def clean(value):
    return "" if value is None else str(value).strip()

# Canonicalize SMILES using RDKit. Returns None if the SMILES is invalid.
def canonical_smiles(smiles):
   

    smiles = clean(smiles)

    if not smiles:
        return None

    try:
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            return None

        return Chem.MolToSmiles(
            mol,
            isomericSmiles=True,
            canonical=True,
        )

    except Exception:
        return None


# Read the input train_stratefied_smiles.csv and create a mapping of MIBiG accession to SMILES
# train_stratefied_smiles.csv format:
# column 0 = MIBiG accession
# columns 1+ = NPAtlas SMILES
def load_input_csv(path):
    result = {}

    with path.open(
        "r",
        newline="",
        encoding="utf-8-sig",
    ) as f:

        reader = csv.reader(f)

        for line_no, row in enumerate(reader, 1):

            if not row:
                continue

            accession = clean(row[0])

            if not accession:
                continue

            if line_no == 1 and accession.lower() in {
                "mibig accession",
                "mibig_accession",
                "accession",
            }:
                continue

            result.setdefault(
                accession,
                [],
            )

            for cell in row[1:]:

                smiles = clean(cell)

                if (
                    smiles
                    and smiles not in result[accession]
                ):
                    result[accession].append(
                        smiles
                    )

    return result


# Read NPAtlas_as.tsv and create a mapping of exact NPAtlas SMILES to compound name. Assigns compound names to NPAtlas SMILES.
# NPAtlas_as.tsv format:
# column 0  = NPAtlas accession
# column 2  = compound name
# column 10 = SMILES
def load_npatlas_tsv(path):
    smiles_to_name = {}

    with path.open(
        "r",
        newline="",
        encoding="utf-8-sig",
    ) as f:

        reader = csv.reader(
            f,
            delimiter="\t",
        )

        for row in reader:

            if not row or len(row) <= 10:
                continue

            name = clean(row[2])
            smiles = clean(row[10])

            if (
                smiles
                and name
                and smiles not in smiles_to_name
            ):
                smiles_to_name[smiles] = name

    return smiles_to_name


# Read COCONUT database coconut_csv-10_2026.csv and map compound names to COCONUT SMILES.
# coconut_csv-10_2026.csv format:
# column 5 = compound name
# column 2 = COCONUT SMILES
def load_coconut_csv(path):
    
    name_to_smiles = {}

    with path.open(
        "r",
        newline="",
        encoding="utf-8-sig",
    ) as f:

        reader = csv.DictReader(f)

        if reader.fieldnames is None:
            raise ValueError(
                "COCONUT CSV does not contain a header row."
            )

        required_columns = {
            "name",
            "canonical_smiles",
        }

        missing = (
            required_columns
            - set(reader.fieldnames)
        )

        if missing:
            raise ValueError(
                "COCONUT CSV is missing required "
                "column(s): "
                + ", ".join(
                    sorted(missing)
                )
            )

        for row in reader:

            name = clean(
                row.get("name")
            )

            smiles = clean(
                row.get("canonical_smiles")
            )

            if not name or not smiles:
                continue

            # CASE-INSENSITIVE MATCHING
            key = name.casefold()

            if key not in name_to_smiles:
                name_to_smiles[key] = smiles

    return name_to_smiles


def find_json(folder, accession):
    exact = folder / f"{accession}.json"

    if exact.is_file():
        return exact

    accession_lower = accession.lower()

    for p in folder.glob("*.json"):

        if p.stem.lower() == accession_lower:
            return p

    return None


# Get DOI and PubMed identifiers from MIBiG database JSON.
def extract_references(obj):
    dois = set()
    pmids = set()

    def walk(x):

        if isinstance(x, dict):

            for value in x.values():
                walk(value)

        elif isinstance(x, list):

            for value in x:
                walk(value)

        elif isinstance(x, str):

            for match in DOI_RE.findall(x):

                dois.add(
                    match.rstrip(
                        ".,;)]}"
                    )
                )

            for match in PMID_RE.findall(x):

                pmids.add(match)

    walk(obj)

    return (
        sorted(dois),
        sorted(pmids),
    )

# Add DOI and PubMed URLs to the output.
def make_reference_links(
    dois,
    pmids,
):
   
    links = []

    for doi in dois:

        links.append(
            f"https://doi.org/{doi}"
        )

    for pmid in pmids:

        links.append(
            f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        )

    return links


def read_mibig_json(path):

    with path.open(
        "r",
        encoding="utf-8",
    ) as f:

        data = json.load(f)

    accession = clean(
        data.get(
            "accession",
            path.stem,
        )
    )

    compounds = []

    for compound in data.get(
        "compounds",
        [],
    ):

        if not isinstance(
            compound,
            dict,
        ):
            continue

        name = clean(
            compound.get("name")
        )

        smiles = clean(
            compound.get("structure")
            or compound.get("smiles")
        )

        compounds.append(
            {
                "name": name,
                "smiles": smiles,
                "canonical": canonical_smiles(
                    smiles
                ),
            }
        )

    dois, pmids = extract_references(
        data
    )

    return (
        accession,
        compounds,
        dois,
        pmids,
    )

# RDKit matching to match identical NPAtlas and MIBiG SMILES.
def compare_accession(
    accession,
    np_smiles,
    mibig_compounds,
    dois,
    pmids,
):

    np_records = []

    for i, smiles in enumerate(
        np_smiles
    ):

        np_records.append(
            {
                "index": i,
                "smiles": smiles,
                "canonical": canonical_smiles(
                    smiles
                ),
                "matched": False,
            }
        )

    matched_np_for_mibig = {}

    for mi, compound in enumerate(
        mibig_compounds
    ):

        mc = compound["canonical"]

        if not mc:
            continue

        for np in np_records:

            if np["matched"]:
                continue

            if np["canonical"] == mc:

                matched_np_for_mibig[mi] = np

                np["matched"] = True

                break

    rows = []

    for mi, compound in enumerate(
        mibig_compounds
    ):

        mibig_smiles = compound[
            "smiles"
        ]

        matched_np = (
            matched_np_for_mibig.get(mi)
        )

        np_smiles_value = (
            matched_np["smiles"]
            if matched_np
            else ""
        )

        manual = (
            not np_smiles_value
            and not mibig_smiles
        )

        rows.append(
            {
                "MIBiG accession": accession,
                "Compound": compound[
                    "name"
                ],
                "NPAtlas SMILES": np_smiles_value,
                "MIBiG SMILES": mibig_smiles,
                "COCONUT SMILES": "",
                "Manual SMILES": (
                    "TRUE"
                    if manual
                    else ""
                ),
                "DOI": "; ".join(dois),
                "PMID": "; ".join(pmids),
                "Reference links": "; ".join(
                    make_reference_links(
                        dois,
                        pmids,
                    )
                ),
            }
        )

    for np in np_records:

        if np["matched"]:
            continue

        rows.append(
            {
                "MIBiG accession": accession,
                "Compound": "",
                "NPAtlas SMILES": np["smiles"],
                "MIBiG SMILES": "",
                "COCONUT SMILES": "",
                "Manual SMILES": "",
                "DOI": "; ".join(dois),
                "PMID": "; ".join(pmids),
                "Reference links": "; ".join(
                    make_reference_links(
                        dois,
                        pmids,
                    )
                ),
            }
        )

    if (
        not mibig_compounds
        and np_smiles
    ):

        for smiles in np_smiles:

            rows.append(
                {
                    "MIBiG accession": accession,
                    "Compound": "",
                    "NPAtlas SMILES": smiles,
                    "MIBiG SMILES": "",
                    "COCONUT SMILES": "",
                    "Manual SMILES": "",
                    "DOI": "; ".join(dois),
                    "PMID": "; ".join(pmids),
                    "Reference links": "; ".join(
                        make_reference_links(
                            dois,
                            pmids,
                        )
                    ),
                }
            )

    return rows


# Fill blank compound names using NPAtlas SMILES matching to compound names with NPAtlas_as.tsv.
def fill_missing_compound_names(
    rows,
    npatlas_smiles_to_name,
):
    
    filled = 0

    for row in rows:

        if clean(
            row.get("Compound")
        ):
            continue

        smiles = clean(
            row.get("NPAtlas SMILES")
        )

        if not smiles:
            continue

        if (
            smiles
            in npatlas_smiles_to_name
        ):

            row["Compound"] = (
                npatlas_smiles_to_name[
                    smiles
                ]
            )

            filled += 1

    return filled


# Merge semicolon-separated values while preserving order and removing duplicates.
def _merge_semicolon_values(
    existing,
    new,
):
    
    values = []

    for source in (
        existing,
        new,
    ):

        for value in source.split(";"):

            value = value.strip()

            if (
                value
                and value not in values
            ):
                values.append(value)

    return "; ".join(values)


# Group rows by MIBiG accession and compound name (case-insensitive).
def group_rows_by_compound_name(rows):

    grouped = {}
    unnamed_rows = []

    for row in rows:

        accession = clean(
            row.get(
                "MIBiG accession"
            )
        )

        name = clean(
            row.get("Compound")
        )

        if not name:

            unnamed_rows.append(row)

            continue

        # CASE-INSENSITIVE GROUPING
        key = (
            accession.casefold(),
            name.casefold(),
        )

        if key not in grouped:

            grouped[key] = row.copy()

            continue

        existing = grouped[key]

        if not existing["Compound"]:

            existing["Compound"] = name

        if (
            row["NPAtlas SMILES"]
            and row["NPAtlas SMILES"]
            != existing["NPAtlas SMILES"]
        ):

            if not existing[
                "NPAtlas SMILES"
            ]:

                existing[
                    "NPAtlas SMILES"
                ] = row[
                    "NPAtlas SMILES"
                ]

        if (
            row["MIBiG SMILES"]
            and row["MIBiG SMILES"]
            != existing["MIBiG SMILES"]
        ):

            if not existing[
                "MIBiG SMILES"
            ]:

                existing[
                    "MIBiG SMILES"
                ] = row[
                    "MIBiG SMILES"
                ]

            elif (
                row["MIBiG SMILES"]
                not in existing[
                    "MIBiG SMILES"
                ].split(" | ")
            ):

                existing[
                    "MIBiG SMILES"
                ] += (
                    " | "
                    + row[
                        "MIBiG SMILES"
                    ]
                )

        existing[
            "Manual SMILES"
        ] = (
            "TRUE"
            if not clean(
                existing[
                    "NPAtlas SMILES"
                ]
            )
            and not clean(
                existing[
                    "MIBiG SMILES"
                ]
            )
            else ""
        )

        for column in [
            "DOI",
            "PMID",
            "Reference links",
        ]:

            existing[column] = (
                _merge_semicolon_values(
                    existing[column],
                    row[column],
                )
            )

    return (
        list(grouped.values())
        + unnamed_rows
    )


# Match (case-insensitive) compound names to COCONUT database compound names and add the corresponding COCONUT canonical SMILES to the output.
def add_coconut_smiles(
    rows,
    coconut_name_to_smiles,
):
    
    matched = 0

    for row in rows:

        compound_name = clean(
            row.get("Compound")
        )

        if not compound_name:
            continue

        # CASE-INSENSITIVE MATCHING
        key = compound_name.casefold()

        coconut_smiles = (
            coconut_name_to_smiles.get(
                key,
                "",
            )
        )

        if coconut_smiles:

            row["COCONUT SMILES"] = (
                coconut_smiles
            )

            # COCONUT now provides a structure,
            # so manual structure determination is
            # no longer required.
            row[
                "Manual SMILES"
            ] = ""

            matched += 1

    return matched


def write_output(
    rows,
    path,
):

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=OUTPUT_COLUMNS,
        )

        writer.writeheader()

        writer.writerows(rows)


def choose_files():

    root = tk.Tk()
    root.withdraw()

    input_csv = (
        filedialog.askopenfilename(
            title=(
                "Select the MIBiG accession and NPAtlas SMILES CSV (train_stratefied_smiles.csv)"
            ),
            filetypes=[
                (
                    "CSV files",
                    "*.csv",
                ),
                (
                    "All files",
                    "*.*",
                ),
            ],
        )
    )

    if not input_csv:

        root.destroy()

        return None

    json_folder = (
        filedialog.askdirectory(
            title=(
                "Select the MIBiG 4.0 JSON folder"
            )
        )
    )

    if not json_folder:

        root.destroy()

        return None

    npatlas_tsv = (
        filedialog.askopenfilename(
            title=(
                "Select NPAtlas_as.tsv"
            ),
            filetypes=[
                (
                    "TSV files",
                    "*.tsv",
                ),
                (
                    "Text files",
                    "*.txt",
                ),
                (
                    "All files",
                    "*.*",
                ),
            ],
        )
    )

    if not npatlas_tsv:

        root.destroy()

        return None

    coconut_csv = (
        filedialog.askopenfilename(
            title=(
                "Select the COCONUT database CSV (coconut_csv-10_2026.csv)"
            ),
            filetypes=[
                (
                    "CSV files",
                    "*.csv",
                ),
                (
                    "All files",
                    "*.*",
                ),
            ],
        )
    )

    if not coconut_csv:

        root.destroy()

        return None

    output_csv = (
        filedialog.asksaveasfilename(
            title=(
                "Save output CSV as"
            ),
            defaultextension=".csv",
            initialfile=(
                "mibig_npatlas_coconut_smiles.csv"
            ),
            filetypes=[
                (
                    "CSV files",
                    "*.csv",
                ),
                (
                    "All files",
                    "*.*",
                ),
            ],
        )
    )

    root.destroy()

    if not output_csv:

        return None

    return (
        Path(input_csv),
        Path(json_folder),
        Path(npatlas_tsv),
        Path(coconut_csv),
        Path(output_csv),
    )

# Print results summary and messages.
def main():

    selected = choose_files()

    if selected is None:

        print("Cancelled.")

        return

    (
        input_csv,
        json_folder,
        npatlas_tsv,
        coconut_csv,
        output_csv,
    ) = selected

    accession_to_smiles = (
        load_input_csv(
            input_csv
        )
    )

    print(
        f"Found "
        f"{len(accession_to_smiles)} "
        f"relevant MIBiG accessions."
    )

    print(
        "\nLoading NPAtlas_as.tsv..."
    )

    npatlas_smiles_to_name = (
        load_npatlas_tsv(
            npatlas_tsv
        )
    )

    print(
        f"Loaded "
        f"{len(npatlas_smiles_to_name)} "
        f"exact SMILES/name mappings."
    )

    print(
        "\nLoading COCONUT CSV..."
    )

    coconut_name_to_smiles = (
        load_coconut_csv(
            coconut_csv
        )
    )

    print(
        f"Loaded "
        f"{len(coconut_name_to_smiles)} "
        f"COCONUT name/SMILES mappings."
    )

    all_rows = []
    missing_json = []

    for number, (
        accession,
        np_smiles,
    ) in enumerate(
        accession_to_smiles.items(),
        1,
    ):

        print(
            f"[{number}/"
            f"{len(accession_to_smiles)}] "
            f"{accession}"
        )

        json_path = find_json(
            json_folder,
            accession,
        )

        if json_path is None:

            print(
                f"    WARNING: "
                f"{accession}.json not found"
            )

            missing_json.append(
                accession
            )

            continue

        try:

            (
                _,
                mibig_compounds,
                dois,
                pmids,
            ) = read_mibig_json(
                json_path
            )

            all_rows.extend(
                compare_accession(
                    accession,
                    np_smiles,
                    mibig_compounds,
                    dois,
                    pmids,
                )
            )

        except Exception as exc:

            print(
                f"    ERROR: {exc}"
            )

    print(
        "\nFilling missing compound names "
        "from NPAtlas_as.tsv..."
    )

    names_filled = (
        fill_missing_compound_names(
            all_rows,
            npatlas_smiles_to_name,
        )
    )

    print(
        f"Compound names filled: "
        f"{names_filled}"
    )

    print(
        "\nGrouping rows by compound name "
        "(case-insensitive)..."
    )

    rows_before_grouping = len(
        all_rows
    )

    all_rows = (
        group_rows_by_compound_name(
            all_rows
        )
    )

    rows_after_grouping = len(
        all_rows
    )

    print(
        f"Rows before grouping: "
        f"{rows_before_grouping}"
    )

    print(
        f"Rows after grouping:  "
        f"{rows_after_grouping}"
    )

    print(
        f"Rows combined:        "
        f"{rows_before_grouping - rows_after_grouping}"
    )

    print(
        "\nMatching final compound names "
        "to COCONUT..."
    )

    coconut_matches = (
        add_coconut_smiles(
            all_rows,
            coconut_name_to_smiles,
        )
    )

    print(
        f"COCONUT matches: "
        f"{coconut_matches}"
    )

    write_output(
        all_rows,
        output_csv,
    )

    print(
        f"Rows written: "
        f"{len(all_rows)}"
    )

    print(
        f"Output: "
        f"{output_csv}"
    )

    print(
        f"COCONUT matches: "
        f"{coconut_matches}"
    )

    if missing_json:

        print(
            "\nMIBiG JSON files not found:"
        )

        for accession in missing_json:

            print(
                f"  {accession}"
            )

    root = tk.Tk()
    root.withdraw()

    root.destroy()


if __name__ == "__main__":
    main()