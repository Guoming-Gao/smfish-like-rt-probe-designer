#!/usr/bin/env python3
# generate_fasta_from_final_csv.py

import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import tkinter as tk
from tkinter import filedialog, messagebox


def read_input_file(input_file):
    """Read either xlsx or csv file and return DataFrame"""
    file_path = Path(input_file)

    if file_path.suffix.lower() == ".xlsx":
        print(f"📊 Reading Excel file: {input_file}")
        return pd.read_excel(input_file)
    elif file_path.suffix.lower() == ".csv":
        print(f"📊 Reading CSV file: {input_file}")
        return pd.read_csv(input_file)
    else:
        raise ValueError(
            f"Unsupported file format: {file_path.suffix}. Please use .xlsx or .csv files."
        )


def xlsx_to_fasta(input_file, output_fasta=None):
    """Convert xlsx or csv file to FASTA format"""

    # Read input file (xlsx or csv)
    df = read_input_file(input_file)
    print(f" Total records: {len(df)}")

    # Set output filename if not provided
    if output_fasta is None:
        input_path = Path(input_file)
        output_fasta = input_path.with_suffix(".fasta")

    # Generate FASTA records
    print(f"🧬 Generating FASTA file...")
    records = []

    for idx, row in df.iterrows():
        # Use ProbeName for the header/ID
        probe_name = row["ProbeName"]

        # Use RTBC_5Prime_Sequence for the sequence
        sequence = row["RTBC_5Prime_Sequence"]

        # Create descriptive FASTA header with additional info
        # You can modify these fields based on what's available in your file
        probe_size = len(sequence)

        # Add other info if columns exist (with fallbacks)
        gene_name = row.get("GeneName", "unknown")
        snp_count = row.get("SNPs_Covered_Count", 0)
        region_type = row.get("RegionType", "unknown")
        pnas_score = row.get("NbOfPNAS", 0)

        # Create descriptive header
        description = f"{gene_name} | RTBC_sequence | size:{probe_size}nt | SNPs:{snp_count} | region:{region_type} | PNAS:{pnas_score}/5"

        # Create SeqRecord
        record = SeqRecord(Seq(sequence), id=f"{probe_name}", description=description)
        records.append(record)

    # Write FASTA file
    SeqIO.write(records, output_fasta, "fasta")
    print(f"💾 Saved FASTA: {output_fasta} ({len(records)} sequences)")

    return output_fasta


def select_input_file():
    """Open file dialog to select input file"""
    # Create root window (hidden)
    root = tk.Tk()
    root.withdraw()  # Hide the main window

    # Configure file dialog
    file_types = [
        ("Excel and CSV files", "*.xlsx *.csv"),
        ("Excel files", "*.xlsx"),
        ("CSV files", "*.csv"),
        ("All files", "*.*"),
    ]

    # Open file dialog
    input_file = filedialog.askopenfilename(
        title="Select input file (xlsx or csv)", filetypes=file_types, initialdir="."
    )

    root.destroy()  # Clean up

    return input_file


def main():
    try:
        # If command line argument provided, use it; otherwise show file picker
        if len(sys.argv) >= 2:
            input_file = sys.argv[1]
            output_file = sys.argv[2] if len(sys.argv) > 2 else None
        else:
            # Show file picker
            input_file = select_input_file()

            if not input_file:  # User cancelled
                print("❌ No file selected. Exiting.")
                sys.exit(1)

            output_file = None

        # Convert file
        result = xlsx_to_fasta(input_file, output_file)

        # Show success message
        print(f"🎉 Conversion completed successfully!")
        print(f" Output: {result}")

        # Optional: Show success popup
        try:
            root = tk.Tk()
            root.withdraw()
            messagebox.showinfo(
                "Success",
                f"Conversion completed successfully!\n\nOutput saved to:\n{result}",
            )
            root.destroy()
        except:
            pass  # Skip popup if there's an issue with GUI

    except Exception as e:
        print(f"❌ Error: {e}")

        # Show error popup
        try:
            root = tk.Tk()
            root.withdraw()
            messagebox.showerror("Error", f"An error occurred:\n\n{str(e)}")
            root.destroy()
        except:
            pass  # Skip popup if there's an issue with GUI

        sys.exit(1)


if __name__ == "__main__":
    main()
