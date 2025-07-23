#!/usr/bin/env python3
# probe_picker_gui.py - GUI-based Top Probe Selection Tool

import pandas as pd
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import threading


class ProbePickerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("FISH-RT Probe Picker - Top Probes per Gene")
        self.root.geometry("800x600")

        # Variables
        self.input_file = tk.StringVar()
        self.output_prefix = tk.StringVar(value="FISH_RT_probes_TOP10")
        self.probes_per_gene = tk.IntVar(value=10)
        self.results_text = tk.StringVar()

        self.setup_gui()

    def setup_gui(self):
        """Setup the GUI components"""

        # Main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Title
        title_label = ttk.Label(
            main_frame, text="FISH-RT Probe Picker", font=("Arial", 16, "bold")
        )
        title_label.grid(row=0, column=0, columnspan=3, pady=(0, 20))

        # File selection section
        file_frame = ttk.LabelFrame(
            main_frame, text="Input File Selection", padding="10"
        )
        file_frame.grid(
            row=1, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 10)
        )

        ttk.Label(file_frame, text="CSV File:").grid(row=0, column=0, sticky=tk.W)
        ttk.Entry(
            file_frame, textvariable=self.input_file, width=60, state="readonly"
        ).grid(row=0, column=1, padx=(10, 5), sticky=(tk.W, tk.E))
        ttk.Button(file_frame, text="Browse...", command=self.browse_file).grid(
            row=0, column=2, padx=(5, 0)
        )

        # Configuration section
        config_frame = ttk.LabelFrame(
            main_frame, text="Selection Parameters", padding="10"
        )
        config_frame.grid(
            row=2, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 10)
        )

        ttk.Label(config_frame, text="Probes per gene:").grid(
            row=0, column=0, sticky=tk.W
        )
        probe_spinbox = ttk.Spinbox(
            config_frame, from_=1, to=50, width=10, textvariable=self.probes_per_gene
        )
        probe_spinbox.grid(row=0, column=1, padx=(10, 0), sticky=tk.W)

        ttk.Label(config_frame, text="Output prefix:").grid(
            row=1, column=0, sticky=tk.W, pady=(10, 0)
        )
        ttk.Entry(config_frame, textvariable=self.output_prefix, width=30).grid(
            row=1, column=1, padx=(10, 0), sticky=tk.W, pady=(10, 0)
        )

        # Sorting criteria info
        info_frame = ttk.LabelFrame(main_frame, text="Selection Criteria", padding="10")
        info_frame.grid(
            row=3, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 10)
        )

        criteria_text = """Probes are sorted by:
1. SNPs_Covered_Count (descending) - Higher SNP coverage first
2. NbOfPNAS (descending) - Better PNAS score first
3. dGScore (descending) - Better thermodynamic score first

Top N probes per gene are then selected."""

        ttk.Label(info_frame, text=criteria_text, justify=tk.LEFT).grid(
            row=0, column=0, sticky=tk.W
        )

        # Action buttons
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=4, column=0, columnspan=3, pady=(10, 0))

        self.process_button = ttk.Button(
            button_frame,
            text="Select Top Probes",
            command=self.process_probes,
            style="Accent.TButton",
        )
        self.process_button.pack(side=tk.LEFT, padx=(0, 10))

        ttk.Button(button_frame, text="Clear Results", command=self.clear_results).pack(
            side=tk.LEFT
        )

        # Results section
        results_frame = ttk.LabelFrame(main_frame, text="Results", padding="10")
        results_frame.grid(
            row=5, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=(10, 0)
        )

        # Results text area with scrollbar
        text_frame = ttk.Frame(results_frame)
        text_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        self.results_textbox = tk.Text(text_frame, height=15, width=80, wrap=tk.WORD)
        scrollbar = ttk.Scrollbar(
            text_frame, orient=tk.VERTICAL, command=self.results_textbox.yview
        )
        self.results_textbox.configure(yscrollcommand=scrollbar.set)

        self.results_textbox.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))

        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(5, weight=1)
        file_frame.columnconfigure(1, weight=1)
        results_frame.columnconfigure(0, weight=1)
        results_frame.rowconfigure(0, weight=1)
        text_frame.columnconfigure(0, weight=1)
        text_frame.rowconfigure(0, weight=1)

    def browse_file(self):
        """Open file dialog to select CSV file"""
        filename = filedialog.askopenfilename(
            title="Select FISH-RT Probes CSV File",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if filename:
            self.input_file.set(filename)

    def process_probes(self):
        """Process the probe selection in a separate thread"""
        if not self.input_file.get():
            messagebox.showerror("Error", "Please select a CSV file first.")
            return

        # Disable button during processing
        self.process_button.configure(state="disabled")
        self.results_textbox.delete(1.0, tk.END)
        self.results_textbox.insert(tk.END, "Processing... Please wait...\n")
        self.root.update()

        # Run processing in separate thread to avoid freezing GUI
        thread = threading.Thread(target=self._process_probes_thread)
        thread.daemon = True
        thread.start()

    def _process_probes_thread(self):
        """Actual probe processing logic"""
        try:
            self._run_probe_selection()
        except Exception as e:
            self.root.after(0, lambda: self._show_error(str(e)))
        finally:
            self.root.after(0, lambda: self.process_button.configure(state="normal"))

    def _run_probe_selection(self):
        """Run the probe selection process"""
        input_file = self.input_file.get()
        output_prefix = self.output_prefix.get()
        probes_per_gene = self.probes_per_gene.get()

        # Read CSV file
        self._update_results("📊 Reading CSV file...\n")
        df = pd.read_csv(input_file)

        self._update_results(f"   Total probes: {len(df)}\n")
        self._update_results(f"   Genes found: {df['GeneName'].nunique()}\n\n")

        # Sort probes by selection criteria
        self._update_results("🔄 Sorting probes by SNP coverage and PNAS score...\n")
        df_sorted = df.sort_values(
            ["SNPs_Covered_Count", "NbOfPNAS", "dGScore"],
            ascending=[False, False, False],
        )

        # Select top probes per gene
        self._update_results(f"🎯 Selecting top {probes_per_gene} probes per gene...\n")
        top_probes = (
            df_sorted.groupby("GeneName").head(probes_per_gene).reset_index(drop=True)
        )

        # Report selection statistics
        gene_counts = top_probes["GeneName"].value_counts().sort_index()
        self._update_results("\n📈 Selected probes per gene:\n")
        for gene, count in gene_counts.items():
            avg_snps = top_probes[top_probes["GeneName"] == gene][
                "SNPs_Covered_Count"
            ].mean()
            avg_pnas = top_probes[top_probes["GeneName"] == gene]["NbOfPNAS"].mean()
            self._update_results(
                f"   {gene}: {count} probes (avg SNPs: {avg_snps:.1f}, avg PNAS: {avg_pnas:.1f}/5)\n"
            )

        self._update_results(
            f"\n✅ Total selected: {len(top_probes)} probes from {top_probes['GeneName'].nunique()} genes\n\n"
        )

        # Save filtered CSV
        input_path = Path(input_file)
        output_dir = input_path.parent
        csv_output = output_dir / f"{output_prefix}.csv"

        # Sort final output by GeneName, then SNPs_Covered_Count (desc), then NbOfPNAS (desc)
        top_probes = top_probes.sort_values(
            ["GeneName", "SNPs_Covered_Count", "NbOfPNAS"],
            ascending=[True, False, False],
        )

        top_probes.to_csv(csv_output, index=False)
        self._update_results(f"💾 Saved CSV: {csv_output}\n")

        # Generate FASTA file
        self._update_results(f"🧬 Generating FASTA file...\n")
        fasta_output = output_dir / f"{output_prefix}.fasta"
        self._generate_fasta_from_csv(top_probes, fasta_output)
        self._update_results(
            f"💾 Saved FASTA: {fasta_output} ({len(top_probes)} sequences)\n\n"
        )

        # Generate summary statistics
        self._update_results("📊 SUMMARY STATISTICS:\n")
        self._update_results(f"   Input probes: {len(df)}\n")
        self._update_results(f"   Selected probes: {len(top_probes)}\n")
        self._update_results(
            f"   Reduction: {(1 - len(top_probes)/len(df))*100:.1f}%\n"
        )
        self._update_results(
            f"   Average SNPs per probe: {top_probes['SNPs_Covered_Count'].mean():.1f}\n"
        )
        self._update_results(
            f"   Average PNAS score: {top_probes['NbOfPNAS'].mean():.1f}/5\n"
        )

        # SNP coverage distribution
        snp_dist = top_probes["SNPs_Covered_Count"].value_counts().sort_index()
        self._update_results(f"\n   SNP coverage distribution:\n")
        for snps, count in snp_dist.items():
            self._update_results(f"     {snps} SNPs: {count} probes\n")

        self._update_results(f"\n🎉 Probe selection completed!\n")
        self._update_results(
            f"   Output files: {csv_output.name}, {fasta_output.name}\n"
        )

        # Show completion message
        self.root.after(
            0,
            lambda: messagebox.showinfo(
                "Success",
                f"Probe selection completed!\n\nSelected {len(top_probes)} probes from {top_probes['GeneName'].nunique()} genes.\n\nFiles saved in: {output_dir}",
            ),
        )

    def _generate_fasta_from_csv(self, df, output_file):
        """Generate FASTA file from CSV data"""
        records = []

        for idx, row in df.iterrows():
            gene_name = row["GeneName"]

            # Use probe sequence only (no RTBC for BLAST)
            sequence = row["Seq"]
            seq_type = "probe_only"

            # Create descriptive FASTA header
            probe_size = row.get("ProbeSize", len(row["Seq"]))
            snp_count = row.get("SNPs_Covered_Count", 0)
            region_type = row.get("RegionType", "unknown")
            pnas_score = row.get("NbOfPNAS", 0)

            description = f"{gene_name} | {seq_type} | size:{probe_size}nt | SNPs:{snp_count} | region:{region_type} | PNAS:{pnas_score}/5"

            record = SeqRecord(
                Seq(sequence), id=f"{gene_name}_probe_{idx}", description=description
            )
            records.append(record)

        SeqIO.write(records, output_file, "fasta")

    def _update_results(self, text):
        """Update results text box in a thread-safe way"""
        self.root.after(0, lambda: self._append_text(text))

    def _append_text(self, text):
        """Append text to results textbox"""
        self.results_textbox.insert(tk.END, text)
        self.results_textbox.see(tk.END)
        self.root.update()

    def _show_error(self, error_msg):
        """Show error message"""
        self.results_textbox.delete(1.0, tk.END)
        self.results_textbox.insert(tk.END, f"❌ Error: {error_msg}\n")
        messagebox.showerror("Error", f"An error occurred:\n\n{error_msg}")

    def clear_results(self):
        """Clear the results text box"""
        self.results_textbox.delete(1.0, tk.END)


def main():
    root = tk.Tk()
    app = ProbePickerGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
