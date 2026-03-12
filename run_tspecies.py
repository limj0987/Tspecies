#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import subprocess
import os
import sys
from Bio import SeqIO

def parse_args():
    """Parse command line arguments provided by the user."""
    parser = argparse.ArgumentParser(description="Tspecies: Automatic pipeline from Genomes to Ne and Divergence Time")
    
    # Input files
    parser.add_argument('--cds1', required=True, help='Path to CDS fasta file of Species A')
    parser.add_argument('--cds2', required=True, help='Path to CDS fasta file of Species B')
    parser.add_argument('--pep1', required=True, help='Path to Protein fasta file of Species A')
    parser.add_argument('--pep2', required=True, help='Path to Protein fasta file of Species B')
    
    # Biological parameters for Tspecies
    parser.add_argument('--mu', type=float, required=True, help='Neutral mutation rate (e.g., 1e-8)')
    parser.add_argument('--g', type=float, default=1.0, help='Generation time in years (default: 1)')
    
    # Computational parameters
    parser.add_argument('--out', default='tspecies_output', help='Output directory name')
    parser.add_argument('--threads', type=int, default=4, help='Number of CPU threads to use')
    
    return parser.parse_args()

def run_blast_rbh(pep1, pep2, out_dir, threads):
    """Step 1: Run bidirectional BLASTP to find Reciprocal Best Hits (RBH)."""
    print("--- Step 1: Running Bidirectional BLASTP (RBH) ---")
    
    # Create BLAST databases
    subprocess.run(f"makeblastdb -in {pep1} -dbtype prot -out {out_dir}/db_A", shell=True, check=True)
    subprocess.run(f"makeblastdb -in {pep2} -dbtype prot -out {out_dir}/db_B", shell=True, check=True)
    
    # Run BLAST A vs B
    print("Running BLAST: Species A vs Species B...")
    blast_A_vs_B = f"blastp -query {pep1} -db {out_dir}/db_B -out {out_dir}/A_vs_B.txt -outfmt 6 -max_target_seqs 1 -num_threads {threads}"
    subprocess.run(blast_A_vs_B, shell=True, check=True)
    
    # Run BLAST B vs A
    print("Running BLAST: Species B vs Species A...")
    blast_B_vs_A = f"blastp -query {pep2} -db {out_dir}/db_A -out {out_dir}/B_vs_A.txt -outfmt 6 -max_target_seqs 1 -num_threads {threads}"
    subprocess.run(blast_B_vs_A, shell=True, check=True)
    
    # Parse BLAST results to find Reciprocal Best Hits
    best_A_to_B = {}
    with open(f"{out_dir}/A_vs_B.txt") as f:
        for line in f:
            cols = line.strip().split('\t')
            qseqid, sseqid = cols[0], cols[1]
            if qseqid not in best_A_to_B:
                best_A_to_B[qseqid] = sseqid
                
    best_B_to_A = {}
    with open(f"{out_dir}/B_vs_A.txt") as f:
        for line in f:
            cols = line.strip().split('\t')
            qseqid, sseqid = cols[0], cols[1]
            if qseqid not in best_B_to_A:
                best_B_to_A[qseqid] = sseqid
                
    # Find intersections (RBH)
    rbh_pairs = []
    for gene_A, gene_B in best_A_to_B.items():
        if best_B_to_A.get(gene_B) == gene_A:
            rbh_pairs.append((gene_A, gene_B))
            
    print(f"Found {len(rbh_pairs)} Reciprocal Best Hit (RBH) ortholog pairs.")
    return rbh_pairs

def align_and_convert_to_axt(rbh_pairs, cds1_path, cds2_path, out_dir):
    """Step 2: Extract CDS, run MAFFT alignment, and convert to AXT format."""
    print("--- Step 2: Aligning sequences and generating AXT format ---")
    
    # Load sequences into memory
    print("Loading CDS sequences...")
    cds1_dict = SeqIO.to_dict(SeqIO.parse(cds1_path, "fasta"))
    cds2_dict = SeqIO.to_dict(SeqIO.parse(cds2_path, "fasta"))
    
    axt_file_path = f"{out_dir}/input.axt"
    temp_fasta = f"{out_dir}/temp_pair.fasta"
    
    valid_pairs_count = 0
    
    with open(axt_file_path, "w") as axt_out:
        for i, (gene_A, gene_B) in enumerate(rbh_pairs):
            # Print progress
            if (i + 1) % 500 == 0:
                print(f"Processed {i + 1} pairs...")
                
            if gene_A in cds1_dict and gene_B in cds2_dict:
                seq_A = cds1_dict[gene_A]
                seq_B = cds2_dict[gene_B]
                
                # Write the pair to a temporary fasta file
                SeqIO.write([seq_A, seq_B], temp_fasta, "fasta")
                
                # Run MAFFT for pairwise alignment
                mafft_cmd = f"mafft --quiet --auto {temp_fasta}"
                result = subprocess.run(mafft_cmd, shell=True, capture_output=True, text=True)
                
                # Parse MAFFT output (which is in FASTA format)
                with open(f"{out_dir}/temp_aligned.fasta", "w") as temp_out:
                    temp_out.write(result.stdout)
                
                aligned_seqs = list(SeqIO.parse(f"{out_dir}/temp_aligned.fasta", "fasta"))
                
                if len(aligned_seqs) == 2:
                    aln_A = str(aligned_seqs[0].seq)
                    aln_B = str(aligned_seqs[1].seq)
                    
                    # Write to AXT format
                    # Format: 
                    # GeneNameA&GeneNameB
                    # ALIGNED_SEQ_A
                    # ALIGNED_SEQ_B
                    # <empty line>
                    axt_out.write(f"{gene_A}&{gene_B}\n{aln_A}\n{aln_B}\n\n")
                    valid_pairs_count += 1

    # Clean up temporary files
    if os.path.exists(temp_fasta): os.remove(temp_fasta)
    if os.path.exists(f"{out_dir}/temp_aligned.fasta"): os.remove(f"{out_dir}/temp_aligned.fasta")
    
    print(f"Successfully generated AXT file for {valid_pairs_count} pairs.")
    return axt_file_path

def run_kaks_and_tspecies(axt_file, mu, g, out_dir):
    """Step 3 & 4: Calculate Ks and run Tspecies R script."""
    print("--- Step 3: Running KaKs_Calculator ---")
    ks_results = f"{out_dir}/kaks_results.txt"
    # Using YN method as a standard model
    kaks_cmd = f"KaKs_Calculator -i {axt_file} -o {ks_results} -m YN"
    subprocess.run(kaks_cmd, shell=True, check=True)
    
    print("--- Step 4: Running Tspecies prediction in R ---")
    r_script_path = f"{out_dir}/run_tspecies_model.R"
    final_csv = f"{out_dir}/tspecies_final_report.csv"
    
    # Generate R script dynamically
    r_code = f"""
# Load the Tspecies package
library(Tspecies)

# Read KaKs_Calculator output
cat("Loading Ks results...\\n")
kaks_data <- read.table("{ks_results}", header=TRUE, sep="\\t", stringsAsFactors=FALSE)

# Filter out invalid Ks values (NA or extremely high values indicating saturation)
valid_ks <- kaks_data$Ks[!is.na(kaks_data$Ks) & kaks_data$Ks < 5]

if(length(valid_ks) < 50) {{
    warning("Too few valid Ks values found. Prediction may be inaccurate.")
}}

# Run Tspecies model
cat("Running GAM prediction...\\n")
# Capture output to file or console
sink("{out_dir}/tspecies_console_log.txt")
Tspecies(Ks = valid_ks, miu = {mu}, g = {g})
sink()

cat("Analysis complete! Check {out_dir}/tspecies_console_log.txt for detailed logs.\\n")
"""
    
    with open(r_script_path, "w") as f:
        f.write(r_code)
        
    subprocess.run(f"Rscript {r_script_path}", shell=True, check=True)
    print(f"--- Pipeline Finished! Results saved in '{out_dir}' directory ---")
    print(f"Please check '{out_dir}/tspecies_console_log.txt' for the Ne and Divergence Time results.")

def main():
    args = parse_args()
    
    # Create output directory
    if not os.path.exists(args.out):
        os.makedirs(args.out)
        
    # Run Pipeline
    rbh_pairs = run_blast_rbh(args.pep1, args.pep2, args.out, args.threads)
    
    if len(rbh_pairs) == 0:
        print("Error: No orthologs found. Check your input fasta files.")
        sys.exit(1)
        
    axt_file = align_and_convert_to_axt(rbh_pairs, args.cds1, args.cds2, args.out)
    
    run_kaks_and_tspecies(axt_file, args.mu, args.g, args.out)

if __name__ == "__main__":
    main()