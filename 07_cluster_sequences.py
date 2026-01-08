#!/usr/bin/env python3
"""
07_cluster_sequences.py - Phase 3: Sequence clustering to reduce redundancy

This script reduces redundancy in the candidate sequence set using:
1. CD-HIT (if available on the system)
2. Pure Python alternative if CD-HIT is not installed

Input: data/processed/candidates_final_seqs.fasta
Output: data/processed/candidates_clustered.fasta
"""

import os
import subprocess
import shutil
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# =============================================================================
# CONFIGURATION
# =============================================================================
INPUT_FASTA = Path("data/processed/candidates_final_seqs.fasta")
OUTPUT_FASTA = Path("data/processed/candidates_clustered.fasta")
IDENTITY_THRESHOLD = 0.8  # 80% identity
WORD_SIZE = 5  # For CD-HIT with -c 0.8, -n 5 is recommended

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def check_cdhit_available():
    """Check if CD-HIT is available on the system."""
    return shutil.which("cd-hit") is not None


def run_cdhit(input_file, output_file, identity=0.9, word_size=5):
    """
    Run CD-HIT for sequence clustering.
    
    Parameters:
        input_file: Input FASTA file
        output_file: Output FASTA file (representatives)
        identity: Identity threshold (0.0 - 1.0)
        word_size: Word size for the algorithm
    """
    cmd = [
        "cd-hit",
        "-i", str(input_file),
        "-o", str(output_file),
        "-c", str(identity),
        "-n", str(word_size),
        "-M", "0",  # No memory limit
        "-T", "0",  # Use all available threads
        "-d", "0"   # Keep full description
    ]
    
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"CD-HIT failed:\n{result.stderr}")
    
    print(result.stdout)
    return True


def calculate_identity(seq1, seq2):
    """
    Calculate sequence identity between two sequences.
    Uses a simple approximate global alignment method.
    """
    # For sequences of very different lengths, identity will be low
    len1, len2 = len(seq1), len(seq2)
    if len1 == 0 or len2 == 0:
        return 0.0
    
    # Calculate length ratio
    length_ratio = min(len1, len2) / max(len1, len2)
    if length_ratio < 0.8:  # If they differ too much in length, they're not similar
        return 0.0
    
    # K-mer method to approximate identity (fast)
    k = 3
    if len1 < k or len2 < k:
        # For very short sequences, compare directly
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        return matches / max(len1, len2)
    
    # Create k-mer set
    kmers1 = set(seq1[i:i+k] for i in range(len(seq1) - k + 1))
    kmers2 = set(seq2[i:i+k] for i in range(len(seq2) - k + 1))
    
    if not kmers1 or not kmers2:
        return 0.0
    
    # Jaccard index approximates identity
    intersection = len(kmers1 & kmers2)
    union = len(kmers1 | kmers2)
    
    return intersection / union if union > 0 else 0.0


def python_clustering(sequences, identity_threshold=0.9):
    """
    Pure Python alternative for sequence clustering.
    Implements a greedy algorithm similar to CD-HIT.
    
    Parameters:
        sequences: List of tuples (id, sequence)
        identity_threshold: Identity threshold for grouping
        
    Returns:
        List of representative sequence IDs
    """
    print("Using pure Python clustering (alternative to CD-HIT)...")
    print(f"Identity threshold: {identity_threshold * 100}%")
    
    # Sort by length (longest to shortest) - CD-HIT strategy
    sorted_seqs = sorted(sequences, key=lambda x: len(x[1]), reverse=True)
    
    representatives = []  # List of (id, sequence) representatives
    clusters = defaultdict(list)  # cluster_id -> list of members
    
    total = len(sorted_seqs)
    print(f"Processing {total} sequences...")
    
    for i, (seq_id, seq) in enumerate(sorted_seqs):
        if (i + 1) % 500 == 0:
            print(f"  Progress: {i + 1}/{total} sequences processed, "
                  f"{len(representatives)} clusters formed")
        
        # Check if it belongs to an existing cluster
        found_cluster = False
        for rep_id, rep_seq in representatives:
            identity = calculate_identity(str(seq), str(rep_seq))
            if identity >= identity_threshold:
                clusters[rep_id].append(seq_id)
                found_cluster = True
                break
        
        # If it doesn't belong to any cluster, create a new one
        if not found_cluster:
            representatives.append((seq_id, seq))
            clusters[seq_id] = [seq_id]
    
    print(f"\nClustering completed:")
    print(f"  - Clusters formed: {len(representatives)}")
    
    # Cluster statistics
    cluster_sizes = [len(members) for members in clusters.values()]
    print(f"  - Average cluster size: {sum(cluster_sizes)/len(cluster_sizes):.2f}")
    print(f"  - Largest cluster: {max(cluster_sizes)} members")
    print(f"  - Singletons (clusters of 1): {sum(1 for s in cluster_sizes if s == 1)}")
    
    return [rep_id for rep_id, _ in representatives]


def main():
    """Main function of the script."""
    print("=" * 70)
    print("PHASE 3 - STEP 1: SEQUENCE CLUSTERING")
    print("=" * 70)
    
    # Verify input file
    if not INPUT_FASTA.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_FASTA}")
    
    # Load sequences
    print(f"\nLoading sequences from: {INPUT_FASTA}")
    sequences = list(SeqIO.parse(INPUT_FASTA, "fasta"))
    initial_count = len(sequences)
    print(f"Initial sequences: {initial_count}")
    
    # Check if CD-HIT is available
    use_cdhit = check_cdhit_available()
    
    if use_cdhit:
        print("\n✓ CD-HIT detected on the system")
        print("-" * 50)
        
        # Run CD-HIT
        run_cdhit(INPUT_FASTA, OUTPUT_FASTA, IDENTITY_THRESHOLD, WORD_SIZE)
        
    else:
        print("\n✗ CD-HIT not found on the system")
        print("  Using pure Python alternative...")
        print("-" * 50)
        
        # Prepare data for Python clustering
        seq_data = [(record.id, record.seq) for record in sequences]
        
        # Run clustering
        representative_ids = python_clustering(seq_data, IDENTITY_THRESHOLD)
        
        # Create dictionary for fast lookup
        seq_dict = {record.id: record for record in sequences}
        
        # Write representatives
        print(f"\nSaving representative sequences to: {OUTPUT_FASTA}")
        representatives = [seq_dict[seq_id] for seq_id in representative_ids]
        SeqIO.write(representatives, OUTPUT_FASTA, "fasta")
    
    # Count final sequences
    final_sequences = list(SeqIO.parse(OUTPUT_FASTA, "fasta"))
    final_count = len(final_sequences)
    
    # Summary
    print("\n" + "=" * 70)
    print("CLUSTERING SUMMARY")
    print("=" * 70)
    print(f"Sequences before clustering:  {initial_count:,}")
    print(f"Sequences after clustering: {final_count:,}")
    print(f"Reduction: {initial_count - final_count:,} sequences ({(1 - final_count/initial_count)*100:.1f}%)")
    print(f"Identity threshold used: {IDENTITY_THRESHOLD * 100}%")
    print(f"\nOutput file: {OUTPUT_FASTA}")
    
    # Warning if there are still many sequences
    if final_count > 500:
        print(f"\n⚠️  WARNING: You still have {final_count} sequences.")
        print("   For 3D prediction, consider reducing the identity threshold:")
        print("   - Edit IDENTITY_THRESHOLD = 0.8 (80%) or lower")
        print("   - Or run: cd-hit -i input.fasta -o output.fasta -c 0.8")
    else:
        print(f"\n✓ Number of sequences ({final_count}) is manageable for 3D prediction")
    
    print("\n" + "=" * 70)
    print("Clustering completed successfully")
    print("=" * 70)


if __name__ == "__main__":
    main()
