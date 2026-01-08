#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
04_filter_homology.py - Homology Filtering using Local Alignment
================================================================================

Description:
    This script performs filtering of candidate sequences based on
    homology with the reference enzyme IsPETase. It uses local alignment
    (Smith-Waterman) to detect conserved regions even in proteins
    that have diverged evolutionarily.

    Filtering strategy:
    - Identity > 20%: Ensures remote homology (likely similar function)
    - Identity < 90%: Ensures novelty (not the same known enzyme)

    Scientific rationale:
    Local alignment is preferable to global alignment for divergent proteins
    because it detects conserved domains without penalizing variable regions
    (loops, insertions, deletions).

Date: December 2025
Project: Deep-PETase-Mining - Phase 2

Usage:
    python 04_filter_homology.py

Dependencies:
    - Biopython (Bio.Align.PairwiseAligner)
    - Pandas
    - tqdm (optional)

Inputs:
    - data/references/ispetase_ref.fasta
    - data/raw/raw_candidates.fasta

Outputs:
    - data/processed/candidates_homology.fasta
    - results/logs/homology_scores.csv
================================================================================
"""

from pathlib import Path
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
from datetime import datetime
import logging
import sys

# Biopython
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import PairwiseAligner, substitution_matrices
except ImportError:
    print("❌ Error: Biopython is not installed.")
    print("   Run: pip install biopython")
    sys.exit(1)

# Pandas
try:
    import pandas as pd
except ImportError:
    print("❌ Error: Pandas is not installed.")
    print("   Run: pip install pandas")
    sys.exit(1)

# tqdm (opcional)
try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False


# ============================================================================
# CONFIGURATION
# ============================================================================

# Project paths
BASE_DIR: Path = Path(__file__).parent.resolve()
REFERENCE_FILE: Path = BASE_DIR / "data" / "references" / "ispetase_ref.fasta"
CANDIDATES_FILE: Path = BASE_DIR / "data" / "raw" / "raw_candidates.fasta"
OUTPUT_FASTA: Path = BASE_DIR / "data" / "processed" / "candidates_homology.fasta"
OUTPUT_CSV: Path = BASE_DIR / "results" / "logs" / "homology_scores.csv"

# Filtering parameters
MIN_IDENTITY: float = 30.0   # Minimum identity to consider homology (raised from 20%)
MAX_IDENTITY: float = 90.0   # Maximum identity to ensure novelty

# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# ============================================================================
# DATA STRUCTURES
# ============================================================================

@dataclass
class AlignmentResult:
    """
    Stores the results of an alignment.
    
    Note for students:
        We use dataclass to create clean data structures.
        It's like a struct in C or a record in Pascal, but Pythonic.
    """
    candidate_id: str
    candidate_description: str
    alignment_score: float
    identity_percent: float
    aligned_length: int
    candidate_length: int
    reference_length: int
    passed_filter: bool


# ============================================================================
# LOADING FUNCTIONS
# ============================================================================

def load_reference(filepath: Path) -> SeqRecord:
    """
    Loads the reference sequence from a FASTA file.
    
    Args:
        filepath (Path): Path to the reference FASTA file.
    
    Returns:
        SeqRecord: Reference sequence record.
    
    Raises:
        FileNotFoundError: If the file does not exist.
    """
    if not filepath.exists():
        logger.error(f"❌ Reference file not found: {filepath}")
        logger.error("   Run first: python 03_get_reference.py")
        raise FileNotFoundError(f"Not found: {filepath}")
    
    record = SeqIO.read(filepath, "fasta")
    logger.info(f"✅ Reference loaded: {record.id} ({len(record.seq)} aa)")
    
    return record


def load_candidates(filepath: Path) -> List[SeqRecord]:
    """
    Loads candidate sequences from a FASTA file.
    
    Args:
        filepath (Path): Path to the candidates FASTA file.
    
    Returns:
        List[SeqRecord]: List of sequence records.
    """
    if not filepath.exists():
        logger.error(f"❌ Candidates file not found: {filepath}")
        logger.error("   Run first: python 01_mining_ncbi.py")
        raise FileNotFoundError(f"Not found: {filepath}")
    
    records = list(SeqIO.parse(filepath, "fasta"))
    logger.info(f"✅ Candidates loaded: {len(records)} sequences")
    
    return records


# ============================================================================
# ALIGNMENT FUNCTIONS
# ============================================================================

def create_aligner() -> PairwiseAligner:
    """
    Creates and configures the aligner for local alignment.
    
    Returns:
        PairwiseAligner: Aligner configured for proteins.
    
    Note for students:
        LOCAL alignment (Smith-Waterman) searches for the best region
        of similarity between two sequences, ignoring the ends.
        It is ideal for detecting conserved domains.
        
        BLOSUM62 is a substitution matrix derived from blocks
        of aligned sequences. The 62 means that sequences with
        ~62% identity were used to build it.
    """
    aligner = PairwiseAligner()
    
    # Configure local mode (Smith-Waterman)
    aligner.mode = 'local'
    
    # Try to use BLOSUM62 matrix (standard for proteins)
    try:
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        logger.info("📊 Using BLOSUM62 substitution matrix")
    except Exception as e:
        logger.warning(f"⚠️ Could not load BLOSUM62, using defaults: {e}")
        # Use default parameters if BLOSUM62 is not available
        aligner.match_score = 2
        aligner.mismatch_score = -1
    
    # Gap penalties
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    
    return aligner


def calculate_identity(seq1: str, seq2: str, alignment) -> Tuple[float, int]:
    """
    Calculates the identity percentage from an alignment.
    
    Args:
        seq1 (str): First sequence.
        seq2 (str): Second sequence.
        alignment: Biopython alignment object.
    
    Returns:
        Tuple[float, int]: (identity percentage, aligned length)
    
    Note for students:
        Identity is calculated as:
        (identical residues / alignment length) * 100
        
        In local alignment, the length is the aligned region,
        not the total length of the sequences.
    """
    # Get the aligned sequences
    aligned_seq1 = alignment[0]
    aligned_seq2 = alignment[1]
    
    # Count identical residues (excluding gaps)
    identical = 0
    aligned_positions = 0
    
    for a, b in zip(str(aligned_seq1), str(aligned_seq2)):
        if a != '-' and b != '-':
            aligned_positions += 1
            if a == b:
                identical += 1
    
    # Calculate percentage
    if aligned_positions > 0:
        identity_percent = (identical / aligned_positions) * 100
    else:
        identity_percent = 0.0
    
    return identity_percent, aligned_positions


def align_sequences(
    reference: SeqRecord,
    candidate: SeqRecord,
    aligner: PairwiseAligner
) -> AlignmentResult:
    """
    Performs alignment between the reference and a candidate.
    
    Args:
        reference (SeqRecord): Reference sequence.
        candidate (SeqRecord): Candidate sequence.
        aligner (PairwiseAligner): Configured aligner.
    
    Returns:
        AlignmentResult: Alignment result with metrics.
    """
    ref_seq = str(reference.seq)
    cand_seq = str(candidate.seq)
    
    try:
        # Perform alignment
        alignments = aligner.align(ref_seq, cand_seq)
        
        # Take the best alignment
        if len(alignments) > 0:
            best_alignment = alignments[0]
            score = best_alignment.score
            
            # Calculate identity
            identity, aligned_length = calculate_identity(
                ref_seq, cand_seq, best_alignment
            )
        else:
            score = 0.0
            identity = 0.0
            aligned_length = 0
        
        # Determine if it passes the filter
        passed = MIN_IDENTITY <= identity <= MAX_IDENTITY
        
        return AlignmentResult(
            candidate_id=candidate.id,
            candidate_description=candidate.description,
            alignment_score=score,
            identity_percent=identity,
            aligned_length=aligned_length,
            candidate_length=len(cand_seq),
            reference_length=len(ref_seq),
            passed_filter=passed
        )
        
    except Exception as e:
        logger.warning(f"⚠️ Error aligning {candidate.id}: {e}")
        return AlignmentResult(
            candidate_id=candidate.id,
            candidate_description=candidate.description,
            alignment_score=0.0,
            identity_percent=0.0,
            aligned_length=0,
            candidate_length=len(cand_seq),
            reference_length=len(ref_seq),
            passed_filter=False
        )


# ============================================================================
# FILTERING AND SAVING FUNCTIONS
# ============================================================================

def filter_candidates(
    reference: SeqRecord,
    candidates: List[SeqRecord],
    aligner: PairwiseAligner
) -> Tuple[List[SeqRecord], List[AlignmentResult]]:
    """
    Filters candidates based on homology with the reference.
    
    Args:
        reference (SeqRecord): Reference sequence.
        candidates (List[SeqRecord]): List of candidates.
        aligner (PairwiseAligner): Configured aligner.
    
    Returns:
        Tuple: (sequences that passed the filter, all results)
    """
    passed_sequences: List[SeqRecord] = []
    all_results: List[AlignmentResult] = []
    
    # Configure iterator with or without progress bar
    if TQDM_AVAILABLE:
        iterator = tqdm(candidates, desc="Aligning", unit="seq")
    else:
        iterator = candidates
        logger.info("Processing alignments (this may take a few minutes)...")
    
    for candidate in iterator:
        result = align_sequences(reference, candidate, aligner)
        all_results.append(result)
        
        if result.passed_filter:
            passed_sequences.append(candidate)
    
    return passed_sequences, all_results


def save_filtered_fasta(sequences: List[SeqRecord], output_path: Path) -> int:
    """
    Saves filtered sequences in FASTA format.
    
    Args:
        sequences (List[SeqRecord]): Sequences to save.
        output_path (Path): Output file path.
    
    Returns:
        int: Number of sequences saved.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    count = SeqIO.write(sequences, output_path, "fasta")
    logger.info(f"💾 {count} sequences saved to: {output_path}")
    
    return count


def save_results_csv(results: List[AlignmentResult], output_path: Path) -> None:
    """
    Saves alignment results in CSV format.
    
    Args:
        results (List[AlignmentResult]): List of results.
        output_path (Path): CSV file path.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Convert to DataFrame
    data = [{
        "ID": r.candidate_id,
        "Description": r.candidate_description[:100],  # Truncate long descriptions
        "Score": round(r.alignment_score, 2),
        "Identity_Percent": round(r.identity_percent, 2),
        "Aligned_Length": r.aligned_length,
        "Candidate_Length": r.candidate_length,
        "Passed_Filter": r.passed_filter
    } for r in results]
    
    df = pd.DataFrame(data)
    df = df.sort_values("Identity_Percent", ascending=False)
    df.to_csv(output_path, index=False)
    
    logger.info(f"💾 Results saved to: {output_path}")


def print_summary(
    total: int,
    passed: int,
    results: List[AlignmentResult]
) -> None:
    """
    Prints a summary of the filtering process.
    
    Args:
        total (int): Total number of candidates.
        passed (int): Number that passed the filter.
        results (List[AlignmentResult]): All results.
    """
    # Calculate statistics
    identities = [r.identity_percent for r in results if r.identity_percent > 0]
    
    print("\n" + "="*70)
    print("📊 HOMOLOGY FILTERING SUMMARY")
    print("="*70)
    
    print(f"\n📈 General statistics:")
    print(f"   • Sequences analyzed:       {total:,}")
    print(f"   • Sequences that passed:    {passed:,}")
    print(f"   • Sequences discarded:      {total - passed:,}")
    print(f"   • Survival rate:            {(passed/total)*100:.1f}%")
    
    if identities:
        print(f"\n📊 Identity distribution:")
        print(f"   • Minimum:  {min(identities):.1f}%")
        print(f"   • Maximum:  {max(identities):.1f}%")
        print(f"   • Average:  {sum(identities)/len(identities):.1f}%")
    
    print(f"\n🎯 Filtering criteria applied:")
    print(f"   • Minimum identity: {MIN_IDENTITY}% (remote homology)")
    print(f"   • Maximum identity: {MAX_IDENTITY}% (ensure novelty)")
    
    # Show some examples of those that passed
    passed_results = [r for r in results if r.passed_filter]
    if passed_results:
        print(f"\n🏆 Top 5 candidates that passed the filter:")
        sorted_passed = sorted(passed_results, key=lambda x: x.identity_percent, reverse=True)[:5]
        for i, r in enumerate(sorted_passed, 1):
            print(f"   {i}. {r.candidate_id}: {r.identity_percent:.1f}% identity")
    
    print("\n" + "="*70)


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main() -> None:
    """
    Main function that orchestrates the homology filtering process.
    """
    print("\n" + "🧬"*30)
    print("   DEEP-PETASE-MINING - Phase 2: Homology Filtering")
    print("🧬"*30 + "\n")
    
    # Step 1: Load reference
    logger.info("Step 1/5: Loading reference sequence...")
    try:
        reference = load_reference(REFERENCE_FILE)
    except FileNotFoundError:
        sys.exit(1)
    
    # Step 2: Load candidates
    logger.info("Step 2/5: Loading candidate sequences...")
    try:
        candidates = load_candidates(CANDIDATES_FILE)
    except FileNotFoundError:
        sys.exit(1)
    
    # Step 3: Configure aligner
    logger.info("Step 3/5: Configuring local aligner (Smith-Waterman)...")
    aligner = create_aligner()
    
    # Step 4: Perform filtering
    logger.info("Step 4/5: Performing alignments and filtering...")
    print(f"\n⚔️ Aligning {len(candidates)} candidates against IsPETase...")
    print(f"   Criterion: {MIN_IDENTITY}% < identity < {MAX_IDENTITY}%\n")
    
    passed_sequences, all_results = filter_candidates(reference, candidates, aligner)
    
    # Step 5: Save results
    logger.info("Step 5/5: Saving results...")
    save_filtered_fasta(passed_sequences, OUTPUT_FASTA)
    save_results_csv(all_results, OUTPUT_CSV)
    
    # Show summary
    print_summary(len(candidates), len(passed_sequences), all_results)
    
    # Final message with recommendations
    print(f"\n📁 Generated files:")
    print(f"   • Filtered FASTA: {OUTPUT_FASTA}")
    print(f"   • Scores CSV:     {OUTPUT_CSV}")
    
    # Recommendations based on results
    if len(passed_sequences) == 0:
        print(f"\n⚠️ WARNING: No sequence passed the filter!")
        print(f"   Suggestions:")
        print(f"   • Reduce MIN_IDENTITY to 15% or 10%")
        print(f"   • Verify that candidate sequences are proteins")
        print(f"   • Consider expanding the NCBI search")
    elif len(passed_sequences) > 50:
        print(f"\n⚠️ You have {len(passed_sequences)} candidates (too many for Phase 3)")
        print(f"   Suggestion: Raise MIN_IDENTITY to 25% or 30%")
    else:
        print(f"\n✅ {len(passed_sequences)} candidates ready for motif analysis")
    
    print("\n👉 Next step: python 05_motif_scan.py")
    print("\n" + "🧬"*30 + "\n")


# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n⚠️ Execution cancelled by user.")
    except Exception as e:
        logger.error(f"❌ Fatal error: {e}")
        raise
