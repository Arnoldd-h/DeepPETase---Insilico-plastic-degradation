#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
05_motif_scan.py - Catalytic Motif Search in Sequences
================================================================================

Description:
    This script searches for essential catalytic motifs in sequences that
    passed the homology filter. It focuses on identifying the characteristic
    motif of serine hydrolases: G-x-S-x-G.

    Scientific rationale:
    PETases belong to the α/β hydrolase superfamily and contain:
    
    1. Serine hydrolase motif (G-x-S-x-G):
       - Glycine (G): Provides backbone flexibility
       - x: Any amino acid
       - Serine (S): Catalytic nucleophile (ESSENTIAL)
       - x: Any amino acid  
       - Glycine (G): Flexibility for the oxyanion turn

    2. Catalytic triad (Ser-His-Asp):
       - Serine: Nucleophile that attacks the ester bond
       - Histidine: General base that activates serine
       - Aspartate: Orients histidine

    Only sequences with the G-x-S-x-G motif have catalytic potential.

Date: December 2025
Project: Deep-PETase-Mining - Phase 2

Usage:
    python 05_motif_scan.py

Dependencies:
    - Biopython
    - re (regular expressions, standard library)

Input:
    - data/processed/candidates_homology.fasta

Output:
    - data/processed/candidates_final_seqs.fasta
    - results/logs/motif_analysis.csv
================================================================================
"""

from pathlib import Path
from typing import List, Dict, Tuple, Optional, Pattern
from dataclasses import dataclass, field
import re
import logging
import sys

# Biopython
try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
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


# ============================================================================
# CONFIGURATION
# ============================================================================

# Project paths
BASE_DIR: Path = Path(__file__).parent.resolve()
INPUT_FILE: Path = BASE_DIR / "data" / "processed" / "candidates_homology.fasta"
OUTPUT_FASTA: Path = BASE_DIR / "data" / "processed" / "candidates_final_seqs.fasta"
OUTPUT_CSV: Path = BASE_DIR / "results" / "logs" / "motif_analysis.csv"

# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# ============================================================================
# MOTIF DEFINITIONS
# ============================================================================

# Catalytic motifs to search for
# Note for students: Regex patterns for proteins use:
#   . = any amino acid
#   [XYZ] = X, Y or Z
#   {n,m} = between n and m repetitions

MOTIFS: Dict[str, Dict] = {
    "GXSXG": {
        "pattern": r"G.S.G",  # Glycine-X-Serine-X-Glycine
        "description": "Serine hydrolase motif (nucleophilic elbow)",
        "essential": True,  # REQUIRED to pass the filter
        "function": "Contains the catalytic serine in a nucleophile turn"
    },
    "CATALYTIC_SER": {
        "pattern": r"[VILM].S[AGTC]G",  # Broader variant of the serine site
        "description": "Extended serine active site",
        "essential": False,
        "function": "Extended context of the serine active site"
    },
    "HIS_BOX": {
        "pattern": r"[LIVMFY].H[GSAC]",  # Histidine box
        "description": "Histidine box (partial catalytic triad)",
        "essential": False,
        "function": "Histidine of the catalytic triad"
    },
    "ASP_BOX": {
        "pattern": r"[LIVMF].[DE][GSAC]",  # Aspartate/glutamate box
        "description": "Acidic residue box (partial catalytic triad)",
        "essential": False,
        "function": "Acidic residue of the catalytic triad"
    },
    "OXYANION_HOLE": {
        "pattern": r"[HY]G.G",  # Oxyanion hole
        "description": "Oxyanion hole motif",
        "essential": False,
        "function": "Stabilizes the tetrahedral intermediate"
    }
}


# ============================================================================
# DATA STRUCTURES
# ============================================================================

@dataclass
class MotifMatch:
    """Represents a found motif match."""
    motif_name: str
    pattern: str
    sequence_matched: str
    start_position: int
    end_position: int


@dataclass
class SequenceAnalysis:
    """Result of motif analysis for a sequence."""
    seq_id: str
    seq_description: str
    seq_length: int
    motifs_found: List[MotifMatch] = field(default_factory=list)
    has_essential_motif: bool = False
    total_motifs: int = 0
    passed_filter: bool = False


# ============================================================================
# ANALYSIS FUNCTIONS
# ============================================================================

def compile_motif_patterns() -> Dict[str, Pattern]:
    """
    Compiles the regular expression patterns for the motifs.
    
    Returns:
        Dict[str, Pattern]: Dictionary of compiled patterns.
    
    Note for students:
        Compiling regex patterns once is more efficient than
        compiling them on each search when you have many sequences.
    """
    compiled = {}
    for name, info in MOTIFS.items():
        try:
            compiled[name] = re.compile(info["pattern"])
            logger.debug(f"Pattern compiled: {name} -> {info['pattern']}")
        except re.error as e:
            logger.error(f"❌ Error compiling pattern {name}: {e}")
    
    return compiled


def find_motifs_in_sequence(
    sequence: str,
    patterns: Dict[str, Pattern]
) -> Tuple[List[MotifMatch], bool]:
    """
    Searches for all defined motifs in a sequence.
    
    Args:
        sequence (str): Amino acid sequence.
        patterns (Dict[str, Pattern]): Compiled patterns.
    
    Returns:
        Tuple: (list of matches, has essential motif)
    """
    matches: List[MotifMatch] = []
    has_essential = False
    
    for motif_name, pattern in patterns.items():
        # Search for all pattern matches
        for match in pattern.finditer(sequence):
            motif_match = MotifMatch(
                motif_name=motif_name,
                pattern=MOTIFS[motif_name]["pattern"],
                sequence_matched=match.group(),
                start_position=match.start() + 1,  # 1-indexed for biology
                end_position=match.end()
            )
            matches.append(motif_match)
            
            # Check if it's an essential motif
            if MOTIFS[motif_name]["essential"]:
                has_essential = True
    
    return matches, has_essential


def analyze_sequence(
    record: SeqRecord,
    patterns: Dict[str, Pattern]
) -> SequenceAnalysis:
    """
    Analyzes a sequence searching for catalytic motifs.
    
    Args:
        record (SeqRecord): Sequence record.
        patterns (Dict[str, Pattern]): Compiled patterns.
    
    Returns:
        SequenceAnalysis: Analysis result.
    """
    sequence = str(record.seq).upper()
    
    # Search for motifs
    matches, has_essential = find_motifs_in_sequence(sequence, patterns)
    
    # Create result
    analysis = SequenceAnalysis(
        seq_id=record.id,
        seq_description=record.description,
        seq_length=len(sequence),
        motifs_found=matches,
        has_essential_motif=has_essential,
        total_motifs=len(matches),
        passed_filter=has_essential  # Only those with GXSXG pass
    )
    
    return analysis


def analyze_all_sequences(
    records: List[SeqRecord],
    patterns: Dict[str, Pattern]
) -> List[SequenceAnalysis]:
    """
    Analyzes all sequences searching for motifs.
    
    Args:
        records (List[SeqRecord]): List of sequences.
        patterns (Dict[str, Pattern]): Compiled patterns.
    
    Returns:
        List[SequenceAnalysis]: Results of all analyses.
    """
    results = []
    
    for record in records:
        analysis = analyze_sequence(record, patterns)
        results.append(analysis)
    
    return results


# ============================================================================
# SAVE FUNCTIONS
# ============================================================================

def save_filtered_sequences(
    records: List[SeqRecord],
    analyses: List[SequenceAnalysis],
    output_path: Path
) -> List[SeqRecord]:
    """
    Saves the sequences that passed the motif filter.
    
    Args:
        records (List[SeqRecord]): Original sequences.
        analyses (List[SequenceAnalysis]): Analysis results.
        output_path (Path): Output path.
    
    Returns:
        List[SeqRecord]: Sequences that passed the filter.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Create analysis dictionary by ID
    analysis_dict = {a.seq_id: a for a in analyses}
    
    # Filter sequences that passed
    passed_records = [
        record for record in records
        if analysis_dict.get(record.id, SequenceAnalysis(
            seq_id="", seq_description="", seq_length=0
        )).passed_filter
    ]
    
    # Save
    if passed_records:
        SeqIO.write(passed_records, output_path, "fasta")
        logger.info(f"💾 {len(passed_records)} sequences saved to: {output_path}")
    else:
        # Create empty file with message
        output_path.write_text("")
        logger.warning(f"⚠️ Empty file created (no sequence passed the filter)")
    
    return passed_records


def save_analysis_csv(
    analyses: List[SequenceAnalysis],
    output_path: Path
) -> None:
    """
    Saves the motif analysis in CSV format.
    
    Args:
        analyses (List[SequenceAnalysis]): Analysis results.
        output_path (Path): Path to the CSV file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Prepare data for CSV
    data = []
    for a in analyses:
        # Summarize found motifs
        motif_summary = "; ".join([
            f"{m.motif_name}:{m.sequence_matched}@{m.start_position}"
            for m in a.motifs_found
        ]) if a.motifs_found else "None"
        
        # Count each motif type
        motif_counts = {}
        for m in a.motifs_found:
            motif_counts[m.motif_name] = motif_counts.get(m.motif_name, 0) + 1
        
        data.append({
            "ID": a.seq_id,
            "Length": a.seq_length,
            "Has_GXSXG": a.has_essential_motif,
            "Total_Motifs": a.total_motifs,
            "GXSXG_Count": motif_counts.get("GXSXG", 0),
            "HIS_BOX_Count": motif_counts.get("HIS_BOX", 0),
            "ASP_BOX_Count": motif_counts.get("ASP_BOX", 0),
            "Passed_Filter": a.passed_filter,
            "Motifs_Detail": motif_summary[:200]  # Truncate if too long
        })
    
    df = pd.DataFrame(data)
    df = df.sort_values("Total_Motifs", ascending=False)
    df.to_csv(output_path, index=False)
    
    logger.info(f"💾 Analysis saved to: {output_path}")


# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

def print_motif_legend() -> None:
    """Prints a legend of the motifs being searched."""
    print("\n" + "-"*70)
    print("📚 CATALYTIC MOTIFS BEING SEARCHED:")
    print("-"*70)
    
    for name, info in MOTIFS.items():
        essential = "⭐ ESSENTIAL" if info["essential"] else "   Optional"
        print(f"\n{essential} | {name}")
        print(f"   Pattern:     {info['pattern']}")
        print(f"   Description: {info['description']}")
        print(f"   Function:    {info['function']}")


def print_summary(
    total: int,
    passed: int,
    analyses: List[SequenceAnalysis]
) -> None:
    """
    Prints a summary of the motif analysis.
    
    Args:
        total (int): Total sequences analyzed.
        passed (int): Sequences that passed the filter.
        analyses (List[SequenceAnalysis]): All analyses.
    """
    print("\n" + "="*70)
    print("📊 MOTIF ANALYSIS SUMMARY")
    print("="*70)
    
    print(f"\n📈 General statistics:")
    print(f"   • Sequences analyzed:         {total:,}")
    print(f"   • With G-x-S-x-G motif:       {passed:,}")
    print(f"   • Without essential motif:    {total - passed:,}")
    print(f"   • Selection rate:             {(passed/total)*100:.1f}%" if total > 0 else "   • Selection rate: N/A")
    
    # Motif statistics
    total_gxsxg = sum(1 for a in analyses if a.has_essential_motif)
    motif_counts = {}
    for a in analyses:
        for m in a.motifs_found:
            motif_counts[m.motif_name] = motif_counts.get(m.motif_name, 0) + 1
    
    print(f"\n🔍 Distribution of found motifs:")
    for name, count in sorted(motif_counts.items(), key=lambda x: x[1], reverse=True):
        essential = "⭐" if MOTIFS[name]["essential"] else "  "
        print(f"   {essential} {name}: {count} matches")
    
    # Show finalists
    passed_analyses = [a for a in analyses if a.passed_filter]
    if passed_analyses:
        print(f"\n🏆 FINALIST SEQUENCES ({len(passed_analyses)}):")
        print("-"*70)
        
        # Sort by number of motifs
        sorted_passed = sorted(passed_analyses, key=lambda x: x.total_motifs, reverse=True)
        
        for i, a in enumerate(sorted_passed, 1):
            gxsxg_matches = [m for m in a.motifs_found if m.motif_name == "GXSXG"]
            gxsxg_str = ", ".join([f"{m.sequence_matched}@{m.start_position}" for m in gxsxg_matches])
            print(f"\n   {i}. {a.seq_id}")
            print(f"      Length: {a.seq_length} aa | Motifs: {a.total_motifs}")
            print(f"      G-x-S-x-G: {gxsxg_str}")
    
    print("\n" + "="*70)


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main() -> None:
    """
    Main function that orchestrates the motif analysis.
    """
    print("\n" + "🧬"*30)
    print("   DEEP-PETASE-MINING - Phase 2: Motif Search")
    print("🧬"*30 + "\n")
    
    # Verify input file
    if not INPUT_FILE.exists():
        logger.error(f"❌ File not found: {INPUT_FILE}")
        logger.error("   Run first: python 04_filter_homology.py")
        sys.exit(1)
    
    # Step 1: Load sequences
    logger.info("Step 1/4: Loading homology-filtered sequences...")
    records = list(SeqIO.parse(INPUT_FILE, "fasta"))
    logger.info(f"✅ Loaded {len(records)} sequences")
    
    if len(records) == 0:
        logger.error("❌ The input file is empty")
        logger.error("   Check the homology filtering or relax the criteria")
        sys.exit(1)
    
    # Show motif legend
    print_motif_legend()
    
    # Step 2: Compile patterns
    logger.info("\nStep 2/4: Compiling search patterns...")
    patterns = compile_motif_patterns()
    
    # Step 3: Analyze sequences
    logger.info("Step 3/4: Analyzing sequences for motifs...")
    print(f"\n🔬 Searching for catalytic motifs in {len(records)} sequences...\n")
    
    analyses = analyze_all_sequences(records, patterns)
    
    # Step 4: Save results
    logger.info("Step 4/4: Saving results...")
    passed_records = save_filtered_sequences(records, analyses, OUTPUT_FASTA)
    save_analysis_csv(analyses, OUTPUT_CSV)
    
    # Show summary
    print_summary(len(records), len(passed_records), analyses)
    
    # Final message
    print(f"\n📁 Generated files:")
    print(f"   • Finalist FASTA: {OUTPUT_FASTA}")
    print(f"   • Analysis CSV:   {OUTPUT_CSV}")
    
    # Recommendations
    if len(passed_records) == 0:
        print(f"\n⚠️ WARNING: No sequence has the G-x-S-x-G motif")
        print(f"   This may mean:")
        print(f"   1. The sequences are not serine hydrolases")
        print(f"   2. The motif is mutated (consider searching for G.S or just S)")
        print(f"   3. You need more candidate sequences")
    elif len(passed_records) > 50:
        print(f"\n⚠️ You have {len(passed_records)} finalists (too many for Phase 3)")
        print(f"   Suggestion: Increase the identity filter in 04_filter_homology.py")
    elif 5 <= len(passed_records) <= 20:
        print(f"\n✅ Excellent! {len(passed_records)} candidates is an ideal number")
        print(f"   Ready for structural prediction in Phase 3")
    else:
        print(f"\n✅ {len(passed_records)} candidates ready for the next phase")
    
    print("\n👉 Next step: python 06_plot_filtering.py")
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
