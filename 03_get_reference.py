#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
03_get_reference.py - Obtaining the Reference Sequence (IsPETase)
================================================================================

Description:
    This script creates the reference file with the sequence of the IsPETase
    enzyme from Ideonella sakaiensis. This enzyme is the "gold standard" for
    PET degradation and will be used as a reference to compare our hypothetical
    candidates.

    IsPETase (PETase from Ideonella sakaiensis):
    - UniProt ID: A0A0K8P6T7
    - PDB: 5XJH, 6EQE
    - Function: Hydrolyzes polyethylene terephthalate (PET)
    - Discovered in: 2016, Osaka, Japan
    - Publication: Yoshida et al., Science 2016

Date: December 2025
Project: Deep-PETase-Mining - Phase 2

Usage:
    python 03_get_reference.py

Output:
    - data/references/ispetase_ref.fasta
================================================================================
"""

from pathlib import Path
from typing import Optional
import logging
import sys

# Biopython for sequence handling
try:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
except ImportError:
    print("❌ Error: Biopython is not installed.")
    print("   Run: pip install biopython")
    sys.exit(1)


# ============================================================================
# CONFIGURATION
# ============================================================================

# Project paths
BASE_DIR: Path = Path(__file__).parent.resolve()
REFERENCES_DIR: Path = BASE_DIR / "data" / "references"
OUTPUT_FILE: Path = REFERENCES_DIR / "ispetase_ref.fasta"

# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# ============================================================================
# REFERENCE SEQUENCE - IsPETase
# ============================================================================

# Amino acid sequence of IsPETase (Ideonella sakaiensis)
# UniProt: A0A0K8P6T7
# This is the mature sequence (without signal peptide) of the PET-degrading enzyme

ISPETASE_SEQUENCE: str = (
    "MNFPRASRLMQAAVLGGLMAVSAAATAQTNPYARGFVYEQNGMKKWGGQGGMQVRTSVRKQLSVPPDQDSWN"
    "GYKGASVQGITPSELISWLPGNPAKEGQSGIFTRAINPSSLKLGQGIGTAGDALIPKEPIYVLSQGGEIYNAS"
    "IGIVGLIHGYPRPLFSWDAIAQQLRKDGAYENIDVEQVCIRSVSWGLNSALGNSLGIAGKGYVAAKSLGLIPV"
    "QIDGLSDHLFGDIDPNWIPDALPLAPGLN"
)

# Sequence information
ISPETASE_ID: str = "A0A0K8P6T7"
ISPETASE_NAME: str = "IsPETase"
ISPETASE_DESCRIPTION: str = (
    "Poly(ethylene terephthalate) hydrolase [Ideonella sakaiensis] "
    "| UniProt: A0A0K8P6T7 | Reference enzyme for PET degradation"
)


# ============================================================================
# FUNCTIONS
# ============================================================================

def create_reference_directory() -> bool:
    """
    Creates the directory for storing reference sequences.
    
    Returns:
        bool: True if created successfully or already exists.
    """
    try:
        REFERENCES_DIR.mkdir(parents=True, exist_ok=True)
        logger.info(f"📁 References directory: {REFERENCES_DIR}")
        return True
    except Exception as e:
        logger.error(f"❌ Error creating directory: {e}")
        return False


def create_reference_record() -> SeqRecord:
    """
    Creates a SeqRecord object with IsPETase information.
    
    Returns:
        SeqRecord: Biopython record with the sequence and metadata.
    
    Note for students:
        SeqRecord is Biopython's main structure for handling sequences.
        It contains the sequence + metadata (ID, description, etc.)
    """
    # Create Seq object (the raw sequence)
    sequence = Seq(ISPETASE_SEQUENCE)
    
    # Create SeqRecord with metadata
    record = SeqRecord(
        seq=sequence,
        id=ISPETASE_ID,
        name=ISPETASE_NAME,
        description=ISPETASE_DESCRIPTION
    )
    
    # Add additional annotations
    record.annotations["molecule_type"] = "protein"
    record.annotations["organism"] = "Ideonella sakaiensis"
    record.annotations["source"] = "UniProt"
    
    return record


def save_reference_fasta(record: SeqRecord, output_path: Path) -> bool:
    """
    Saves the reference sequence in FASTA format.
    
    Args:
        record (SeqRecord): Record with the sequence.
        output_path (Path): Output file path.
    
    Returns:
        bool: True if saved successfully.
    """
    try:
        SeqIO.write(record, output_path, "fasta")
        logger.info(f"💾 Sequence saved to: {output_path}")
        return True
    except Exception as e:
        logger.error(f"❌ Error saving FASTA: {e}")
        return False


def print_sequence_info(record: SeqRecord) -> None:
    """
    Prints detailed information about the reference sequence.
    
    Args:
        record (SeqRecord): Record with the sequence.
    """
    print("\n" + "="*70)
    print("🧬 REFERENCE SEQUENCE - IsPETase")
    print("="*70)
    
    print(f"\n📋 Information:")
    print(f"   • UniProt ID:  {record.id}")
    print(f"   • Name:        {record.name}")
    print(f"   • Organism:    Ideonella sakaiensis")
    print(f"   • Length:      {len(record.seq)} amino acids")
    
    print(f"\n🔬 Enzymatic characteristics:")
    print(f"   • Family:      α/β hydrolase")
    print(f"   • Activity:    PET hydrolase (EC 3.1.1.-)")
    print(f"   • Active site: Catalytic triad Ser-His-Asp")
    print(f"   • Motif:       G-x-S-x-G (serine hydrolase)")
    
    # Display the first and last amino acids
    seq_str = str(record.seq)
    print(f"\n📝 Sequence (preview):")
    print(f"   N-term: {seq_str[:50]}...")
    print(f"   C-term: ...{seq_str[-50:]}")
    
    # Search for the catalytic motif in the sequence
    import re
    motif_pattern = r'G.S.G'
    matches = list(re.finditer(motif_pattern, seq_str))
    
    print(f"\n🎯 G-x-S-x-G motifs found: {len(matches)}")
    for match in matches:
        start, end = match.span()
        print(f"   • Position {start+1}-{end}: {match.group()}")
    
    print("\n" + "="*70)


def main() -> None:
    """
    Main function that orchestrates the creation of the reference file.
    """
    print("\n" + "🧬"*30)
    print("   DEEP-PETASE-MINING - Phase 2: Get Reference")
    print("🧬"*30 + "\n")
    
    # Step 1: Create directory
    logger.info("Step 1/3: Creating references directory...")
    if not create_reference_directory():
        sys.exit(1)
    
    # Step 2: Create sequence record
    logger.info("Step 2/3: Creating sequence record...")
    record = create_reference_record()
    
    # Step 3: Save FASTA file
    logger.info("Step 3/3: Saving FASTA file...")
    if not save_reference_fasta(record, OUTPUT_FILE):
        sys.exit(1)
    
    # Display information
    print_sequence_info(record)
    
    # Final message
    print("\n✅ Reference sequence created successfully!")
    print(f"📄 File: {OUTPUT_FILE}")
    print("\n👉 Next step: python 04_filter_homology.py")
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
