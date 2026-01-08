#!/usr/bin/env python3
"""
08_predict_structures_esm.py - Phase 3: 3D Structure Prediction with ESMFold

This script uses the public ESMFold API to predict the 3D structure
of representative sequences obtained from clustering.

Input: data/processed/candidates_clustered.fasta
Output: results/structures/*.pdb
"""

import os
import time
import requests
from pathlib import Path
from Bio import SeqIO
from tqdm import tqdm
import json

# =============================================================================
# CONFIGURATION
# =============================================================================
INPUT_FASTA = Path("data/processed/candidates_clustered.fasta")
OUTPUT_DIR = Path("results/structures")
LOG_FILE = Path("results/logs/esm_predictions.csv")

# ESMFold API
ESMATLAS_API = "https://api.esmatlas.com/foldSequence/v1/pdb/"

# API Parameters
MAX_SEQUENCE_LENGTH = 400  # ESMFold has a length limit
REQUEST_TIMEOUT = 600  # 10 minutes timeout (increased for long sequences)
CONNECT_TIMEOUT = 30  # Initial connection timeout
DELAY_BETWEEN_REQUESTS = 3  # Seconds between requests (be nice to the API)
MAX_RETRIES = 5  # Maximum number of retries per sequence

# =============================================================================
# AUXILIARY FUNCTIONS
# =============================================================================

def ensure_directories():
    """Create necessary directories if they don't exist."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)


def clean_sequence(sequence):
    """
    Clean the sequence for the ESMFold API.
    - Removes non-standard characters
    - Replaces ambiguous amino acids
    """
    # Standard amino acids
    standard_aa = set("ACDEFGHIKLMNPQRSTVWY")
    
    # Clean sequence
    cleaned = ""
    for aa in str(sequence).upper():
        if aa in standard_aa:
            cleaned += aa
        elif aa == 'X':
            cleaned += 'A'  # Replace X with Alanine
        elif aa == 'U':
            cleaned += 'C'  # Selenocysteine -> Cysteine
        elif aa == 'O':
            cleaned += 'K'  # Pyrrolysine -> Lysine
        # Ignore other characters (gaps, etc.)
    
    return cleaned


def predict_structure_esm(sequence, seq_id, retries=MAX_RETRIES):
    """
    Predict the 3D structure of a sequence using the ESMFold API.
    
    Parameters:
        sequence: Amino acid sequence
        seq_id: Sequence identifier
        retries: Number of retries in case of error
        
    Returns:
        tuple: (success: bool, pdb_content or error_message: str)
    """
    # Clean sequence
    clean_seq = clean_sequence(sequence)
    
    # Check length
    if len(clean_seq) > MAX_SEQUENCE_LENGTH:
        return False, f"Sequence too long ({len(clean_seq)} > {MAX_SEQUENCE_LENGTH} aa)"
    
    if len(clean_seq) < 10:
        return False, f"Sequence too short ({len(clean_seq)} aa)"
    
    for attempt in range(retries):
        try:
            # Make POST request to API with separate timeout for connection and reading
            response = requests.post(
                ESMATLAS_API,
                data=clean_seq,
                headers={"Content-Type": "text/plain"},
                timeout=(CONNECT_TIMEOUT, REQUEST_TIMEOUT)  # (connection, reading)
            )
            
            # Check response
            if response.status_code == 200:
                return True, response.text
            
            elif response.status_code == 429:
                # Rate limiting - wait longer
                wait_time = 60 * (attempt + 1)  # Increased: 60, 120, 180s...
                print(f"    Rate limit reached. Waiting {wait_time}s...")
                time.sleep(wait_time)
                continue
                
            elif response.status_code == 503:
                # Service unavailable - retry
                wait_time = 30 * (attempt + 1)  # Increased
                print(f"    Service unavailable. Retrying in {wait_time}s...")
                time.sleep(wait_time)
                continue
            
            elif response.status_code == 504:
                # Gateway timeout - prediction is taking too long
                wait_time = 60 * (attempt + 1)
                print(f"    Gateway timeout. Retrying in {wait_time}s...")
                time.sleep(wait_time)
                continue
                
            else:
                return False, f"HTTP Error {response.status_code}: {response.text[:200]}"
                
        except requests.exceptions.Timeout:
            if attempt < retries - 1:
                wait_time = 30 * (attempt + 1)
                print(f"    Timeout. Waiting {wait_time}s and retrying ({attempt + 2}/{retries})...")
                time.sleep(wait_time)
                continue
            return False, "Timeout after multiple attempts"
        
        except requests.exceptions.ConnectionError:
            if attempt < retries - 1:
                wait_time = 15 * (attempt + 1)
                print(f"    Connection error. Waiting {wait_time}s and retrying ({attempt + 2}/{retries})...")
                time.sleep(wait_time)
                continue
            return False, "Connection error after multiple attempts"
            
        except requests.exceptions.RequestException as e:
            if attempt < retries - 1:
                wait_time = 10 * (attempt + 1)
                print(f"    Network error. Waiting {wait_time}s and retrying ({attempt + 2}/{retries})...")
                time.sleep(wait_time)
                continue
            return False, f"Connection error: {str(e)}"
    
    return False, "Maximum number of retries reached"


def main():
    """Main function of the script."""
    print("=" * 70)
    print("PHASE 3 - STEP 2: 3D STRUCTURE PREDICTION (ESMFold)")
    print("=" * 70)
    
    # Create directories
    ensure_directories()
    
    # Check input file
    if not INPUT_FASTA.exists():
        raise FileNotFoundError(
            f"Input file not found: {INPUT_FASTA}\n"
            "Did you run 07_cluster_sequences.py first?"
        )
    
    # Load sequences
    print(f"\nLoading sequences from: {INPUT_FASTA}")
    sequences = list(SeqIO.parse(INPUT_FASTA, "fasta"))
    total_sequences = len(sequences)
    print(f"Sequences to process: {total_sequences}")
    
    # Check existing sequences (to resume interrupted predictions)
    existing_pdbs = set(f.stem for f in OUTPUT_DIR.glob("*.pdb"))
    if existing_pdbs:
        print(f"Already existing structures: {len(existing_pdbs)}")
        sequences = [s for s in sequences if s.id not in existing_pdbs]
        print(f"Pending sequences: {len(sequences)}")
    
    if not sequences:
        print("\n✓ All sequences have already been processed")
        return
    
    # Statistics
    stats = {
        "successful": 0,
        "failed": 0,
        "too_long": 0,
        "errors": []
    }
    
    # Log file
    log_exists = LOG_FILE.exists()
    log_file = open(LOG_FILE, "a", encoding="utf-8")
    if not log_exists:
        log_file.write("Sequence_ID,Length,Status,Message\n")
    
    print(f"\nStarting predictions...")
    print(f"API: {ESMATLAS_API}")
    print(f"Delay between requests: {DELAY_BETWEEN_REQUESTS}s")
    print("-" * 70)
    
    # Process sequences with progress bar
    for record in tqdm(sequences, desc="Predicting structures", unit="seq"):
        seq_id = record.id
        sequence = str(record.seq)
        seq_length = len(sequence)
        
        # Clean ID for filename
        safe_id = "".join(c if c.isalnum() or c in "._-" else "_" for c in seq_id)
        output_file = OUTPUT_DIR / f"{safe_id}.pdb"
        
        # Predict structure
        success, result = predict_structure_esm(sequence, seq_id)
        
        if success:
            # Save PDB file
            with open(output_file, "w") as f:
                f.write(result)
            stats["successful"] += 1
            log_file.write(f"{seq_id},{seq_length},Success,Saved to {output_file.name}\n")
        else:
            stats["failed"] += 1
            if "too long" in result:
                stats["too_long"] += 1
            stats["errors"].append((seq_id, result))
            log_file.write(f'{seq_id},{seq_length},Failed,"{result}"\n')
            tqdm.write(f"  ✗ {seq_id}: {result}")
        
        # Flush log
        log_file.flush()
        
        # Wait between requests
        time.sleep(DELAY_BETWEEN_REQUESTS)
    
    log_file.close()
    
    # Final summary
    print("\n" + "=" * 70)
    print("STRUCTURE PREDICTION SUMMARY")
    print("=" * 70)
    print(f"Total sequences processed: {len(sequences)}")
    print(f"  ✓ Successful: {stats['successful']}")
    print(f"  ✗ Failed: {stats['failed']}")
    if stats["too_long"] > 0:
        print(f"    - Too long (>{MAX_SEQUENCE_LENGTH} aa): {stats['too_long']}")
    
    # Count total PDB files
    total_pdbs = len(list(OUTPUT_DIR.glob("*.pdb")))
    print(f"\nPDB files generated: {total_pdbs}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Prediction log: {LOG_FILE}")
    
    if stats["errors"]:
        print(f"\nFirst errors found:")
        for seq_id, error in stats["errors"][:5]:
            print(f"  - {seq_id}: {error[:50]}...")
    
    print("\n" + "=" * 70)
    if stats["successful"] > 0:
        print("Structure prediction completed")
        print("Next step: Run 09_filter_structures_plddt.py to filter by quality")
    else:
        print("⚠️ No structures were generated successfully")
        print("Check your internet connection and API availability")
    print("=" * 70)


if __name__ == "__main__":
    main()
