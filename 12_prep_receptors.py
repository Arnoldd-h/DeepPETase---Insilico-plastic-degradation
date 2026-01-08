#!/usr/bin/env python3
"""
12_prep_receptors.py - Phase 4: Batch Receptor (Protein) Preparation

This script converts all high-confidence PDB files to PDBQT format
for use with AutoDock Vina.

Input: results/structures/*.pdb (excluding low_confidence/)
Output: results/structures/pdbqt/*.pdbqt
"""

import os
import subprocess
import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False
    print("⚠️ tqdm not installed. Install with: pip install tqdm")

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR = Path(__file__).parent.resolve()
PDB_DIR = BASE_DIR / "results" / "structures"
PDBQT_DIR = BASE_DIR / "results" / "structures" / "pdbqt"

# Number of CPUs for parallel processing
NUM_CPUS = max(1, multiprocessing.cpu_count() - 1)  # Leave 1 CPU free


# =============================================================================
# FUNCTIONS
# =============================================================================

def check_obabel():
    """Verify that OpenBabel is installed."""
    if shutil.which("obabel") is None:
        print("❌ Error: OpenBabel (obabel) is not installed.")
        print("   Install with: sudo apt install openbabel")
        return False
    print("✓ OpenBabel detected")
    return True


def create_directories():
    """Create the PDBQT directory if it does not exist."""
    PDBQT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"✓ Output directory: {PDBQT_DIR}")


def get_pdb_files():
    """
    Get the list of PDB files to convert.
    Excludes the low_confidence/ folder.
    """
    pdb_files = []
    
    for pdb_file in PDB_DIR.glob("*.pdb"):
        # Only files in the main directory (not subdirectories)
        if pdb_file.is_file():
            pdb_files.append(pdb_file)
    
    return sorted(pdb_files)


def get_already_converted():
    """Get the list of already converted proteins (for resume)."""
    converted = set()
    
    for pdbqt_file in PDBQT_DIR.glob("*.pdbqt"):
        protein_id = pdbqt_file.stem
        converted.add(protein_id)
    
    return converted


def convert_pdb_to_pdbqt(pdb_file):
    """
    Convert a PDB file to PDBQT using OpenBabel.
    
    Args:
        pdb_file: Path to the PDB file
        
    Returns:
        tuple: (protein_id, success, message)
    """
    protein_id = pdb_file.stem
    output_file = PDBQT_DIR / f"{protein_id}.pdbqt"
    
    # If already exists, skip
    if output_file.exists() and output_file.stat().st_size > 0:
        return (protein_id, True, "Already existed")
    
    cmd = [
        "obabel",
        "-ipdb", str(pdb_file),
        "-opdbqt",
        "-O", str(output_file),
        "-xr",                          # Rigid molecule (receptor)
        "-h",                           # Add hydrogens
        "--partialcharge", "gasteiger"  # Gasteiger charges
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120  # 2 minutes maximum per protein
        )
        
        if result.returncode == 0 and output_file.exists():
            if output_file.stat().st_size > 0:
                return (protein_id, True, "Converted")
            else:
                return (protein_id, False, "Empty file")
        else:
            error_msg = result.stderr[:100] if result.stderr else "Unknown error"
            return (protein_id, False, error_msg)
            
    except subprocess.TimeoutExpired:
        return (protein_id, False, "Timeout")
    except Exception as e:
        return (protein_id, False, str(e)[:100])


def convert_single_worker(pdb_file):
    """Worker for parallel processing."""
    return convert_pdb_to_pdbqt(pdb_file)


def main():
    """Main function."""
    print("=" * 70)
    print("PHASE 4 - STEP 2: BATCH RECEPTOR PREPARATION")
    print("=" * 70)
    
    # Verificar OpenBabel
    if not check_obabel():
        return False
    
    # Crear directorios
    create_directories()
    
    # Get PDB files
    pdb_files = get_pdb_files()
    total_files = len(pdb_files)
    print(f"\n📂 PDB files found: {total_files}")
    
    if total_files == 0:
        print("❌ No PDB files found")
        return False
    
    # Check how many are already converted
    already_converted = get_already_converted()
    print(f"✓ Previously converted: {len(already_converted)}")
    
    # Filter those that still need to be converted
    pending_files = [f for f in pdb_files if f.stem not in already_converted]
    pending_count = len(pending_files)
    
    if pending_count == 0:
        print("\n✅ All files are already converted")
        return True
    
    print(f"📋 Pending conversion: {pending_count}")
    print(f"🖥️  Using {NUM_CPUS} CPUs for parallel conversion")
    
    # Statistics
    success_count = len(already_converted)
    failed_count = 0
    failed_proteins = []
    
    print("\n" + "-" * 70)
    print("Converting PDB to PDBQT...")
    print("-" * 70)
    
    # Parallel processing
    with ProcessPoolExecutor(max_workers=NUM_CPUS) as executor:
        futures = {executor.submit(convert_single_worker, pdb): pdb for pdb in pending_files}
        
        if TQDM_AVAILABLE:
            iterator = tqdm(as_completed(futures), total=pending_count, 
                          desc="Converting", unit="prot")
        else:
            iterator = as_completed(futures)
            print(f"Processing {pending_count} proteins...")
        
        for future in iterator:
            try:
                protein_id, success, message = future.result()
                
                if success:
                    success_count += 1
                else:
                    failed_count += 1
                    failed_proteins.append((protein_id, message))
                    
            except Exception as e:
                failed_count += 1
                failed_proteins.append(("Unknown", str(e)[:50]))
    
    # Final summary
    print("\n" + "=" * 70)
    print("CONVERSION SUMMARY")
    print("=" * 70)
    print(f"\nTotal PDB files: {total_files}")
    print(f"✓ Successfully converted: {success_count}")
    print(f"✗ Failed: {failed_count}")
    
    if failed_proteins:
        print(f"\nProteins with errors (first 10):")
        for protein_id, error in failed_proteins[:10]:
            print(f"  - {protein_id}: {error}")
    
    # Verify generated files
    generated = list(PDBQT_DIR.glob("*.pdbqt"))
    print(f"\nPDBQT files generated: {len(generated)}")
    print(f"Location: {PDBQT_DIR}")
    
    print("\n" + "=" * 70)
    if failed_count == 0:
        print("✅ ALL RECEPTORS PREPARED SUCCESSFULLY")
    else:
        print(f"⚠️ CONVERSION COMPLETED WITH {failed_count} ERRORS")
    print("=" * 70)
    print("\nNext step: Run 13_run_virtual_screening.py")
    
    return failed_count == 0


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
