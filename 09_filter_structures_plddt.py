#!/usr/bin/env python3
"""
09_filter_structures_plddt.py - Phase 3: Filtering structures by quality (pLDDT)

This script evaluates the quality of structures predicted by ESMFold
using the pLDDT score (predicted Local Distance Difference Test).

The pLDDT is stored in the B-factor field of PDB files.
Reference values:
  - > 90: Very high confidence
  - 70-90: Reasonable confidence
  - 50-70: Low confidence
  - < 50: Very low confidence

Input: results/structures/*.pdb
Output: 
  - High quality structures: results/structures/ (not moved)
  - Low quality structures: results/structures/low_confidence/
  - Quality CSV: results/logs/structure_quality.csv
"""

import os
import shutil
from pathlib import Path
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
import csv

# Suppress Biopython warnings
warnings.filterwarnings("ignore", category=PDBConstructionWarning)

# =============================================================================
# CONFIGURATION
# =============================================================================
STRUCTURES_DIR = Path("results/structures")
LOW_CONFIDENCE_DIR = STRUCTURES_DIR / "low_confidence"
OUTPUT_CSV = Path("results/logs/structure_quality.csv")

# pLDDT threshold to consider a structure as high quality
PLDDT_THRESHOLD = 70.0

# =============================================================================
# AUXILIARY FUNCTIONS
# =============================================================================

def ensure_directories():
    """Create necessary directories if they don't exist."""
    LOW_CONFIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)


def normalize_plddt(plddt_values):
    """
    Detect the pLDDT scale and normalize to 0-100 scale.
    
    ESMFold can return pLDDT in two scales:
    - Scale 0-1 (normalized)
    - Scale 0-100 (standard)
    
    This function automatically detects the scale and converts to 0-100.
    """
    if not plddt_values:
        return plddt_values
    
    max_val = max(plddt_values)
    
    # If the maximum value is <= 1.0, we assume 0-1 scale
    if max_val <= 1.0:
        return [v * 100 for v in plddt_values]
    else:
        # Already in 0-100 scale
        return plddt_values


def calculate_average_plddt(pdb_file):
    """
    Calculate the average pLDDT of a PDB structure.
    
    In ESMFold/AlphaFold PDB files, pLDDT is stored
    in the B-factor field of each atom.
    
    Parameters:
        pdb_file: Path to the PDB file
        
    Returns:
        tuple: (average_plddt, num_residues, plddt_per_residue)
    """
    parser = PDBParser(QUIET=True)
    
    try:
        structure = parser.get_structure("protein", pdb_file)
    except Exception as e:
        return None, 0, [], str(e)
    
    # Collect pLDDT per residue (using CA atom of each residue)
    plddt_values = []
    residue_plddt = {}
    
    for model in structure:
        for chain in model:
            for residue in chain:
                # Try to get the CA atom (alpha carbon)
                if "CA" in residue:
                    ca_atom = residue["CA"]
                    plddt = ca_atom.get_bfactor()
                    plddt_values.append(plddt)
                    res_id = f"{chain.id}:{residue.get_resname()}{residue.get_id()[1]}"
                    residue_plddt[res_id] = plddt
                else:
                    # If no CA, use the first atom of the residue
                    atoms = list(residue.get_atoms())
                    if atoms:
                        plddt = atoms[0].get_bfactor()
                        plddt_values.append(plddt)
    
    if not plddt_values:
        return None, 0, [], "No valid atoms found"
    
    # Normalize pLDDT to 0-100 scale if it's in 0-1 scale
    plddt_values = normalize_plddt(plddt_values)
    
    average_plddt = sum(plddt_values) / len(plddt_values)
    return average_plddt, len(plddt_values), plddt_values, None


def get_plddt_category(plddt):
    """Categorize pLDDT according to standard ranges."""
    if plddt >= 90:
        return "Very High"
    elif plddt >= 70:
        return "High"
    elif plddt >= 50:
        return "Low"
    else:
        return "Very Low"


def main():
    """Main function of the script."""
    print("=" * 70)
    print("PHASE 3 - STEP 3: FILTERING STRUCTURES BY QUALITY (pLDDT)")
    print("=" * 70)
    
    # Create directories
    ensure_directories()
    
    # Restore files from low_confidence if they exist (for re-running the script)
    low_conf_files = list(LOW_CONFIDENCE_DIR.glob("*.pdb"))
    if low_conf_files:
        print(f"\nRestoring {len(low_conf_files)} files from low_confidence/...")
        for pdb_file in low_conf_files:
            dest_file = STRUCTURES_DIR / pdb_file.name
            shutil.move(str(pdb_file), str(dest_file))
        print("Files restored successfully.")
    
    # Search for PDB files
    pdb_files = list(STRUCTURES_DIR.glob("*.pdb"))
    
    if not pdb_files:
        raise FileNotFoundError(
            f"No PDB files found in: {STRUCTURES_DIR}\n"
            "Did you run 08_predict_structures_esm.py first?"
        )
    
    print(f"\nPDB files found: {len(pdb_files)}")
    print(f"pLDDT threshold: {PLDDT_THRESHOLD}")
    print("-" * 70)
    
    # Statistics
    results = []
    passed_count = 0
    failed_count = 0
    error_count = 0
    
    # Quality categories
    quality_distribution = {
        "Very High": 0,  # >= 90
        "High": 0,       # 70-90
        "Low": 0,        # 50-70
        "Very Low": 0    # < 50
    }
    
    print("\nAnalyzing structures...")
    
    for pdb_file in pdb_files:
        seq_id = pdb_file.stem
        
        # Calculate average pLDDT
        avg_plddt, num_residues, plddt_values, error = calculate_average_plddt(pdb_file)
        
        if error:
            print(f"  ✗ {seq_id}: Error - {error}")
            results.append({
                "ID": seq_id,
                "Average_pLDDT": "N/A",
                "Num_Residues": 0,
                "Min_pLDDT": "N/A",
                "Max_pLDDT": "N/A",
                "Category": "Error",
                "Status": "Error"
            })
            error_count += 1
            continue
        
        # Additional statistics
        min_plddt = min(plddt_values)
        max_plddt = max(plddt_values)
        category = get_plddt_category(avg_plddt)
        quality_distribution[category] += 1
        
        # Determine if it passes the filter
        if avg_plddt >= PLDDT_THRESHOLD:
            status = "Passed"
            passed_count += 1
            print(f"  ✓ {seq_id}: pLDDT = {avg_plddt:.1f} ({category})")
        else:
            status = "Failed"
            failed_count += 1
            # Move to low confidence folder
            dest_file = LOW_CONFIDENCE_DIR / pdb_file.name
            shutil.move(str(pdb_file), str(dest_file))
            print(f"  ✗ {seq_id}: pLDDT = {avg_plddt:.1f} ({category}) -> Moved to low_confidence/")
        
        results.append({
            "ID": seq_id,
            "Average_pLDDT": round(avg_plddt, 2),
            "Num_Residues": num_residues,
            "Min_pLDDT": round(min_plddt, 2),
            "Max_pLDDT": round(max_plddt, 2),
            "Category": category,
            "Status": status
        })
    
    # Save CSV
    print(f"\nSaving results to: {OUTPUT_CSV}")
    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        fieldnames = ["ID", "Average_pLDDT", "Num_Residues", "Min_pLDDT", "Max_pLDDT", "Category", "Status"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    
    # Final summary
    print("\n" + "=" * 70)
    print("QUALITY FILTERING SUMMARY")
    print("=" * 70)
    print(f"\nTotal structures analyzed: {len(pdb_files)}")
    print(f"  ✓ Passed the filter (pLDDT >= {PLDDT_THRESHOLD}): {passed_count}")
    print(f"  ✗ Did not pass the filter: {failed_count}")
    if error_count > 0:
        print(f"  ⚠ Errors: {error_count}")
    
    print(f"\nDistribution by quality category:")
    print(f"  - Very High (>= 90): {quality_distribution['Very High']}")
    print(f"  - High (70-90):      {quality_distribution['High']}")
    print(f"  - Low (50-70):       {quality_distribution['Low']}")
    print(f"  - Very Low (< 50):   {quality_distribution['Very Low']}")
    
    print(f"\nOutput files:")
    print(f"  - High quality structures: {STRUCTURES_DIR}/ ({passed_count} files)")
    print(f"  - Low quality structures: {LOW_CONFIDENCE_DIR}/ ({failed_count} files)")
    print(f"  - Quality report: {OUTPUT_CSV}")
    
    # Calculate overall average pLDDT
    valid_plddt = [r["Average_pLDDT"] for r in results if r["Average_pLDDT"] != "N/A"]
    if valid_plddt:
        overall_avg = sum(valid_plddt) / len(valid_plddt)
        print(f"\nOverall average pLDDT: {overall_avg:.2f}")
    
    print("\n" + "=" * 70)
    if passed_count > 0:
        print(f"Filtering completed. {passed_count} structures ready for molecular docking.")
        print("Next step: Run 10_plot_plddt.py to visualize the distribution")
    else:
        print("⚠️ No structure passed the quality filter")
        print("Consider reducing the pLDDT threshold or reviewing the input sequences")
    print("=" * 70)


if __name__ == "__main__":
    main()
