#!/usr/bin/env python3
"""
13_run_virtual_screening.py - Phase 4: Massive Molecular Docking with AutoDock Vina

This script performs virtual screening (Blind Docking) of all candidate proteins
against the MHET ligand.

Strategy: Blind Docking
- A Grid Box covering the entire protein + 5Å margin is calculated
- This allows finding binding sites without prior knowledge of the active site

Input: 
    - results/structures/pdbqt/*.pdbqt (receptors)
    - data/ligands/mhet.pdbqt (ligand)
    
Output: 
    - results/docking_results/*.pdbqt (poses)
    - results/final_screening.csv (scores)

Features:
    - Checkpointing: Automatically resumes if interrupted
    - Optional parallel processing
    - Pilot mode for testing (--pilot 5)
"""

import os
import sys
import subprocess
import shutil
import argparse
import csv
import re
from pathlib import Path
from datetime import datetime
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR = Path(__file__).parent.resolve()
PDBQT_DIR = BASE_DIR / "results" / "structures" / "pdbqt"
LIGAND_FILE = BASE_DIR / "data" / "ligands" / "mhet.pdbqt"
DOCKING_OUTPUT_DIR = BASE_DIR / "results" / "docking_results"
RESULTS_CSV = BASE_DIR / "results" / "final_screening.csv"

# Vina Parameters
EXHAUSTIVENESS = 8       # 8 = standard, 4 = fast, 16 = more accurate
NUM_MODES = 9            # Number of poses to generate
ENERGY_RANGE = 3         # Energy range for poses (kcal/mol)
GRID_PADDING = 5.0       # Extra margin for Grid Box (Angstroms)

# Parallelism
NUM_WORKERS = 1          # Vina is already parallel, better to run sequentially


# =============================================================================
# FUNCTIONS
# =============================================================================

def check_vina():
    """Verifies that AutoDock Vina is installed."""
    if shutil.which("vina") is None:
        print("❌ Error: AutoDock Vina is not installed.")
        print("   Install with: sudo apt install autodock-vina")
        return False
    
    # Verify version
    result = subprocess.run(["vina", "--version"], capture_output=True, text=True)
    version = result.stdout.strip() if result.stdout else result.stderr.strip()
    print(f"✓ AutoDock Vina detected: {version}")
    return True


def create_directories():
    """Creates necessary directories."""
    DOCKING_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"✓ Output directory: {DOCKING_OUTPUT_DIR}")


def calculate_grid_box(pdbqt_file):
    """
    Calculates the center and size of the Grid Box to cover the entire protein.
    
    Args:
        pdbqt_file: Path to the receptor PDBQT file
        
    Returns:
        dict: {center_x, center_y, center_z, size_x, size_y, size_z}
    """
    coords = []
    
    with open(pdbqt_file, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
                except (ValueError, IndexError):
                    continue
    
    if not coords:
        return None
    
    # Calculate limits
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)
    
    # Protein center
    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2
    center_z = (min_z + max_z) / 2
    
    # Size with padding
    size_x = (max_x - min_x) + 2 * GRID_PADDING
    size_y = (max_y - min_y) + 2 * GRID_PADDING
    size_z = (max_z - min_z) + 2 * GRID_PADDING
    
    # Vina has a limit of 126 Å per dimension
    size_x = min(size_x, 126.0)
    size_y = min(size_y, 126.0)
    size_z = min(size_z, 126.0)
    
    return {
        'center_x': round(center_x, 3),
        'center_y': round(center_y, 3),
        'center_z': round(center_z, 3),
        'size_x': round(size_x, 3),
        'size_y': round(size_y, 3),
        'size_z': round(size_z, 3)
    }


def parse_vina_output(stdout):
    """
    Parses the Vina output to extract the best affinity.
    
    Returns:
        float: Best binding affinity (kcal/mol), or None if parsing failed
    """
    # Search for the results table
    # Typical format:
    #    1     -7.2      0.000      0.000
    #    2     -6.8      1.234      2.345
    
    lines = stdout.split('\n')
    
    for line in lines:
        # The first line with mode "1" has the best affinity
        match = re.match(r'\s*1\s+(-?\d+\.?\d*)', line)
        if match:
            return float(match.group(1))
    
    # Alternative: search for "Affinity:" in newer versions
    for line in lines:
        if "Affinity:" in line:
            match = re.search(r'(-?\d+\.?\d*)', line)
            if match:
                return float(match.group(1))
    
    return None


def run_docking(receptor_file, ligand_file, output_file, grid_box):
    """
    Executes AutoDock Vina for a receptor-ligand pair.
    
    Returns:
        tuple: (success, affinity, message)
    """
    cmd = [
        "vina",
        "--receptor", str(receptor_file),
        "--ligand", str(ligand_file),
        "--out", str(output_file),
        "--center_x", str(grid_box['center_x']),
        "--center_y", str(grid_box['center_y']),
        "--center_z", str(grid_box['center_z']),
        "--size_x", str(grid_box['size_x']),
        "--size_y", str(grid_box['size_y']),
        "--size_z", str(grid_box['size_z']),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", str(NUM_MODES),
        "--energy_range", str(ENERGY_RANGE),
        "--cpu", "1"  # One CPU per process
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600  # 10 minutes maximum
        )
        
        if result.returncode == 0:
            affinity = parse_vina_output(result.stdout)
            if affinity is not None:
                return (True, affinity, "OK")
            else:
                return (False, None, "Could not parse affinity")
        else:
            error = result.stderr[:100] if result.stderr else "Unknown error"
            return (False, None, error)
            
    except subprocess.TimeoutExpired:
        return (False, None, "Timeout (>10 min)")
    except Exception as e:
        return (False, None, str(e)[:100])


def load_checkpoint():
    """Loads previous results from CSV to resume."""
    completed = {}
    
    if RESULTS_CSV.exists():
        with open(RESULTS_CSV, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                protein_id = row.get('Protein_ID', '')
                if protein_id:
                    completed[protein_id] = row
    
    return completed


def save_result(result_row, write_header=False):
    """Saves a result to CSV."""
    fieldnames = ['Protein_ID', 'Binding_Affinity', 'Grid_Center', 'Grid_Size', 
                  'Status', 'Message', 'Timestamp']
    
    mode = 'w' if write_header else 'a'
    
    with open(RESULTS_CSV, mode, newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerow(result_row)


def dock_single_protein(args):
    """Worker for docking a single protein."""
    receptor_file, ligand_file = args
    protein_id = receptor_file.stem
    
    # Calculate Grid Box
    grid_box = calculate_grid_box(receptor_file)
    if grid_box is None:
        return {
            'Protein_ID': protein_id,
            'Binding_Affinity': None,
            'Grid_Center': None,
            'Grid_Size': None,
            'Status': 'Failed',
            'Message': 'Could not calculate Grid Box',
            'Timestamp': datetime.now().isoformat()
        }
    
    # Output file
    output_file = DOCKING_OUTPUT_DIR / f"{protein_id}_docked.pdbqt"
    
    # Execute docking
    success, affinity, message = run_docking(
        receptor_file, ligand_file, output_file, grid_box
    )
    
    grid_center = f"({grid_box['center_x']}, {grid_box['center_y']}, {grid_box['center_z']})"
    grid_size = f"({grid_box['size_x']}, {grid_box['size_y']}, {grid_box['size_z']})"
    
    return {
        'Protein_ID': protein_id,
        'Binding_Affinity': affinity,
        'Grid_Center': grid_center,
        'Grid_Size': grid_size,
        'Status': 'Success' if success else 'Failed',
        'Message': message,
        'Timestamp': datetime.now().isoformat()
    }


def main():
    """Main function."""
    global EXHAUSTIVENESS
    
    # Command line arguments
    parser = argparse.ArgumentParser(description='Virtual Screening with AutoDock Vina')
    parser.add_argument('--pilot', type=int, default=0,
                       help='Pilot mode: process only N proteins (e.g.: --pilot 5)')
    parser.add_argument('--exhaustiveness', type=int, default=EXHAUSTIVENESS,
                       help=f'Search exhaustiveness (default: {EXHAUSTIVENESS})')
    args = parser.parse_args()
    
    EXHAUSTIVENESS = args.exhaustiveness
    
    print("=" * 70)
    print("PHASE 4 - STEP 3: VIRTUAL SCREENING (MOLECULAR DOCKING)")
    print("=" * 70)
    
    # Verify Vina
    if not check_vina():
        return False
    
    # Verify input files
    if not LIGAND_FILE.exists():
        print(f"❌ Error: Ligand not found: {LIGAND_FILE}")
        print("   Run first: python 11_prep_ligand.py")
        return False
    print(f"✓ Ligand: {LIGAND_FILE}")
    
    if not PDBQT_DIR.exists():
        print(f"❌ Error: PDBQT directory not found: {PDBQT_DIR}")
        print("   Run first: python 12_prep_receptors.py")
        return False
    
    # Create directories
    create_directories()
    
    # Get list of receptors
    receptor_files = sorted(PDBQT_DIR.glob("*.pdbqt"))
    total_receptors = len(receptor_files)
    print(f"📂 PDBQT receptors found: {total_receptors}")
    
    if total_receptors == 0:
        print("❌ No PDBQT files found")
        return False
    
    # Pilot mode
    if args.pilot > 0:
        receptor_files = receptor_files[:args.pilot]
        print(f"🧪 PILOT MODE: Processing only {args.pilot} proteins")
    
    # Load checkpoint
    completed = load_checkpoint()
    print(f"✓ Previously processed: {len(completed)}")
    
    # Filter pending
    pending = [f for f in receptor_files if f.stem not in completed]
    pending_count = len(pending)
    
    if pending_count == 0:
        print("\n✅ All proteins have already been processed")
        return True
    
    print(f"📋 Pending to process: {pending_count}")
    print(f"⚙️  Exhaustiveness: {EXHAUSTIVENESS}")
    
    # Estimate time
    estimated_time = pending_count * 2  # ~2 min per protein
    print(f"⏱️  Estimated time: ~{estimated_time} minutes ({estimated_time/60:.1f} hours)")
    
    print("\n" + "-" * 70)
    print("Starting Virtual Screening...")
    print("-" * 70)
    
    # Write header if it's a new file
    write_header = not RESULTS_CSV.exists() or len(completed) == 0
    if write_header:
        save_result({}, write_header=True)
    
    # Statistics
    success_count = 0
    failed_count = 0
    affinities = []
    
    # Process sequentially (Vina is already parallel internally)
    if TQDM_AVAILABLE:
        iterator = tqdm(pending, desc="Docking", unit="prot")
    else:
        iterator = pending
        print(f"Processing {pending_count} proteins...")
    
    for receptor_file in iterator:
        result = dock_single_protein((receptor_file, LIGAND_FILE))
        
        # Save result immediately (checkpointing)
        save_result(result, write_header=False)
        
        if result['Status'] == 'Success':
            success_count += 1
            if result['Binding_Affinity'] is not None:
                affinities.append(result['Binding_Affinity'])
                
            if TQDM_AVAILABLE:
                iterator.set_postfix({
                    'Affinity': f"{result['Binding_Affinity']:.1f}",
                    'OK': success_count
                })
        else:
            failed_count += 1
            if not TQDM_AVAILABLE:
                print(f"  ✗ {result['Protein_ID']}: {result['Message']}")
    
    # Final summary
    print("\n" + "=" * 70)
    print("VIRTUAL SCREENING SUMMARY")
    print("=" * 70)
    
    total_processed = success_count + failed_count + len(completed)
    print(f"\nTotal proteins processed: {total_processed}")
    print(f"✓ Successful docking: {success_count + len([c for c in completed.values() if c.get('Status') == 'Success'])}")
    print(f"✗ Failed: {failed_count}")
    
    if affinities:
        best_affinity = min(affinities)
        avg_affinity = sum(affinities) / len(affinities)
        print(f"\nAffinity statistics:")
        print(f"  - Best affinity: {best_affinity:.2f} kcal/mol")
        print(f"  - Average: {avg_affinity:.2f} kcal/mol")
        print(f"  - Hits (< -7.0): {sum(1 for a in affinities if a < -7.0)}")
    
    print(f"\nResults saved to: {RESULTS_CSV}")
    print(f"Docking poses in: {DOCKING_OUTPUT_DIR}")
    
    print("\n" + "=" * 70)
    print("✅ VIRTUAL SCREENING COMPLETED")
    print("=" * 70)
    print("\nNext step: Run 14_select_top_hits.py")
    
    return True


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
