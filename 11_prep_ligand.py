#!/usr/bin/env python3
"""
11_prep_ligand.py - Phase 4: MHET Ligand Preparation for Docking

This script downloads and prepares the MHET (Mono-(2-hydroxyethyl) terephthalate)
ligand for molecular docking with AutoDock Vina.

MHET is the main intermediate product of PET degradation.
If an enzyme can bind to MHET, it likely has PETase activity.

Input: Download from PubChem (CID 1550473)
Output: data/ligands/mhet.pdbqt
"""

import os
import subprocess
import shutil
from pathlib import Path
import requests
import time

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR = Path(__file__).parent.resolve()
LIGANDS_DIR = BASE_DIR / "data" / "ligands"

# PubChem CID for MHET
MHET_CID = "1550473"
MHET_SDF = LIGANDS_DIR / "mhet.sdf"
MHET_PDBQT = LIGANDS_DIR / "mhet.pdbqt"

# PubChem URL to download 3D SDF
PUBCHEM_URL = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{MHET_CID}/SDF?record_type=3d"
PUBCHEM_URL_2D = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{MHET_CID}/SDF"


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
    """Create the ligands directory if it doesn't exist."""
    LIGANDS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"✓ Directory created: {LIGANDS_DIR}")


def download_mhet_sdf():
    """
    Download the MHET structure from PubChem.
    Tries to get the 3D version first, if not available, uses 2D.
    """
    print(f"\n📥 Downloading MHET (CID: {MHET_CID}) from PubChem...")
    
    # Try to download 3D version first
    try:
        print("   Trying 3D structure...")
        response = requests.get(PUBCHEM_URL, timeout=30)
        
        if response.status_code == 200 and len(response.content) > 100:
            with open(MHET_SDF, 'wb') as f:
                f.write(response.content)
            print(f"✓ 3D structure downloaded: {MHET_SDF}")
            return True
    except Exception as e:
        print(f"   ⚠️ Could not obtain 3D structure: {e}")
    
    # If no 3D available, download 2D (OpenBabel will generate 3D)
    try:
        print("   Downloading 2D structure...")
        response = requests.get(PUBCHEM_URL_2D, timeout=30)
        
        if response.status_code == 200:
            with open(MHET_SDF, 'wb') as f:
                f.write(response.content)
            print(f"✓ 2D structure downloaded: {MHET_SDF}")
            print("   (OpenBabel will generate 3D coordinates)")
            return True
        else:
            print(f"❌ HTTP Error: {response.status_code}")
            return False
            
    except requests.exceptions.RequestException as e:
        print(f"❌ Connection error: {e}")
        return False


def convert_sdf_to_pdbqt():
    """
    Convert the SDF file to PDBQT using OpenBabel.
    
    Important flags:
    --gen3d: Generate 3D coordinates if they don't exist
    -h: Add hydrogens
    -p 7.4: Protonate according to physiological pH
    """
    print(f"\n🔄 Converting SDF to PDBQT...")
    
    cmd = [
        "obabel",
        "-isdf", str(MHET_SDF),
        "-opdbqt",
        "-O", str(MHET_PDBQT),
        "--gen3d",           # Generate 3D coordinates
        "-h",                # Add hydrogens
        "-p", "7.4",         # Physiological pH
        "--partialcharge", "gasteiger"  # Gasteiger charges
    ]
    
    print(f"   Comando: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60
        )
        
        if result.returncode == 0 and MHET_PDBQT.exists():
            # Verify the file is not empty
            if MHET_PDBQT.stat().st_size > 0:
                print(f"✓ Conversion successful: {MHET_PDBQT}")
                return True
            else:
                print("❌ Error: Empty PDBQT file")
                return False
        else:
            print(f"❌ Conversion error:")
            print(f"   stdout: {result.stdout}")
            print(f"   stderr: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        print("❌ Timeout during conversion")
        return False
    except Exception as e:
        print(f"❌ Error: {e}")
        return False


def verify_pdbqt():
    """Verify that the PDBQT file is valid."""
    print(f"\n🔍 Verifying PDBQT file...")
    
    if not MHET_PDBQT.exists():
        print("❌ PDBQT file not found")
        return False
    
    with open(MHET_PDBQT, 'r') as f:
        content = f.read()
    
    # Count atoms
    atom_count = content.count("ATOM") + content.count("HETATM")
    
    # Verify it has the basic structure
    has_atoms = atom_count > 0
    has_charges = "REMARK" in content or atom_count > 0
    
    print(f"   Atoms found: {atom_count}")
    print(f"   File size: {MHET_PDBQT.stat().st_size} bytes")
    
    if has_atoms:
        print("✓ Valid PDBQT file")
        
        # Show first lines
        print("\n   First lines of the file:")
        for line in content.split('\n')[:10]:
            print(f"   {line}")
        
        return True
    else:
        print("❌ Invalid or empty PDBQT file")
        return False


def main():
    """Main function."""
    print("=" * 70)
    print("PHASE 4 - STEP 1: MHET LIGAND PREPARATION")
    print("=" * 70)
    
    # Verify OpenBabel
    if not check_obabel():
        return False
    
    # Create directories
    create_directories()
    
    # Download MHET
    if not download_mhet_sdf():
        return False
    
    # Convert to PDBQT
    if not convert_sdf_to_pdbqt():
        return False
    
    # Verify result
    if not verify_pdbqt():
        return False
    
    print("\n" + "=" * 70)
    print("✅ LIGAND PREPARED SUCCESSFULLY")
    print("=" * 70)
    print(f"\nFile ready for docking: {MHET_PDBQT}")
    print("\nNext step: Run 12_prep_receptors.py")
    
    return True


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
