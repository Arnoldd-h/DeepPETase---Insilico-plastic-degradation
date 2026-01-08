#!/usr/bin/env python3
"""
14_select_top_hits.py - Phase 4: Top Candidate Selection

This script analyzes the virtual screening results and selects
the best candidates based on their binding affinity.

Selection criteria:
    - Affinity threshold: < -7.0 kcal/mol (more negative = better binding)
    - Alternative: Top 1% by affinity

Input:
    - results/final_screening.csv (docking results)
    - results/structures/*.pdb (original structures)
    
Output:
    - results/top_hits_pdb/ (Top 10 PDBs copied)
    - results/top_hits_summary.csv (Candidate summary)
    - results/figures/affinity_distribution.png
    - results/figures/top_hits_barplot.png
"""

import os
import sys
import shutil
import argparse
import csv
from pathlib import Path
from datetime import datetime

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False
    print("⚠️  Pandas not available, using basic csv")

try:
    import matplotlib
    matplotlib.use('Agg')  # Non-GUI backend
    import matplotlib.pyplot as plt
    import seaborn as sns
    PLOTTING_AVAILABLE = True
except ImportError:
    PLOTTING_AVAILABLE = False
    print("⚠️  Matplotlib/Seaborn not available, plots disabled")

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR = Path(__file__).parent.resolve()
RESULTS_CSV = BASE_DIR / "results" / "final_screening.csv"
STRUCTURES_DIR = BASE_DIR / "results" / "structures"
TOP_HITS_DIR = BASE_DIR / "results" / "top_hits_pdb"
FIGURES_DIR = BASE_DIR / "results" / "figures"
SUMMARY_CSV = BASE_DIR / "results" / "top_hits_summary.csv"

# Selection criteria
AFFINITY_THRESHOLD = -7.0  # kcal/mol (more negative = better)
TOP_N = 10                  # Number of candidates to copy
PERCENTILE = 1              # Alternative top 1%


# =============================================================================
# FUNCTIONS
# =============================================================================

def load_results_csv():
    """Load virtual screening results."""
    if not RESULTS_CSV.exists():
        print(f"❌ Error: File not found {RESULTS_CSV}")
        print("   Run first: python 13_run_virtual_screening.py")
        return None
    
    if PANDAS_AVAILABLE:
        df = pd.read_csv(RESULTS_CSV)
        return df
    else:
        # Load with basic csv
        results = []
        with open(RESULTS_CSV, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                results.append(row)
        return results


def filter_successful(data):
    """Filter only successful results with valid affinity."""
    if PANDAS_AVAILABLE:
        df = data[data['Status'] == 'Success'].copy()
        df['Binding_Affinity'] = pd.to_numeric(df['Binding_Affinity'], errors='coerce')
        df = df.dropna(subset=['Binding_Affinity'])
        return df
    else:
        return [r for r in data 
                if r.get('Status') == 'Success' 
                and r.get('Binding_Affinity', '')]


def select_top_hits(data, threshold=AFFINITY_THRESHOLD, top_n=TOP_N, use_percentile=False):
    """
    Select the best candidates.
    
    Args:
        data: DataFrame or list with results
        threshold: Affinity threshold (kcal/mol)
        top_n: Number of candidates to select
        use_percentile: If True, use percentile instead of fixed threshold
    """
    if PANDAS_AVAILABLE:
        df = data.sort_values('Binding_Affinity', ascending=True)
        
        if use_percentile:
            cutoff = df['Binding_Affinity'].quantile(PERCENTILE / 100)
            hits = df[df['Binding_Affinity'] <= cutoff]
            print(f"✓ Using percentile {PERCENTILE}%: cutoff = {cutoff:.2f} kcal/mol")
        else:
            hits = df[df['Binding_Affinity'] < threshold]
            print(f"✓ Using fixed threshold: < {threshold} kcal/mol")
        
        # Sort by affinity (best first)
        hits = hits.sort_values('Binding_Affinity', ascending=True)
        
        # Top N
        top = df.head(top_n)
        
        return hits, top
    else:
        # Sort manually
        sorted_data = sorted(data, key=lambda x: float(x['Binding_Affinity']))
        
        # Filter by threshold
        hits = [r for r in sorted_data if float(r['Binding_Affinity']) < threshold]
        top = sorted_data[:top_n]
        
        return hits, top


def copy_top_pdb_files(top_hits, max_files=TOP_N):
    """Copy PDB files from top hits to a separate directory."""
    TOP_HITS_DIR.mkdir(parents=True, exist_ok=True)
    
    copied = 0
    
    if PANDAS_AVAILABLE:
        for _, row in top_hits.head(max_files).iterrows():
            protein_id = row['Protein_ID']
            pdb_file = STRUCTURES_DIR / f"{protein_id}.pdb"
            
            if pdb_file.exists():
                dest = TOP_HITS_DIR / f"{copied+1:02d}_{protein_id}.pdb"
                shutil.copy2(pdb_file, dest)
                copied += 1
                print(f"  ✓ Copied: {protein_id} (Affinity: {row['Binding_Affinity']:.2f})")
            else:
                print(f"  ⚠️  Not found: {pdb_file}")
    else:
        for i, row in enumerate(top_hits[:max_files]):
            protein_id = row['Protein_ID']
            pdb_file = STRUCTURES_DIR / f"{protein_id}.pdb"
            
            if pdb_file.exists():
                dest = TOP_HITS_DIR / f"{i+1:02d}_{protein_id}.pdb"
                shutil.copy2(pdb_file, dest)
                copied += 1
    
    print(f"\n✓ {copied} structures copied to {TOP_HITS_DIR}")
    return copied


def save_summary(hits, top):
    """Save hits summary."""
    if PANDAS_AVAILABLE:
        hits.to_csv(SUMMARY_CSV, index=False)
        print(f"✓ Summary saved: {SUMMARY_CSV}")
    else:
        with open(SUMMARY_CSV, 'w', newline='') as f:
            if hits:
                writer = csv.DictWriter(f, fieldnames=hits[0].keys())
                writer.writeheader()
                writer.writerows(hits)


def plot_affinity_distribution(data, hits):
    """Generate a histogram of affinity distribution."""
    if not PLOTTING_AVAILABLE:
        print("⚠️  Plots not available (matplotlib not installed)")
        return
    
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Panel 1: Histograma con distribución
    ax1 = axes[0]
    affinities = data['Binding_Affinity'].values
    
    ax1.hist(affinities, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    ax1.axvline(x=AFFINITY_THRESHOLD, color='red', linestyle='--', 
                linewidth=2, label=f'Threshold: {AFFINITY_THRESHOLD} kcal/mol')
    
    # Mark hits zone
    hit_affinities = hits['Binding_Affinity'].values if len(hits) > 0 else []
    if len(hit_affinities) > 0:
        ax1.axvspan(min(affinities), AFFINITY_THRESHOLD, alpha=0.2, color='green', 
                    label=f'Hits: {len(hits)}')
    
    ax1.set_xlabel('Binding Affinity (kcal/mol)', fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    ax1.set_title('Distribution of Binding Affinities', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: BoxPlot
    ax2 = axes[1]
    
    # Categorize
    categories = ['All Candidates']
    box_data = [affinities]
    
    if len(hit_affinities) > 0:
        categories.append('Hits (< -7.0)')
        box_data.append(hit_affinities)
    
    bp = ax2.boxplot(box_data, labels=categories, patch_artist=True)
    
    colors = ['lightblue', 'lightgreen']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    
    ax2.axhline(y=AFFINITY_THRESHOLD, color='red', linestyle='--', 
                linewidth=2, label=f'Threshold')
    ax2.set_ylabel('Binding Affinity (kcal/mol)', fontsize=12)
    ax2.set_title('Affinity Comparison', fontsize=14)
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    output_path = FIGURES_DIR / "affinity_distribution.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Distribution plot: {output_path}")


def plot_top_hits_barplot(top_hits, max_show=20):
    """Generate a barplot of top hits."""
    if not PLOTTING_AVAILABLE:
        return
    
    if len(top_hits) == 0:
        return
    
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    
    # Take only the first N
    if PANDAS_AVAILABLE:
        plot_data = top_hits.head(max_show)
        protein_ids = plot_data['Protein_ID'].values
        affinities = plot_data['Binding_Affinity'].values
    else:
        plot_data = top_hits[:max_show]
        protein_ids = [r['Protein_ID'] for r in plot_data]
        affinities = [float(r['Binding_Affinity']) for r in plot_data]
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Colors based on affinity
    colors = ['darkgreen' if a < -8.0 else 'green' if a < -7.5 else 'limegreen' 
              for a in affinities]
    
    bars = ax.barh(range(len(protein_ids)), affinities, color=colors, edgecolor='black')
    
    ax.set_yticks(range(len(protein_ids)))
    ax.set_yticklabels(protein_ids, fontsize=9)
    ax.invert_yaxis()  # Best on top
    
    ax.axvline(x=AFFINITY_THRESHOLD, color='red', linestyle='--', 
               linewidth=2, label=f'Threshold: {AFFINITY_THRESHOLD}')
    
    ax.set_xlabel('Binding Affinity (kcal/mol)', fontsize=12)
    ax.set_title(f'Top {len(protein_ids)} Candidates - Virtual Screening', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='x')
    
    # Add values on bars
    for bar, aff in zip(bars, affinities):
        ax.text(aff - 0.3, bar.get_y() + bar.get_height()/2, 
                f'{aff:.2f}', va='center', ha='right', fontsize=8, color='white')
    
    plt.tight_layout()
    
    output_path = FIGURES_DIR / "top_hits_barplot.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Barplot de Top Hits: {output_path}")


def print_stats(data, hits, top):
    """Print detailed statistics."""
    print("\n" + "=" * 70)
    print("VIRTUAL SCREENING STATISTICS")
    print("=" * 70)
    
    if PANDAS_AVAILABLE:
        affinities = data['Binding_Affinity']
        
        print(f"\n📊 Affinity distribution:")
        print(f"  - Total proteins analyzed: {len(data)}")
        print(f"  - Average affinity: {affinities.mean():.2f} kcal/mol")
        print(f"  - Standard deviation: {affinities.std():.2f} kcal/mol")
        print(f"  - Minimum (best): {affinities.min():.2f} kcal/mol")
        print(f"  - Maximum (worst): {affinities.max():.2f} kcal/mol")
        print(f"  - Median: {affinities.median():.2f} kcal/mol")
        
        print(f"\n🎯 Identified hits (Affinity < {AFFINITY_THRESHOLD}):")
        print(f"  - Total hits: {len(hits)}")
        print(f"  - Percentage: {100*len(hits)/len(data):.2f}%")
        
        if len(hits) > 0:
            print(f"\n🏆 Top {TOP_N} best candidates:")
            for i, (_, row) in enumerate(top.head(TOP_N).iterrows(), 1):
                print(f"  {i:2d}. {row['Protein_ID']}: {row['Binding_Affinity']:.2f} kcal/mol")
    else:
        affinities = [float(r['Binding_Affinity']) for r in data]
        print(f"\n📊 Total proteins: {len(data)}")
        print(f"  - Best affinity: {min(affinities):.2f} kcal/mol")
        print(f"  - Hits (< {AFFINITY_THRESHOLD}): {len(hits)}")


def main():
    """Main function."""
    global AFFINITY_THRESHOLD, TOP_N
    
    # Arguments
    parser = argparse.ArgumentParser(description='Top Hits Selection')
    parser.add_argument('--threshold', type=float, default=AFFINITY_THRESHOLD,
                       help=f'Affinity threshold (default: {AFFINITY_THRESHOLD})')
    parser.add_argument('--top', type=int, default=TOP_N,
                       help=f'Number of top hits to copy (default: {TOP_N})')
    parser.add_argument('--percentile', action='store_true',
                       help='Use percentile instead of fixed threshold')
    args = parser.parse_args()
    
    AFFINITY_THRESHOLD = args.threshold
    TOP_N = args.top
    
    print("=" * 70)
    print("PHASE 4 - STEP 4: TOP CANDIDATE SELECTION")
    print("=" * 70)
    
    # Load results
    print(f"\n📂 Loading results from: {RESULTS_CSV}")
    data = load_results_csv()
    
    if data is None:
        return False
    
    if PANDAS_AVAILABLE:
        print(f"✓ Loaded {len(data)} records")
    else:
        print(f"✓ Loaded {len(data)} records")
    
    # Filter successful
    print("\n📋 Filtering successful results...")
    successful = filter_successful(data)
    
    if PANDAS_AVAILABLE:
        print(f"✓ Successful docking: {len(successful)}")
    else:
        print(f"✓ Successful docking: {len(successful)}")
    
    if (PANDAS_AVAILABLE and len(successful) == 0) or (not PANDAS_AVAILABLE and len(successful) == 0):
        print("❌ No successful results to analyze")
        return False
    
    # Select top hits
    print(f"\n🔍 Selecting candidates...")
    hits, top = select_top_hits(successful, args.threshold, args.top, args.percentile)
    
    # Statistics
    print_stats(successful, hits, top)
    
    # Copy PDB files
    print(f"\n📁 Copying Top {TOP_N} structures...")
    copy_top_pdb_files(top, TOP_N)
    
    # Save summary
    print(f"\n💾 Saving summary...")
    save_summary(hits, top)
    
    # Generate plots
    print(f"\n📊 Generating visualizations...")
    if PANDAS_AVAILABLE:
        plot_affinity_distribution(successful, hits)
        plot_top_hits_barplot(top)
    
    # Final summary
    print("\n" + "=" * 70)
    print("✅ TOP HITS SELECTION COMPLETED")
    print("=" * 70)
    
    print(f"\n📂 Generated files:")
    print(f"  - Summary: {SUMMARY_CSV}")
    print(f"  - Top PDBs: {TOP_HITS_DIR}/")
    if PLOTTING_AVAILABLE:
        print(f"  - Plots: {FIGURES_DIR}/")
    
    print(f"\n🎯 Next step:")
    print("  1. Review candidates in PyMOL or Chimera")
    print("  2. Validate binding site with structural analysis")
    print("  3. Consider synthesis and experimental characterization")
    
    return True


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
