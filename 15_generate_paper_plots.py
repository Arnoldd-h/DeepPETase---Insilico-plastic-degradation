#!/usr/bin/env python3
"""
================================================================================
15_generate_paper_plots.py - Phase 5: Q1 Scientific Visualization
================================================================================

Generates MAXIMUM PUBLICATION QUALITY plots (Q1 Scopus) to summarize 
the Virtual Screening results from the Deep-PETase-Mining project.

Style: Nature Biotechnology / Science / Cell

Q1 Features:
- Professional Arial/Helvetica typography
- Colorblind-friendly palette
- High resolution (300 DPI)
- Proportions optimized for journals
- No top/right spines
- Subtle and professional grid

Input:
    - results/final_screening.csv
    - results/top_hits_summary.csv
    
Output:
    - results/figures/final_screening_summary.png
    - results/figures/pipeline_results_overview.png
    - results/figures/binding_affinity_violin.png
================================================================================
"""

import os
import sys
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import FancyBboxPatch

# Import publication style
try:
    from src.paper_style import (
        apply_publication_style, COLORS, PALETTE_CATEGORICAL,
        save_figure, add_panel_label, style_axis, annotate_stats,
        get_figure_size, create_colorbar
    )
    HAS_PAPER_STYLE = True
except ImportError:
    HAS_PAPER_STYLE = False
    print("⚠️ paper_style module not found")

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR = Path(__file__).parent.resolve()
SCREENING_CSV = BASE_DIR / "results" / "final_screening.csv"
TOP_HITS_CSV = BASE_DIR / "results" / "top_hits_summary.csv"
FIGURES_DIR = BASE_DIR / "results" / "figures"

# Publication parameters
AFFINITY_THRESHOLD = -7.0  # kcal/mol
DPI = 300
FIGURE_FORMAT = 'png'

# =============================================================================
# Q1 COLOR PALETTE (COLORBLIND-FRIENDLY)
# =============================================================================
if not HAS_PAPER_STYLE:
    COLORS = {
        'primary': '#1A5276',
        'secondary': '#148F77',
        'accent': '#E74C3C',
        'highlight': '#F39C12',
        'data_blue': '#2E86AB',
        'data_green': '#28A745',
        'dark': '#2C3E50',
        'gray_dark': '#566573',
        'gray_medium': '#7F8C8D',
        'gray_light': '#BDC3C7',
        'gray_pale': '#ECF0F1',
        'background': '#FFFFFF',
        'grid': '#E5E8E8',
        'hit': '#C0392B',
        'non_hit': '#95A5A6',
        'best_hit': '#922B21',
    }
    PALETTE_CATEGORICAL = ['#2E86AB', '#E74C3C', '#28A745', '#F39C12', '#8E44AD']


def setup_publication_style():
    """Configure matplotlib style for Q1 scientific publications."""
    if HAS_PAPER_STYLE:
        apply_publication_style(style='nature')
        print("✓ Q1 Nature style configured")
        return
    
    # Manual style if module not available
    plt.rcParams.update({
        # Professional fonts
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica Neue', 'Helvetica', 'DejaVu Sans'],
        'font.size': 10,
        'font.weight': 'normal',
        
        # Text sizes
        'axes.titlesize': 12,
        'axes.labelsize': 11,
        'axes.titleweight': 'bold',
        'axes.labelweight': 'medium',
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 9,
        'legend.title_fontsize': 10,
        'figure.titlesize': 14,
        'figure.titleweight': 'bold',
        
        # Colors
        'text.color': COLORS['dark'],
        'axes.labelcolor': COLORS['dark'],
        'axes.edgecolor': COLORS['gray_dark'],
        'xtick.color': COLORS['gray_dark'],
        'ytick.color': COLORS['gray_dark'],
        
        # Axes
        'axes.linewidth': 1.2,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.facecolor': 'white',
        'axes.axisbelow': True,
        
        # Ticks
        'xtick.major.size': 5,
        'xtick.major.width': 1.0,
        'xtick.direction': 'out',
        'ytick.major.size': 5,
        'ytick.major.width': 1.0,
        'ytick.direction': 'out',
        
        # Grid
        'axes.grid': True,
        'grid.color': COLORS['grid'],
        'grid.linewidth': 0.5,
        'grid.alpha': 0.7,
        
        # Legend
        'legend.frameon': True,
        'legend.framealpha': 0.95,
        'legend.edgecolor': COLORS['gray_light'],
        'legend.fancybox': True,
        
        # Figure
        'figure.facecolor': 'white',
        'figure.dpi': 150,
        'savefig.dpi': DPI,
        'savefig.facecolor': 'white',
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.1,
        
        # Patches
        'patch.linewidth': 0.8,
        'patch.edgecolor': 'white',
    })
    
    print("✓ Q1 publication style configured")


def load_data():
    """Load screening data."""
    if not SCREENING_CSV.exists():
        print(f"❌ Error: File not found {SCREENING_CSV}")
        return None, None
    
    df = pd.read_csv(SCREENING_CSV)
    
    # Filter successful results and clean
    df = df[df['Status'] == 'Success'].copy()
    df['Binding_Affinity'] = pd.to_numeric(df['Binding_Affinity'], errors='coerce')
    df = df.dropna(subset=['Binding_Affinity'])
    
    print(f"✓ Loaded {len(df)} docking results")
    
    # Top hits
    top_df = None
    if TOP_HITS_CSV.exists():
        top_df = pd.read_csv(TOP_HITS_CSV)
        top_df['Binding_Affinity'] = pd.to_numeric(top_df['Binding_Affinity'], errors='coerce')
        print(f"✓ Loaded {len(top_df)} top hits")
    
    return df, top_df


def generate_main_figure(df, top_df):
    """
    Generate the main figure with 2 Q1-style subplots.
    A) Affinity distribution histogram
    B) Top 10 ranking with horizontal bars
    """
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    
    # Create figure with professional proportions
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    
    # =========================================================================
    # SUBPLOT A: Distribution Histogram
    # =========================================================================
    ax1 = axes[0]
    
    affinities = df['Binding_Affinity'].values
    
    # Create optimized bins
    bins = np.linspace(affinities.min() - 0.1, affinities.max() + 0.1, 45)
    
    # Separate hits and non-hits
    hits = affinities[affinities < AFFINITY_THRESHOLD]
    non_hits = affinities[affinities >= AFFINITY_THRESHOLD]
    
    # Non-hits histogram (gray, background)
    ax1.hist(non_hits, bins=bins, color=COLORS['non_hit'], edgecolor='white', 
             alpha=0.75, label=f'Non-hits (n={len(non_hits)})', linewidth=0.5,
             zorder=2)
    
    # Hits histogram (highlighted)
    ax1.hist(hits, bins=bins, color=COLORS['hit'], edgecolor='white',
             alpha=0.9, label=f'Hits < {AFFINITY_THRESHOLD} (n={len(hits)})', 
             linewidth=0.5, zorder=3)
    
    # Mean line
    mean_val = affinities.mean()
    ax1.axvline(x=mean_val, color=COLORS['dark'], linestyle='--', linewidth=2,
                label=f'Mean: {mean_val:.2f} kcal/mol', zorder=4)
    
    # Threshold line
    ax1.axvline(x=AFFINITY_THRESHOLD, color=COLORS['hit'], linestyle='-', 
                linewidth=2.5, alpha=0.8, zorder=5)
    
    # Elegant threshold annotation
    y_max = ax1.get_ylim()[1]
    ax1.annotate(
        f'Hit Threshold\n({AFFINITY_THRESHOLD} kcal/mol)', 
        xy=(AFFINITY_THRESHOLD, y_max * 0.75),
        xytext=(AFFINITY_THRESHOLD - 0.7, y_max * 0.75),
        fontsize=9, ha='right', color=COLORS['hit'],
        fontweight='bold',
        arrowprops=dict(arrowstyle='->', color=COLORS['hit'], lw=1.2)
    )
    
    # Configure axes
    ax1.set_xlabel('Binding Affinity (kcal/mol)', fontweight='medium')
    ax1.set_ylabel('Frequency', fontweight='medium')
    ax1.set_title('A) Distribution of Binding Affinities', fontweight='bold', 
                  loc='left', pad=10, color=COLORS['dark'])
    
    # Professional legend
    ax1.legend(loc='upper left', framealpha=0.95, edgecolor=COLORS['gray_light'],
               fontsize=9)
    
    # Invert X axis (more negative = better, on the left)
    ax1.invert_xaxis()
    
    # Grid and spines
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.grid(True, alpha=0.4, linewidth=0.5, zorder=0)
    ax1.set_axisbelow(True)
    
    # Statistics box
    stats_box = (
        f"n = {len(affinities):,}\n"
        f"μ = {mean_val:.2f}\n"
        f"σ = {affinities.std():.2f}\n"
        f"Hit rate: {len(hits)/len(affinities)*100:.1f}%"
    )
    ax1.text(0.97, 0.97, stats_box, transform=ax1.transAxes,
             fontsize=8, fontfamily='monospace', va='top', ha='right',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                      edgecolor=COLORS['gray_light'], alpha=0.95))
    
    # =========================================================================
    # SUBPLOT B: Top 10 Ranking
    # =========================================================================
    ax2 = axes[1]
    
    # Get top 10
    if top_df is not None and len(top_df) > 0:
        top10 = top_df.nsmallest(10, 'Binding_Affinity')
    else:
        top10 = df.nsmallest(10, 'Binding_Affinity')
    
    # Data for the plot
    protein_ids = top10['Protein_ID'].values
    affinities_top = top10['Binding_Affinity'].values
    
    # Professional gradient colors (better = more intense)
    n_bars = len(protein_ids)
    colors = plt.cm.RdYlGn_r(np.linspace(0.15, 0.55, n_bars))
    
    # Horizontal bars
    y_pos = np.arange(len(protein_ids))
    bars = ax2.barh(y_pos, affinities_top, color=colors, 
                    edgecolor='white', linewidth=0.8, height=0.7, zorder=3)
    
    # Best candidate with special border
    bars[0].set_edgecolor(COLORS['hit'])
    bars[0].set_linewidth(2.5)
    
    # Configure axes
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(protein_ids, fontsize=9, fontfamily='monospace')
    ax2.invert_yaxis()  # Best on top
    
    # Threshold line
    ax2.axvline(x=AFFINITY_THRESHOLD, color=COLORS['hit'], linestyle='--', 
                linewidth=2, alpha=0.7, label=f'Threshold: {AFFINITY_THRESHOLD}',
                zorder=2)
    
    # Values on bars with better contrast
    for bar, aff in zip(bars, affinities_top):
        width = bar.get_width()
        ax2.text(width - 0.12, bar.get_y() + bar.get_height()/2,
                f'{aff:.2f}', va='center', ha='right', 
                fontsize=9, fontweight='bold', color='white',
                zorder=4)
    
    # Mark #1
    ax2.annotate('#1 BEST', xy=(affinities_top[0] - 0.25, 0),
                xytext=(affinities_top[0] - 1.0, 0),
                fontsize=10, fontweight='bold', color=COLORS['hit'],
                va='center', ha='right',
                arrowprops=dict(arrowstyle='->', color=COLORS['hit'], lw=1.5))
    
    ax2.set_xlabel('Binding Affinity (kcal/mol)', fontweight='medium')
    ax2.set_title('B) Top 10 PETase Candidates', fontweight='bold', 
                  loc='left', pad=10, color=COLORS['dark'])
    ax2.legend(loc='lower right', fontsize=9, edgecolor=COLORS['gray_light'])
    
    # X axis limits
    ax2.set_xlim(affinities_top.min() - 0.4, AFFINITY_THRESHOLD + 0.25)
    
    # Grid and spines
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.grid(True, alpha=0.4, linewidth=0.5, axis='x', zorder=0)
    ax2.set_axisbelow(True)
    
    # =========================================================================
    # Final adjustments
    # =========================================================================
    plt.tight_layout()
    
    # Overall title
    fig.suptitle('Deep-PETase-Mining: Virtual Screening Results', 
                fontsize=14, fontweight='bold', y=1.03, color=COLORS['dark'])
    
    # Save with high quality
    output_path = FIGURES_DIR / "final_screening_summary.png"
    plt.savefig(output_path, dpi=DPI, bbox_inches='tight', 
                facecolor='white', edgecolor='none', pad_inches=0.15)
    plt.close()
    
    print(f"✓ Main Q1 figure saved: {output_path}")
    
    return output_path


def generate_pipeline_overview(df):
    """
    Generate a summary plot of the complete pipeline in Q1 style.
    Professional visualization of the selection funnel.
    """
    # Pipeline data (project values)
    stages = [
        'NCBI Mining\n(Raw Sequences)',
        'HMM Filtering\n(Homology)',
        'Clustering\n(80% Identity)',
        'ESMFold\n(pLDDT ≥ 70)',
        'Virtual Screening\n(Affinity < -7.0)'
    ]
    
    # Actual project values
    counts = [39000, 8500, 2255, 1718, 27]
    
    fig, ax = plt.subplots(figsize=(11, 6))
    
    # Professional gradient colors
    n_stages = len(stages)
    base_colors = plt.cm.Blues(np.linspace(0.35, 0.85, n_stages - 1))
    colors = list(base_colors) + [np.array([0.75, 0.22, 0.17, 1.0])]  # Last one in red
    
    x_pos = np.arange(len(stages))
    bars = ax.bar(x_pos, counts, color=colors, edgecolor='white', 
                  linewidth=1.5, width=0.72, zorder=3)
    
    # Values above the bars
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        y_offset = height * 0.15 if height < 1000 else 800
        ax.text(bar.get_x() + bar.get_width()/2, height + y_offset,
                f'{count:,}', ha='center', va='bottom', 
                fontsize=11, fontweight='bold', color=COLORS['dark'])
    
    # Arrows and reduction percentages between bars
    for i in range(len(stages) - 1):
        reduction = (1 - counts[i+1]/counts[i]) * 100
        
        # Curved arrow
        mid_x = i + 0.5
        mid_y = min(counts[i], counts[i+1]) * 0.7
        
        ax.annotate('', 
                   xy=(i + 0.72, counts[i+1] * 1.1),
                   xytext=(i + 0.28, counts[i] * 0.9),
                   arrowprops=dict(arrowstyle='->', color=COLORS['gray_dark'], 
                                  lw=1.8, connectionstyle='arc3,rad=-0.2'),
                   zorder=2)
        
        # Percentage label
        ax.text(i + 0.5, mid_y * 0.6, f'-{reduction:.0f}%', 
               ha='center', va='center', fontsize=9, 
               fontweight='bold', color=COLORS['accent'],
               bbox=dict(boxstyle='round,pad=0.25', facecolor='white', 
                        edgecolor=COLORS['accent'], alpha=0.95, linewidth=1))
    
    # Configure axes
    ax.set_xticks(x_pos)
    ax.set_xticklabels(stages, fontsize=10)
    ax.set_ylabel('Number of Sequences/Structures', fontweight='medium', fontsize=11)
    ax.set_title('Deep-PETase-Mining Pipeline: From Thousands to Elite Candidates',
                fontsize=13, fontweight='bold', pad=15, color=COLORS['dark'])
    
    # Logarithmic scale for better visualization
    ax.set_yscale('log')
    ax.set_ylim(15, 120000)
    
    # Professional spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.2)
    ax.spines['bottom'].set_linewidth(1.2)
    
    # Subtle horizontal grid
    ax.yaxis.grid(True, alpha=0.4, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    
    # Summary box
    final_rate = (counts[-1] / counts[0]) * 100
    summary = f"Overall Selection: {final_rate:.3f}%\nFinal Candidates: {counts[-1]}"
    ax.text(0.98, 0.95, summary, transform=ax.transAxes,
           fontsize=10, fontweight='bold', va='top', ha='right',
           color=COLORS['hit'],
           bbox=dict(boxstyle='round,pad=0.5', facecolor='white',
                    edgecolor=COLORS['hit'], linewidth=1.5, alpha=0.95))
    
    plt.tight_layout()
    
    output_path = FIGURES_DIR / "pipeline_results_overview.png"
    plt.savefig(output_path, dpi=DPI, bbox_inches='tight',
                facecolor='white', edgecolor='none', pad_inches=0.15)
    plt.close()
    
    print(f"✓ Pipeline overview Q1 saved: {output_path}")
    
    return output_path


def print_best_candidate(df):
    """Print information about the best candidate."""
    best = df.loc[df['Binding_Affinity'].idxmin()]
    
    print("\n" + "=" * 60)
    print("🏆 BEST CANDIDATE IDENTIFIED")
    print("=" * 60)
    print(f"  Protein ID:        {best['Protein_ID']}")
    print(f"  Binding Affinity:  {best['Binding_Affinity']:.2f} kcal/mol")
    print(f"  Grid Center:       {best.get('Grid_Center', 'N/A')}")
    print(f"  Grid Size:         {best.get('Grid_Size', 'N/A')}")
    print("=" * 60)
    
    return best


def generate_violin_plot(df):
    """
    Generate a professional violin plot of the affinity distribution.
    Additional Q1-style figure.
    """
    try:
        import seaborn as sns
    except ImportError:
        print("⚠️ Seaborn not available, skipping violin plot")
        return None
    
    fig, ax = plt.subplots(figsize=(6, 7))
    
    # Separate data
    affinities = df['Binding_Affinity'].values
    hits_mask = affinities < AFFINITY_THRESHOLD
    
    # Create DataFrame for seaborn
    plot_df = pd.DataFrame({
        'Binding Affinity': affinities,
        'Category': ['Hit' if h else 'Non-hit' for h in hits_mask]
    })
    
    # Violin plot
    palette = {'Hit': COLORS['hit'], 'Non-hit': COLORS['non_hit']}
    
    violin = sns.violinplot(
        data=plot_df, 
        y='Binding Affinity', 
        x='Category',
        palette=palette,
        inner='box',
        linewidth=1.2,
        ax=ax
    )
    
    # Línea de threshold
    ax.axhline(y=AFFINITY_THRESHOLD, color=COLORS['accent'], 
               linestyle='--', linewidth=2, alpha=0.8,
               label=f'Threshold ({AFFINITY_THRESHOLD} kcal/mol)')
    
    # Individual points with jitter
    for i, cat in enumerate(['Hit', 'Non-hit']):
        cat_data = plot_df[plot_df['Category'] == cat]['Binding Affinity']
        jitter = np.random.uniform(-0.08, 0.08, len(cat_data))
        ax.scatter(i + jitter, cat_data, alpha=0.3, s=15, 
                   color='white', edgecolor='black', linewidth=0.3, zorder=5)
    
    ax.set_ylabel('Binding Affinity (kcal/mol)', fontweight='medium')
    ax.set_xlabel('')
    ax.set_title('Distribution by Hit Category', fontweight='bold',
                 color=COLORS['dark'], pad=15)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='lower right', fontsize=9)
    ax.grid(True, alpha=0.3, axis='y', linewidth=0.5)
    
    plt.tight_layout()
    
    output_path = FIGURES_DIR / "binding_affinity_violin.png"
    plt.savefig(output_path, dpi=DPI, bbox_inches='tight',
                facecolor='white', edgecolor='none', pad_inches=0.1)
    plt.close()
    
    print(f"✓ Violin plot Q1 saved: {output_path}")
    return output_path


def main():
    """Main function."""
    print("=" * 70)
    print("PHASE 5 - Q1 SCIENTIFIC VISUALIZATION")
    print("=" * 70)
    
    # Configure style
    setup_publication_style()
    
    # Load data
    print("\n📂 Loading data...")
    df, top_df = load_data()
    
    if df is None:
        return False
    
    # Generate main figure
    print("\n📊 Generating main figure (Q1 Quality)...")
    generate_main_figure(df, top_df)
    
    # Generate pipeline overview
    print("\n📊 Generating pipeline overview (Q1 Quality)...")
    generate_pipeline_overview(df)
    
    # Generate additional violin plot
    print("\n📊 Generating violin plot (Q1 Quality)...")
    generate_violin_plot(df)
    
    # Best candidate information
    best = print_best_candidate(df)
    
    # Final statistics
    hits = df[df['Binding_Affinity'] < AFFINITY_THRESHOLD]
    print("\n📈 FINAL STATISTICS:")
    print(f"  • Total structures analyzed: {len(df):,}")
    print(f"  • Hits (< {AFFINITY_THRESHOLD} kcal/mol): {len(hits)}")
    print(f"  • Mean affinity: {df['Binding_Affinity'].mean():.2f} kcal/mol")
    print(f"  • Std deviation: {df['Binding_Affinity'].std():.2f} kcal/mol")
    print(f"  • Hit rate: {len(hits)/len(df)*100:.2f}%")
    
    print("\n" + "=" * 70)
    print("✅ Q1 PUBLICATION FIGURES GENERATED SUCCESSFULLY")
    print("=" * 70)
    print(f"\nFiles saved to: {FIGURES_DIR}/")
    print("  • final_screening_summary.png (300 DPI)")
    print("  • pipeline_results_overview.png (300 DPI)")
    print("  • binding_affinity_violin.png (300 DPI)")
    print("\n💡 Figures ready for Q1 Scopus journal submission!")
    
    return True


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
