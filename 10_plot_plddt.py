#!/usr/bin/env python3
"""
================================================================================
10_plot_plddt.py - Phase 3: pLDDT Visualization (Q1 Quality)
================================================================================

Generates PUBLICATION QUALITY Q1 plots for the distribution of 
pLDDT scores from ESMFold predicted structures.

Style: Nature Biotechnology / Cell (Q1 Scopus)

Input: results/logs/structure_quality.csv
Output: 
    - results/figures/plddt_distribution.png
    - results/figures/phase3_summary.png
================================================================================
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator
import numpy as np
from pathlib import Path

# Importar estilo de publicación
try:
    from src.paper_style import (
        apply_publication_style, COLORS, save_figure, 
        add_panel_label, style_axis, annotate_stats
    )
    HAS_PAPER_STYLE = True
except ImportError:
    HAS_PAPER_STYLE = False
    print("⚠️ paper_style module not found, using default style")

# =============================================================================
# CONFIGURATION
# =============================================================================
INPUT_CSV = Path("results/logs/structure_quality.csv")
OUTPUT_FIGURE = Path("results/figures/plddt_distribution.png")
PLDDT_THRESHOLD = 70.0

# =============================================================================
# Q1 COLOR PALETTE
# =============================================================================
if not HAS_PAPER_STYLE:
    COLORS = {
        'primary': '#1A5276',
        'secondary': '#148F77',
        'accent': '#E74C3C',
        'data_green': '#28A745',
        'dark': '#2C3E50',
        'gray_dark': '#566573',
        'gray_light': '#BDC3C7',
        'background': '#FFFFFF',
        'grid': '#E5E8E8',
        'quality_high': '#27AE60',
        'quality_medium': '#F1C40F',
        'quality_low': '#E67E22',
        'quality_very_low': '#E74C3C',
    }

# Colors for pLDDT categories (colorblind-friendly)
PLDDT_COLORS = {
    'very_high': '#1E8449',   # Dark green - Very high confidence (>90)
    'high': '#27AE60',        # Green - High confidence (70-90)
    'low': '#F39C12',         # Orange - Low confidence (50-70)
    'very_low': '#C0392B'     # Red - Very low confidence (<50)
}

FIGURE_DPI = 300

# =============================================================================
# Q1 STYLE CONFIGURATION
# =============================================================================

def setup_q1_style():
    """Configures Q1 Scopus publication style."""
    if HAS_PAPER_STYLE:
        apply_publication_style(style='nature')
    else:
        plt.rcParams.update({
            'font.family': 'sans-serif',
            'font.sans-serif': ['Arial', 'Helvetica Neue', 'DejaVu Sans'],
            'font.size': 10,
            'axes.titlesize': 12,
            'axes.labelsize': 11,
            'axes.titleweight': 'bold',
            'axes.labelweight': 'medium',
            'axes.linewidth': 1.2,
            'axes.spines.top': False,
            'axes.spines.right': False,
            'axes.edgecolor': '#2C3E50',
            'axes.labelcolor': '#2C3E50',
            'axes.facecolor': 'white',
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'legend.fontsize': 9,
            'legend.framealpha': 0.95,
            'legend.edgecolor': '#BDC3C7',
            'figure.facecolor': 'white',
            'figure.dpi': 150,
            'savefig.dpi': 300,
            'savefig.bbox': 'tight',
            'grid.alpha': 0.5,
            'grid.linewidth': 0.5,
        })

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def ensure_directories():
    """Creates necessary directories if they don't exist."""
    OUTPUT_FIGURE.parent.mkdir(parents=True, exist_ok=True)


def main():
    """Main function of the script."""
    print("=" * 70)
    print("PHASE 3 - STEP 4: pLDDT QUALITY VISUALIZATION (Q1 Style)")
    print("=" * 70)
    
    # Configure Q1 style
    setup_q1_style()
    print("✓ Q1 Scopus style configured")
    
    # Create directories
    ensure_directories()
    
    # Verify input file
    if not INPUT_CSV.exists():
        raise FileNotFoundError(
            f"File not found: {INPUT_CSV}\n"
            "Did you run 09_filter_structures_plddt.py first?"
        )
    
    # Load data
    print(f"\nLoading data from: {INPUT_CSV}")
    df = pd.read_csv(INPUT_CSV)
    
    # Filter valid pLDDT values
    df_valid = df[df["Average_pLDDT"] != "N/A"].copy()
    df_valid["Average_pLDDT"] = pd.to_numeric(df_valid["Average_pLDDT"])
    
    if df_valid.empty:
        print("⚠️ No valid pLDDT data to visualize")
        return
    
    print(f"Structures with valid pLDDT: {len(df_valid)}")
    
    # Statistics
    plddt_values = df_valid["Average_pLDDT"]
    passed = (plddt_values >= PLDDT_THRESHOLD).sum()
    failed = (plddt_values < PLDDT_THRESHOLD).sum()
    
    print(f"\npLDDT Statistics:")
    print(f"  - Mean: {plddt_values.mean():.2f}")
    print(f"  - Median: {plddt_values.median():.2f}")
    print(f"  - Standard deviation: {plddt_values.std():.2f}")
    print(f"  - Minimum: {plddt_values.min():.2f}")
    print(f"  - Maximum: {plddt_values.max():.2f}")
    print(f"  - Passed filter (>= {PLDDT_THRESHOLD}): {passed}")
    print(f"  - Failed filter: {failed}"
    
    # =========================================================================
    # CREATE MAIN FIGURE (Q1 Quality)
    # =========================================================================
    fig, ax = plt.subplots(figsize=(9, 6))
    
    # Define bins for histogram
    bins = np.arange(0, 102, 2.5)  # Finer bins for better resolution
    
    # Create histogram
    n, bin_edges, patches = ax.hist(
        plddt_values, 
        bins=bins, 
        edgecolor='white',
        linewidth=0.6,
        alpha=0.9,
        zorder=3
    )
    
    # Color bars according to threshold with professional gradient
    for i, patch in enumerate(patches):
        bin_center = (bin_edges[i] + bin_edges[i+1]) / 2
        if bin_center >= 90:
            patch.set_facecolor(PLDDT_COLORS['very_high'])
        elif bin_center >= 70:
            patch.set_facecolor(PLDDT_COLORS['high'])
        elif bin_center >= 50:
            patch.set_facecolor(PLDDT_COLORS['low'])
        else:
            patch.set_facecolor(PLDDT_COLORS['very_low'])
    
    # Vertical threshold line with elegant style
    ax.axvline(
        x=PLDDT_THRESHOLD, 
        color=COLORS['accent'], 
        linestyle='--', 
        linewidth=2.5,
        zorder=5,
        alpha=0.9
    )
    
    # Threshold annotation
    ax.annotate(
        f'Quality Threshold\n(pLDDT = {PLDDT_THRESHOLD})',
        xy=(PLDDT_THRESHOLD, ax.get_ylim()[1] * 0.85),
        xytext=(PLDDT_THRESHOLD - 12, ax.get_ylim()[1] * 0.85),
        fontsize=9,
        fontweight='bold',
        color=COLORS['accent'],
        ha='right',
        va='center',
        arrowprops=dict(
            arrowstyle='->', 
            color=COLORS['accent'],
            lw=1.5
        )
    )
    
    # Statistics in professional box
    stats_text = (
        f"Total structures: {len(plddt_values):,}\n"
        f"Mean: {plddt_values.mean():.1f}\n"
        f"Median: {plddt_values.median():.1f}\n"
        f"Passed: {passed:,} ({passed/len(plddt_values)*100:.1f}%)"
    )
    ax.text(
        0.03, 0.97, 
        stats_text,
        transform=ax.transAxes,
        fontsize=9,
        fontfamily='monospace',
        verticalalignment='top',
        bbox=dict(
            boxstyle='round,pad=0.5',
            facecolor='white',
            edgecolor=COLORS['gray_light'],
            alpha=0.95,
            linewidth=0.8
        ),
        zorder=6
    )
    
    # Configure axes with Q1 style
    ax.set_xlabel('Average pLDDT Score', fontsize=11, fontweight='medium')
    ax.set_ylabel('Number of Structures', fontsize=11, fontweight='medium')
    ax.set_title(
        'Distribution of ESMFold Prediction Confidence\n'
        'Deep-PETase-Mining Structure Prediction',
        fontsize=12, 
        fontweight='bold',
        pad=15,
        color=COLORS['dark']
    )
    
    # Professional limits and ticks
    ax.set_xlim(0, 100)
    ax.set_xticks(range(0, 101, 10))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=8))
    
    # Elegant spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.2)
    ax.spines['bottom'].set_linewidth(1.2)
    
    # Subtle grid
    ax.grid(True, alpha=0.4, linewidth=0.5, linestyle='-', zorder=0)
    ax.set_axisbelow(True)
    
    # Professional legend
    legend_elements = [
        mpatches.Patch(color=PLDDT_COLORS['very_high'], 
                       label='Very High Confidence (≥90)'),
        mpatches.Patch(color=PLDDT_COLORS['high'], 
                       label='High Confidence (70-90)'),
        mpatches.Patch(color=PLDDT_COLORS['low'], 
                       label='Low Confidence (50-70)'),
        mpatches.Patch(color=PLDDT_COLORS['very_low'], 
                       label='Very Low Confidence (<50)'),
    ]
    ax.legend(
        handles=legend_elements,
        loc='upper right',
        fontsize=8,
        framealpha=0.95,
        edgecolor=COLORS['gray_light'],
        title='pLDDT Category',
        title_fontsize=9
    )
    
    # Save figure with high quality
    plt.tight_layout()
    print(f"\nSaving Q1 figure to: {OUTPUT_FIGURE}")
    plt.savefig(OUTPUT_FIGURE, dpi=FIGURE_DPI, bbox_inches='tight', 
                facecolor='white', edgecolor='none', pad_inches=0.1)
    plt.close()
    
    # Create Phase 3 summary plot
    create_phase3_summary_plot(df_valid)
    
    print("\n" + "=" * 70)
    print("Q1 VISUALIZATION COMPLETED")
    print("=" * 70)
    print(f"\nGenerated files (300 DPI):")
    print(f"  • {OUTPUT_FIGURE}")
    print(f"  • {OUTPUT_FIGURE.parent / 'phase3_summary.png'}")
    print("\n" + "=" * 70)


def create_phase3_summary_plot(df):
    """
    Creates a professional summary plot of the entire Phase 3.
    Q1 Scopus style with multiple panels.
    """
    output_file = Path("results/figures/phase3_summary.png")
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Phase 3: Structure Prediction Quality Summary', 
                 fontsize=14, fontweight='bold', y=0.98, color=COLORS['dark'])
    
    plddt_values = df["Average_pLDDT"]
    
    # =========================================================================
    # Panel A: pLDDT Histogram
    # =========================================================================
    ax1 = axes[0, 0]
    bins = np.arange(0, 102, 5)
    n, bin_edges, patches = ax1.hist(plddt_values, bins=bins, 
                                      color=COLORS['primary'], 
                                      edgecolor='white', 
                                      alpha=0.85, linewidth=0.6)
    
    # Color by category
    for i, patch in enumerate(patches):
        bin_center = (bin_edges[i] + bin_edges[i+1]) / 2
        if bin_center >= 90:
            patch.set_facecolor(PLDDT_COLORS['very_high'])
        elif bin_center >= 70:
            patch.set_facecolor(PLDDT_COLORS['high'])
        elif bin_center >= 50:
            patch.set_facecolor(PLDDT_COLORS['low'])
        else:
            patch.set_facecolor(PLDDT_COLORS['very_low'])
    
    ax1.axvline(x=PLDDT_THRESHOLD, color=COLORS['accent'], 
                linestyle='--', linewidth=2)
    ax1.set_xlabel('Average pLDDT', fontweight='medium')
    ax1.set_ylabel('Frequency', fontweight='medium')
    ax1.set_title('A) pLDDT Distribution', fontweight='bold', loc='left')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.grid(True, alpha=0.3, linewidth=0.5)
    
    # =========================================================================
    # Panel B: Category Pie Chart
    # =========================================================================
    ax2 = axes[0, 1]
    category_counts = df['Category'].value_counts()
    
    colors_map = {
        'Very High': PLDDT_COLORS['very_high'],
        'High': PLDDT_COLORS['high'], 
        'Low': PLDDT_COLORS['low'],
        'Very Low': PLDDT_COLORS['very_low']
    }
    colors = [colors_map.get(cat, COLORS['gray_light']) for cat in category_counts.index]
    
    wedges, texts, autotexts = ax2.pie(
        category_counts.values, 
        labels=None,  # Sin labels externos
        colors=colors,
        autopct='%1.1f%%',
        startangle=90,
        explode=[0.03] * len(category_counts),
        pctdistance=0.75,
        wedgeprops=dict(linewidth=2, edgecolor='white')
    )
    
    # Percentage format
    for autotext in autotexts:
        autotext.set_fontsize(9)
        autotext.set_fontweight('bold')
        autotext.set_color('white')
    
    ax2.legend(wedges, category_counts.index, title='Category',
               loc='center left', bbox_to_anchor=(1, 0.5),
               fontsize=9, title_fontsize=10)
    ax2.set_title('B) Quality Category Distribution', fontweight='bold', loc='left')
    
    # =========================================================================
    # Panel C: pLDDT Boxplot by Status
    # =========================================================================
    ax3 = axes[1, 0]
    status_data = []
    labels = []
    box_colors = []
    
    if 'Status' in df.columns:
        passed_data = df[df['Status'] == 'Passed']['Average_pLDDT'].values
        failed_data = df[df['Status'] == 'Failed']['Average_pLDDT'].values
        
        if len(passed_data) > 0:
            status_data.append(passed_data)
            labels.append(f'Passed\n(n={len(passed_data):,})')
            box_colors.append(PLDDT_COLORS['high'])
        if len(failed_data) > 0:
            status_data.append(failed_data)
            labels.append(f'Failed\n(n={len(failed_data):,})')
            box_colors.append(PLDDT_COLORS['very_low'])
    
    if status_data:
        bp = ax3.boxplot(status_data, tick_labels=labels, patch_artist=True,
                         widths=0.6, notch=False,
                         flierprops=dict(marker='o', markersize=4, alpha=0.5))
        
        for patch, color in zip(bp['boxes'], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
            patch.set_linewidth(1.5)
        
        for median in bp['medians']:
            median.set_color('white')
            median.set_linewidth(2)
        
        for whisker in bp['whiskers']:
            whisker.set_linewidth(1.2)
        for cap in bp['caps']:
            cap.set_linewidth(1.2)
    
    ax3.axhline(y=PLDDT_THRESHOLD, color=COLORS['accent'], 
                linestyle='--', linewidth=2, alpha=0.8,
                label=f'Threshold ({PLDDT_THRESHOLD})')
    ax3.set_ylabel('Average pLDDT', fontweight='medium')
    ax3.set_title('C) pLDDT by Filter Status', fontweight='bold', loc='left')
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.legend(loc='lower right', fontsize=8)
    ax3.grid(True, alpha=0.3, axis='y', linewidth=0.5)
    
    # =========================================================================
    # Panel D: Length vs pLDDT Scatter
    # =========================================================================
    ax4 = axes[1, 1]
    scatter = ax4.scatter(
        df['Num_Residues'], 
        df['Average_pLDDT'],
        c=df['Average_pLDDT'],
        cmap='RdYlGn',
        alpha=0.6,
        s=30,
        edgecolors='white',
        linewidth=0.3,
        vmin=40, vmax=95
    )
    
    ax4.axhline(y=PLDDT_THRESHOLD, color=COLORS['accent'], 
                linestyle='--', linewidth=2, alpha=0.8,
                label=f'Threshold ({PLDDT_THRESHOLD})')
    
    ax4.set_xlabel('Sequence Length (residues)', fontweight='medium')
    ax4.set_ylabel('Average pLDDT', fontweight='medium')
    ax4.set_title('D) Sequence Length vs Prediction Quality', fontweight='bold', loc='left')
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.legend(loc='lower right', fontsize=8)
    ax4.grid(True, alpha=0.3, linewidth=0.5)
    
    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax4, shrink=0.8, aspect=20, pad=0.02)
    cbar.set_label('pLDDT', fontsize=9, fontweight='medium')
    cbar.ax.tick_params(labelsize=8)
    cbar.outline.set_linewidth(0.8)
    
    # Save without tight_layout to avoid conflict with colorbar
    plt.savefig(output_file, dpi=FIGURE_DPI, bbox_inches='tight', 
                facecolor='white', edgecolor='none', pad_inches=0.15)
    plt.close()
    
    print(f"✓ Phase 3 summary (Q1) saved: {output_file}")


if __name__ == "__main__":
    main()
