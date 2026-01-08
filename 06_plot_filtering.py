#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
06_plot_filtering.py - Filtering Process Visualization (Q1 Quality)
================================================================================

Description:
    This script generates Q1 PUBLICATION QUALITY visualizations that 
    document the Phase 2 filtering process. 

    Generated plots:
    - Identity distribution with cutoff lines
    - Filtering funnel diagram
    - Score vs identity scatter plot

Date: January 2026
Project: Deep-PETase-Mining - Phase 2

Style: Nature Biotechnology / Cell style (Q1 Scopus)

Input:
    - results/logs/homology_scores.csv
    - results/logs/motif_analysis.csv (optional)

Output:
    - results/figures/identity_distribution.png
    - results/figures/filtering_funnel.png
    - results/figures/score_vs_identity.png
================================================================================
"""

from pathlib import Path
from typing import Optional, Tuple
import logging
import sys

# Pandas
try:
    import pandas as pd
except ImportError:
    print("❌ Error: Pandas is not installed.")
    print("   Run: pip install pandas")
    sys.exit(1)

# Matplotlib
try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.ticker import PercentFormatter, MaxNLocator
    from matplotlib.patches import FancyBboxPatch
except ImportError:
    print("❌ Error: Matplotlib is not installed.")
    print("   Run: pip install matplotlib")
    sys.exit(1)

# Seaborn
try:
    import seaborn as sns
except ImportError:
    print("❌ Error: Seaborn is not installed.")
    print("   Run: pip install seaborn")
    sys.exit(1)

# NumPy
try:
    import numpy as np
except ImportError:
    print("❌ Error: NumPy is not installed.")
    print("   Run: pip install numpy")
    sys.exit(1)

# Import publication style
try:
    from src.paper_style import (
        apply_publication_style, COLORS, PALETTE_CATEGORICAL,
        save_figure, add_panel_label, style_axis, annotate_stats,
        get_figure_size
    )
    HAS_PAPER_STYLE = True
except ImportError:
    HAS_PAPER_STYLE = False
    print("⚠️ paper_style module not found, using default style")


# ============================================================================
# CONFIGURATION
# ============================================================================

# Project paths
BASE_DIR: Path = Path(__file__).parent.resolve()
HOMOLOGY_CSV: Path = BASE_DIR / "results" / "logs" / "homology_scores.csv"
MOTIF_CSV: Path = BASE_DIR / "results" / "logs" / "motif_analysis.csv"
FIGURES_DIR: Path = BASE_DIR / "results" / "figures"

# Filtering parameters (must match 04_filter_homology.py)
MIN_IDENTITY: float = 20.0
MAX_IDENTITY: float = 90.0

# ============================================================================
# Q1 COLOR PALETTE (if paper_style is not available)
# ============================================================================
if not HAS_PAPER_STYLE:
    COLORS = {
        'primary': '#1A5276',
        'secondary': '#148F77',
        'accent': '#E74C3C',
        'data_blue': '#2E86AB',
        'data_green': '#28A745',
        'dark': '#2C3E50',
        'gray_dark': '#566573',
        'gray_light': '#BDC3C7',
        'background': '#FFFFFF',
        'grid': '#E5E8E8',
        'hit': '#E74C3C',
        'non_hit': '#95A5A6',
        'quality_high': '#27AE60',
    }
    PALETTE_CATEGORICAL = ['#2E86AB', '#E74C3C', '#28A745', '#F39C12', '#8E44AD']

# ============================================================================
# Q1 STYLE CONFIGURATION
# ============================================================================

def setup_q1_style():
    """Configure Q1 Scopus publication style."""
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
            'xtick.color': '#566573',
            'ytick.color': '#566573',
            'xtick.major.width': 1.0,
            'ytick.major.width': 1.0,
            'legend.fontsize': 9,
            'legend.framealpha': 0.95,
            'legend.edgecolor': '#BDC3C7',
            'figure.facecolor': 'white',
            'figure.dpi': 150,
            'savefig.dpi': 300,
            'savefig.bbox': 'tight',
            'savefig.facecolor': 'white',
            'grid.alpha': 0.5,
            'grid.linewidth': 0.5,
        })

# Logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# ============================================================================
# LOADING FUNCTIONS
# ============================================================================

def load_homology_data(filepath: Path) -> pd.DataFrame:
    """
    Load homology data from CSV.
    
    Args:
        filepath (Path): Path to the CSV file.
    
    Returns:
        pd.DataFrame: DataFrame with homology scores.
    """
    if not filepath.exists():
        logger.error(f"❌ File not found: {filepath}")
        logger.error("   Run first: python 04_filter_homology.py")
        raise FileNotFoundError(f"File not found: {filepath}")
    
    df = pd.read_csv(filepath)
    logger.info(f"✅ Homology data loaded: {len(df)} records")
    
    return df


def load_motif_data(filepath: Path) -> Optional[pd.DataFrame]:
    """
    Load motif data from CSV (optional).
    
    Args:
        filepath (Path): Path to the CSV file.
    
    Returns:
        Optional[pd.DataFrame]: DataFrame with motif analysis or None.
    """
    if not filepath.exists():
        logger.warning(f"⚠️ Motif file not found: {filepath}")
        return None
    
    df = pd.read_csv(filepath)
    logger.info(f"✅ Motif data loaded: {len(df)} records")
    
    return df


# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

def plot_identity_distribution(
    df: pd.DataFrame,
    output_path: Path
) -> None:
    """
    Generate a professional histogram of identity distribution.
    Q1 Scopus style with Nature/Cell design.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Crear figura con proporciones áureas
    fig, ax = plt.subplots(figsize=(8, 5.5))
    
    # Datos de identidad
    identity = df['Identity_Percent']
    passed = df[df['Passed_Filter'] == True]
    failed = df[df['Passed_Filter'] == False]
    
    # Histograma principal con estilo profesional
    n, bins, patches = ax.hist(
        identity,
        bins=40,
        color=COLORS['primary'],
        edgecolor='white',
        alpha=0.85,
        linewidth=0.8,
        zorder=3
    )
    
    # Colorear barras según región
    for i, patch in enumerate(patches):
        bin_center = (bins[i] + bins[i+1]) / 2
        if MIN_IDENTITY <= bin_center <= MAX_IDENTITY:
            patch.set_facecolor(COLORS['data_green'])
            patch.set_alpha(0.9)
        else:
            patch.set_facecolor(COLORS['non_hit'])
            patch.set_alpha(0.7)
    
    # Líneas de corte con estilo elegante
    for threshold, label in [(MIN_IDENTITY, f'Min: {MIN_IDENTITY}%'), 
                             (MAX_IDENTITY, f'Max: {MAX_IDENTITY}%')]:
        ax.axvline(
            threshold, 
            color=COLORS['accent'], 
            linestyle='--', 
            linewidth=2.0,
            zorder=4,
            alpha=0.9
        )
        # Anotación elegante
        y_pos = ax.get_ylim()[1] * 0.92
        ax.annotate(
            label,
            xy=(threshold, y_pos),
            xytext=(threshold + (3 if threshold == MIN_IDENTITY else -3), y_pos),
            fontsize=9,
            fontweight='bold',
            color=COLORS['accent'],
            ha='left' if threshold == MIN_IDENTITY else 'right',
            va='center'
        )
    
    # Sombrear región de selección
    ax.axvspan(
        MIN_IDENTITY, 
        MAX_IDENTITY, 
        alpha=0.08, 
        color=COLORS['data_green'],
        zorder=1
    )
    
    # Estadísticas en caja elegante
    stats_text = (
        f"Total sequences: {len(df):,}\n"
        f"Selected: {len(passed):,} ({len(passed)/len(df)*100:.1f}%)\n"
        f"Rejected: {len(failed):,}\n"
        f"Mean identity: {identity.mean():.1f}%"
    )
    
    ax.text(
        0.97, 0.95, 
        stats_text,
        transform=ax.transAxes,
        fontsize=9,
        fontfamily='monospace',
        verticalalignment='top',
        horizontalalignment='right',
        bbox=dict(
            boxstyle='round,pad=0.5',
            facecolor='white',
            edgecolor=COLORS['gray_light'],
            alpha=0.95,
            linewidth=0.8
        ),
        zorder=5
    )
    
    # Configurar ejes con estilo Q1
    ax.set_xlabel('Sequence Identity to IsPETase (%)', fontsize=11, fontweight='medium')
    ax.set_ylabel('Number of Sequences', fontsize=11, fontweight='medium')
    ax.set_title(
        'Distribution of Sequence Identity\n'
        'Deep-PETase-Mining Homology Filtering',
        fontsize=12, 
        fontweight='bold',
        pad=15,
        color=COLORS['dark']
    )
    
    # Límites y ticks
    ax.set_xlim(0, 100)
    ax.xaxis.set_major_locator(MaxNLocator(11))
    ax.yaxis.set_major_locator(MaxNLocator(8))
    
    # Spines elegantes
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.2)
    ax.spines['bottom'].set_linewidth(1.2)
    
    # Grid sutil
    ax.grid(True, alpha=0.4, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    
    # Leyenda profesional
    legend_elements = [
        mpatches.Patch(color=COLORS['data_green'], alpha=0.9, label='Selected region'),
        mpatches.Patch(color=COLORS['non_hit'], alpha=0.7, label='Rejected'),
        plt.Line2D([0], [0], color=COLORS['accent'], linestyle='--', 
                   linewidth=2, label='Identity thresholds')
    ]
    ax.legend(
        handles=legend_elements,
        loc='upper left',
        framealpha=0.95,
        edgecolor=COLORS['gray_light'],
        fontsize=9
    )
    
    # Save with high quality
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, facecolor='white', edgecolor='none',
                bbox_inches='tight', pad_inches=0.1)
    plt.close()
    
    logger.info(f"📊 Q1 plot saved: {output_path}")


def plot_filtering_funnel(
    df_homology: pd.DataFrame,
    df_motif: Optional[pd.DataFrame],
    output_path: Path
) -> None:
    """
    Generate a professional Q1-style funnel diagram.
    Elegant visualization of the progressive filtering process.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Calculate numbers for each stage
    total_initial = len(df_homology)
    passed_homology = len(df_homology[df_homology['Passed_Filter'] == True])
    
    if df_motif is not None:
        passed_motif = len(df_motif[df_motif['Passed_Filter'] == True])
    else:
        passed_motif = passed_homology
    
    # Data for the funnel
    stages = ['NCBI Mining\n(Raw Candidates)', 
              'Homology Filter\n(20-90% identity)', 
              'Motif Validation\n(G-x-S-x-G)']
    values = [total_initial, passed_homology, passed_motif]
    
    # Professional colors with gradient
    colors = [COLORS['primary'], COLORS['secondary'], COLORS['data_green']]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Draw centered horizontal bars (funnel)
    y_positions = [2, 1, 0]
    max_width = 0.9
    
    for i, (stage, value, y, color) in enumerate(zip(stages, values, y_positions, colors)):
        # Width proportional to value
        width = (value / total_initial) * max_width
        left = (1 - width) / 2
        
        # Bar with rounded edges (simulated)
        bar = ax.barh(
            y, 
            width, 
            height=0.55,
            left=left,
            color=color,
            edgecolor='white',
            linewidth=2,
            alpha=0.9,
            zorder=3
        )
        
        # Main number
        percentage = (value / total_initial) * 100
        ax.text(
            0.5, y + 0.05,
            f'{value:,}',
            ha='center', va='center',
            fontsize=16, fontweight='bold',
            color='white',
            zorder=4
        )
        
        # Percentage below
        if i > 0:
            ax.text(
                0.5, y - 0.18,
                f'({percentage:.1f}% retained)',
                ha='center', va='center',
                fontsize=9,
                color='white',
                alpha=0.9,
                zorder=4
            )
        
        # Stage label
        ax.text(
            0.02, y,
            stage,
            ha='left', va='center',
            fontsize=10,
            fontweight='medium',
            color=COLORS['dark'],
            zorder=4
        )
        
        # Reduction arrow
        if i < len(values) - 1:
            reduction = (1 - values[i+1]/values[i]) * 100
            ax.annotate(
                '',
                xy=(0.5, y - 0.35),
                xytext=(0.5, y - 0.65),
                arrowprops=dict(
                    arrowstyle='->',
                    color=COLORS['gray_dark'],
                    lw=2,
                    shrinkA=0,
                    shrinkB=0
                ),
                zorder=2
            )
            ax.text(
                0.58, y - 0.5,
                f'-{reduction:.0f}%',
                ha='left', va='center',
                fontsize=9,
                fontweight='bold',
                color=COLORS['accent'],
                zorder=4
            )
    
    # Configure axes
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.7, 2.7)
    ax.axis('off')
    
    # Professional title
    ax.set_title(
        'Filtering Funnel: Progressive Candidate Selection\n'
        'Deep-PETase-Mining Pipeline',
        fontsize=13,
        fontweight='bold',
        pad=20,
        color=COLORS['dark']
    )
    
    # Legend
    legend_elements = [
        mpatches.Patch(color=colors[0], alpha=0.9, label='Initial candidates'),
        mpatches.Patch(color=colors[1], alpha=0.9, label='Homology-filtered'),
        mpatches.Patch(color=colors[2], alpha=0.9, label='Motif-validated')
    ]
    ax.legend(
        handles=legend_elements,
        loc='upper right',
        framealpha=0.95,
        edgecolor=COLORS['gray_light'],
        fontsize=9
    )
    
    # Final summary box
    final_rate = (passed_motif / total_initial) * 100
    summary = f"Final Selection Rate: {final_rate:.1f}%"
    ax.text(
        0.5, -0.5,
        summary,
        ha='center', va='center',
        fontsize=11,
        fontweight='bold',
        color=COLORS['data_green'],
        bbox=dict(
            boxstyle='round,pad=0.4',
            facecolor='white',
            edgecolor=COLORS['data_green'],
            linewidth=1.5
        )
    )
    
    # Save
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, facecolor='white', edgecolor='none',
                bbox_inches='tight', pad_inches=0.15)
    plt.close()
    
    logger.info(f"📊 Q1 plot saved: {output_path}")


def plot_score_vs_identity(
    df: pd.DataFrame,
    output_path: Path
) -> None:
    """
    Generate a professional scatter plot of score vs identity.
    Q1 style with optimized colors, sizes and transparency.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Separate by filter
    passed = df[df['Passed_Filter'] == True]
    failed = df[df['Passed_Filter'] == False]
    
    # Scatter plot with professional style
    # First the rejected ones (background)
    ax.scatter(
        failed['Identity_Percent'], 
        failed['Score'],
        c=COLORS['non_hit'],
        alpha=0.35,
        s=35,
        label=f'Rejected (n={len(failed):,})',
        edgecolors='white',
        linewidth=0.3,
        zorder=2
    )
    
    # Then the selected ones (foreground)
    scatter = ax.scatter(
        passed['Identity_Percent'], 
        passed['Score'],
        c=passed['Identity_Percent'],
        cmap='viridis',
        alpha=0.75,
        s=55,
        label=f'Selected (n={len(passed):,})',
        edgecolors='white',
        linewidth=0.5,
        zorder=3
    )
    
    # Cutoff lines with elegant style
    ax.axvline(MIN_IDENTITY, color=COLORS['accent'], linestyle='--', 
               linewidth=1.8, alpha=0.8, zorder=4)
    ax.axvline(MAX_IDENTITY, color=COLORS['accent'], linestyle='--', 
               linewidth=1.8, alpha=0.8, zorder=4)
    
    # Shade selection region
    ax.axvspan(MIN_IDENTITY, MAX_IDENTITY, alpha=0.06, 
               color=COLORS['data_green'], zorder=1)
    
    # Colorbar for identity
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8, aspect=25, pad=0.02)
    cbar.set_label('Sequence Identity (%)', fontsize=10, fontweight='medium')
    cbar.ax.tick_params(labelsize=9)
    cbar.outline.set_linewidth(0.8)
    
    # Configure axes
    ax.set_xlabel('Sequence Identity (%)', fontsize=11, fontweight='medium')
    ax.set_ylabel('Alignment Score (bits)', fontsize=11, fontweight='medium')
    ax.set_title(
        'Alignment Quality vs Sequence Identity\n'
        'Candidate Selection Overview',
        fontsize=12,
        fontweight='bold',
        pad=15,
        color=COLORS['dark']
    )
    
    # Professional spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.2)
    ax.spines['bottom'].set_linewidth(1.2)
    
    # Legend
    ax.legend(
        loc='upper left',
        framealpha=0.95,
        edgecolor=COLORS['gray_light'],
        fontsize=9
    )
    
    # Subtle grid
    ax.grid(True, alpha=0.4, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    
    # Statistics
    stats_text = (
        f"Correlation: {passed['Identity_Percent'].corr(passed['Score']):.3f}\n"
        f"Mean score: {passed['Score'].mean():.1f}"
    )
    ax.text(
        0.97, 0.05,
        stats_text,
        transform=ax.transAxes,
        fontsize=9,
        fontfamily='monospace',
        va='bottom', ha='right',
        bbox=dict(
            boxstyle='round,pad=0.4',
            facecolor='white',
            edgecolor=COLORS['gray_light'],
            alpha=0.95
        )
    )
    
    # Save without tight_layout to avoid conflict with colorbar
    plt.savefig(output_path, dpi=300, facecolor='white', edgecolor='none',
                bbox_inches='tight', pad_inches=0.1)
    plt.close()
    
    logger.info(f"📊 Q1 plot saved: {output_path}")


def print_summary(
    df_homology: pd.DataFrame,
    df_motif: Optional[pd.DataFrame]
) -> None:
    """
    Print a summary of the visualized data.
    
    Args:
        df_homology (pd.DataFrame): Homology data.
        df_motif (Optional[pd.DataFrame]): Motif data.
    """
    print("\n" + "="*70)
    print("📊 VISUALIZATION SUMMARY")
    print("="*70)
    
    # Homology statistics
    identity = df_homology['Identity_Percent']
    passed = df_homology[df_homology['Passed_Filter'] == True]
    
    print(f"\n📈 Identity Distribution:")
    print(f"   • Minimum:  {identity.min():.1f}%")
    print(f"   • Maximum:  {identity.max():.1f}%")
    print(f"   • Mean:     {identity.mean():.1f}%")
    print(f"   • Median:   {identity.median():.1f}%")
    
    print(f"\n🎯 Homology Filtering Results:")
    print(f"   • Total analyzed:     {len(df_homology):,}")
    print(f"   • Within range:       {len(passed):,}")
    print(f"   • Selection rate:     {(len(passed)/len(df_homology))*100:.1f}%")
    
    if df_motif is not None:
        motif_passed = df_motif[df_motif['Passed_Filter'] == True]
        print(f"\n🧬 Motif Filtering Results:")
        print(f"   • Analyzed:           {len(df_motif):,}")
        print(f"   • With G-x-S-x-G:     {len(motif_passed):,}")
        print(f"   • Selection rate:     {(len(motif_passed)/len(df_motif))*100:.1f}%")
    
    print("\n" + "="*70)


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main() -> None:
    """
    Main function that generates all Q1 visualizations.
    """
    print("\n" + "="*70)
    print("   DEEP-PETASE-MINING - Phase 2: Q1 Publication Figures")
    print("="*70 + "\n")
    
    # Configure Q1 style
    setup_q1_style()
    logger.info("✓ Q1 Scopus style configured")
    
    # Step 1: Load data
    logger.info("Step 1/4: Loading data...")
    try:
        df_homology = load_homology_data(HOMOLOGY_CSV)
    except FileNotFoundError:
        sys.exit(1)
    
    df_motif = load_motif_data(MOTIF_CSV)
    
    # Step 2: Identity distribution plot
    logger.info("Step 2/4: Generating identity histogram (Q1)...")
    plot_identity_distribution(
        df_homology,
        FIGURES_DIR / "identity_distribution.png"
    )
    
    # Step 3: Funnel diagram
    logger.info("Step 3/4: Generating funnel diagram (Q1)...")
    plot_filtering_funnel(
        df_homology,
        df_motif,
        FIGURES_DIR / "filtering_funnel.png"
    )
    
    # Step 4: Score vs identity scatter plot
    logger.info("Step 4/4: Generating scatter plot (Q1)...")
    plot_score_vs_identity(
        df_homology,
        FIGURES_DIR / "score_vs_identity.png"
    )
    
    # Summary
    print_summary(df_homology, df_motif)
    
    # Final message
    print(f"\n📁 Q1 plots generated in: {FIGURES_DIR}")
    print(f"   • identity_distribution.png (300 DPI)")
    print(f"   • filtering_funnel.png (300 DPI)")
    print(f"   • score_vs_identity.png (300 DPI)")
    
    print("\n✅ Phase 2 Q1 figures completed!")
    print("\n" + "="*70 + "\n")


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
