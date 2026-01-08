#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
02_exploratory_analysis.py - Exploratory Data Analysis (EDA)
================================================================================

Description:
    This script performs an exploratory analysis of the sequences downloaded
    from NCBI. The goal is to validate data quality and generate
    informative visualizations for the portfolio.

    Analyses performed:
    1. Descriptive statistics (mean, median, std. deviation)
    2. Protein length distribution
    3. Organism frequency analysis
    4. Duplicate detection
    5. Report chart generation

Date: December 2025
Project: Deep-PETase-Mining - Phase 1

Usage:
    python 02_exploratory_analysis.py

Dependencies:
    - Python 3.9+
    - Pandas
    - Matplotlib
    - Seaborn

Outputs:
    - results/figures/distribucion_longitud.png
    - results/figures/top_organismos.png
    - Statistics printed to console and saved in log
================================================================================
"""

# ============================================================================
# IMPORTS
# ============================================================================
from pathlib import Path
from typing import Tuple, Optional
from datetime import datetime
import logging
import sys

# Pandas for data analysis
try:
    import pandas as pd
except ImportError:
    print("❌ Error: Pandas is not installed.")
    print("   Run: pip install pandas")
    sys.exit(1)

# Matplotlib for visualization
try:
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
except ImportError:
    print("❌ Error: Matplotlib is not installed.")
    print("   Run: pip install matplotlib")
    sys.exit(1)

# Seaborn for elegant statistical charts
try:
    import seaborn as sns
except ImportError:
    print("❌ Error: Seaborn is not installed.")
    print("   Run: pip install seaborn")
    sys.exit(1)

# Numpy for statistical calculations
try:
    import numpy as np
except ImportError:
    print("❌ Error: NumPy is not installed.")
    print("   Run: pip install numpy")
    sys.exit(1)


# ============================================================================
# GLOBAL CONFIGURATION
# ============================================================================

# Project paths
BASE_DIR: Path = Path(__file__).parent.resolve()
METADATA_DIR: Path = BASE_DIR / "data" / "metadata"
FIGURES_DIR: Path = BASE_DIR / "results" / "figures"
LOGS_DIR: Path = BASE_DIR / "results" / "logs"

# Input file
INPUT_CSV: Path = METADATA_DIR / "candidates_info.csv"

# Output files
OUTPUT_HIST: Path = FIGURES_DIR / "distribucion_longitud.png"
OUTPUT_BAR: Path = FIGURES_DIR / "top_organismos.png"

# Chart style configuration
# We use a professional and clean style
plt.style.use('seaborn-v0_8-whitegrid')  # Clean style with grid
sns.set_palette("husl")  # Vibrant but professional color palette

# Font configuration for charts
plt.rcParams.update({
    'font.size': 11,
    'font.family': 'sans-serif',
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 16,
    'figure.dpi': 150,  # High resolution for portfolio
    'savefig.dpi': 300,  # Even higher resolution for saving
    'savefig.bbox': 'tight'
})


# ============================================================================
# LOGGING CONFIGURATION
# ============================================================================
def setup_logging() -> logging.Logger:
    """
    Configures the logging system with console and file output.
    
    Returns:
        logging.Logger: Configured logger.
    """
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    log_filename = LOGS_DIR / f"eda_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_filename, encoding='utf-8'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    logger = logging.getLogger(__name__)
    logger.info(f"📝 Log saved to: {log_filename}")
    
    return logger


# ============================================================================
# DATA LOADING FUNCTIONS
# ============================================================================
def load_data(filepath: Path) -> pd.DataFrame:
    """
    Loads the metadata CSV file into a DataFrame.
    
    Args:
        filepath (Path): Path to the CSV file.
    
    Returns:
        pd.DataFrame: DataFrame with the metadata.
    
    Raises:
        FileNotFoundError: If the file does not exist.
    
    Note for students:
        It is always good practice to verify that the file exists
        before attempting to load it, and to provide clear error messages.
    """
    if not filepath.exists():
        logger.error(f"❌ File not found: {filepath}")
        logger.error("   Make sure to run first: python 01_mining_ncbi.py")
        raise FileNotFoundError(f"File not found: {filepath}")
    
    try:
        df = pd.read_csv(filepath)
        logger.info(f"✅ Data loaded: {len(df)} records from {filepath.name}")
        return df
    except Exception as e:
        logger.error(f"❌ Error reading CSV: {e}")
        raise


# ============================================================================
# ANALYSIS FUNCTIONS
# ============================================================================
def calculate_statistics(df: pd.DataFrame) -> dict:
    """
    Calculates basic descriptive statistics for the dataset.
    
    Args:
        df (pd.DataFrame): DataFrame with the data.
    
    Returns:
        dict: Dictionary with all calculated statistics.
    
    Note for students:
        Descriptive statistics are the first step in any EDA.
        They give us an overview of the data before visualizing it.
    """
    stats = {
        "total_sequences": len(df),
        "mean_length": df["length"].mean(),
        "median_length": df["length"].median(),
        "std_length": df["length"].std(),
        "min_length": df["length"].min(),
        "max_length": df["length"].max(),
        "q25_length": df["length"].quantile(0.25),
        "q75_length": df["length"].quantile(0.75),
        "unique_organisms": df["organism"].nunique(),
        "most_common_organism": df["organism"].mode().iloc[0] if not df["organism"].mode().empty else "N/A",
        "most_common_count": df["organism"].value_counts().iloc[0] if len(df) > 0 else 0
    }
    
    return stats


def check_duplicates(df: pd.DataFrame) -> Tuple[int, pd.DataFrame]:
    """
    Checks if duplicates exist in the dataset based on ID.
    
    Args:
        df (pd.DataFrame): DataFrame with the data.
    
    Returns:
        Tuple[int, pd.DataFrame]: 
            - Number of duplicates found
            - DataFrame with duplicate records (if any)
    
    Note for students:
        Duplicates can bias statistical analyses.
        It is always important to detect them before continuing.
    """
    # Check duplicates by ID
    duplicates = df[df.duplicated(subset=["id"], keep=False)]
    num_duplicates = len(duplicates) // 2  # Divide by 2 because both copies are counted
    
    return num_duplicates, duplicates


def print_statistics(stats: dict, duplicates: int) -> None:
    """
    Prints statistics in a formatted and visually appealing way.
    
    Args:
        stats (dict): Dictionary with the statistics.
        duplicates (int): Number of duplicates found.
    """
    print("\n" + "="*70)
    print("📊 DATASET DESCRIPTIVE STATISTICS")
    print("="*70)
    
    print("\n📈 GENERAL SUMMARY:")
    print(f"   • Total sequences: {stats['total_sequences']:,}")
    print(f"   • Unique organisms: {stats['unique_organisms']:,}")
    
    print("\n📏 LENGTH DISTRIBUTION (amino acids):")
    print(f"   • Mean:               {stats['mean_length']:.1f} aa")
    print(f"   • Median:             {stats['median_length']:.1f} aa")
    print(f"   • Standard deviation: {stats['std_length']:.1f} aa")
    print(f"   • Minimum:            {stats['min_length']} aa")
    print(f"   • Maximum:            {stats['max_length']} aa")
    print(f"   • 25th Percentile:    {stats['q25_length']:.1f} aa")
    print(f"   • 75th Percentile:    {stats['q75_length']:.1f} aa")
    
    print("\n🦠 MOST FREQUENT ORGANISM:")
    print(f"   • {stats['most_common_organism']}")
    print(f"   • Occurrences: {stats['most_common_count']:,} sequences")
    
    print("\n🔍 QUALITY CONTROL:")
    if duplicates > 0:
        print(f"   ⚠️ Found {duplicates} pairs of duplicates by ID")
        print(f"      Consider removing them before further analysis")
    else:
        print(f"   ✅ No duplicates found by ID")
    
    print("\n" + "="*70)


# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================
def plot_length_distribution(df: pd.DataFrame, output_path: Path) -> None:
    """
    Generates a histogram of protein length distribution.
    
    Args:
        df (pd.DataFrame): DataFrame with the data.
        output_path (Path): Path where to save the figure.
    
    Note for students:
        A histogram is essential to understand the distribution of your data.
        For proteins, we expect an approximately normal distribution,
        but skewed towards smaller proteins (there are more short genes).
    """
    # Ensure the directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Create figure with optimized size
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Calculate optimal number of bins using Freedman-Diaconis rule
    q75, q25 = np.percentile(df["length"], [75, 25])
    iqr = q75 - q25
    bin_width = 2 * iqr / (len(df) ** (1/3))
    n_bins = int((df["length"].max() - df["length"].min()) / bin_width) if bin_width > 0 else 30
    n_bins = max(20, min(n_bins, 50))  # Between 20 and 50 bins
    
    # Histogram with KDE (Kernel Density Estimation)
    sns.histplot(
        data=df,
        x="length",
        bins=n_bins,
        kde=True,
        color="#3498db",
        edgecolor="white",
        linewidth=0.5,
        alpha=0.7,
        ax=ax
    )
    
    # Vertical lines for key statistics
    mean_val = df["length"].mean()
    median_val = df["length"].median()
    
    ax.axvline(mean_val, color="#e74c3c", linestyle="--", linewidth=2, 
               label=f"Mean: {mean_val:.0f} aa")
    ax.axvline(median_val, color="#2ecc71", linestyle="-.", linewidth=2, 
               label=f"Median: {median_val:.0f} aa")
    
    # Shade interquartile range (IQR) region
    ax.axvspan(q25, q75, alpha=0.15, color="#9b59b6", 
               label=f"IQR: {q25:.0f}-{q75:.0f} aa")
    
    # Configure axes and titles
    ax.set_xlabel("Protein length (amino acids)", fontsize=12, fontweight='bold')
    ax.set_ylabel("Frequency (number of sequences)", fontsize=12, fontweight='bold')
    ax.set_title("Length Distribution of Candidate Hypothetical Proteins\n"
                 "Deep-PETase-Mining - Phase 1", 
                 fontsize=14, fontweight='bold', pad=20)
    
    # Legend
    ax.legend(loc="upper right", framealpha=0.95, edgecolor="gray")
    
    # Add annotation with statistics
    stats_text = (f"n = {len(df):,}\n"
                  f"σ = {df['length'].std():.1f} aa\n"
                  f"Range: {df['length'].min()}-{df['length'].max()} aa")
    
    ax.text(0.98, 0.75, stats_text,
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment='top',
            horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # X axis format
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'{int(x):,}'))
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, facecolor='white', edgecolor='none')
    plt.close()
    
    logger.info(f"📊 Histogram saved to: {output_path}")


def plot_top_organisms(df: pd.DataFrame, output_path: Path, top_n: int = 10) -> None:
    """
    Generates a horizontal bar chart with the most frequent organisms.
    
    Args:
        df (pd.DataFrame): DataFrame with the data.
        output_path (Path): Path where to save the figure.
        top_n (int): Number of top organisms to show.
    
    Note for students:
        Horizontal bars are better than vertical ones when
        category names are long (like organism names).
        This improves chart readability.
    """
    # Ensure the directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Calculate frequencies
    organism_counts = df["organism"].value_counts().head(top_n)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Gradient color palette
    colors = sns.color_palette("viridis", n_colors=top_n)[::-1]
    
    # Horizontal bar chart
    bars = ax.barh(
        y=range(len(organism_counts)),
        width=organism_counts.values,
        color=colors,
        edgecolor="white",
        linewidth=0.5
    )
    
    # Labels on bars (number of sequences)
    for i, (bar, count) in enumerate(zip(bars, organism_counts.values)):
        # Calculate percentage
        percentage = (count / len(df)) * 100
        
        # Text position (inside or outside bar depending on size)
        if bar.get_width() > organism_counts.max() * 0.3:
            ax.text(bar.get_width() - organism_counts.max() * 0.02, bar.get_y() + bar.get_height()/2,
                   f'{count:,} ({percentage:.1f}%)',
                   ha='right', va='center', color='white', fontweight='bold', fontsize=10)
        else:
            ax.text(bar.get_width() + organism_counts.max() * 0.02, bar.get_y() + bar.get_height()/2,
                   f'{count:,} ({percentage:.1f}%)',
                   ha='left', va='center', color='black', fontsize=10)
    
    # Configure axes
    ax.set_yticks(range(len(organism_counts)))
    
    # Truncate long organism names for better visualization
    labels = [name[:50] + "..." if len(name) > 50 else name 
              for name in organism_counts.index]
    ax.set_yticklabels(labels, fontsize=10)
    
    # Invert Y axis so the most frequent is at the top
    ax.invert_yaxis()
    
    # Titles and labels
    ax.set_xlabel("Number of sequences", fontsize=12, fontweight='bold')
    ax.set_ylabel("Organism / Metagenomic source", fontsize=12, fontweight='bold')
    ax.set_title(f"Top {top_n} Most Frequent Organisms/Sources\n"
                 "Hypothetical Proteins from Plastic-Related Metagenomes",
                 fontsize=14, fontweight='bold', pad=20)
    
    # Add subtle vertical grid
    ax.xaxis.grid(True, linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)
    
    # X axis format
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'{int(x):,}'))
    
    # Additional information
    total_text = f"Total sequences: {len(df):,}\nUnique organisms: {df['organism'].nunique():,}"
    ax.text(0.98, 0.02, total_text,
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment='bottom',
            horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, facecolor='white', edgecolor='none')
    plt.close()
    
    logger.info(f"📊 Bar chart saved to: {output_path}")


def generate_summary_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Generates a summary table grouped by organism.
    
    Args:
        df (pd.DataFrame): DataFrame with the data.
    
    Returns:
        pd.DataFrame: Summary table with statistics per organism.
    """
    summary = df.groupby("organism").agg(
        count=("id", "count"),
        mean_length=("length", "mean"),
        min_length=("length", "min"),
        max_length=("length", "max")
    ).round(1).sort_values("count", ascending=False)
    
    return summary


# ============================================================================
# MAIN FUNCTION
# ============================================================================
def main() -> None:
    """
    Main function that executes the complete exploratory analysis.
    
    Steps:
    1. Load data
    2. Calculate statistics
    3. Detect duplicates
    4. Generate visualizations
    5. Show summary
    """
    print("\n" + "🧬"*30)
    print("   DEEP-PETASE-MINING - Exploratory Analysis (EDA)")
    print("🧬"*30 + "\n")
    
    # Step 1: Load data
    logger.info("📂 Step 1/4: Loading data...")
    try:
        df = load_data(INPUT_CSV)
    except FileNotFoundError:
        sys.exit(1)
    
    # Show first rows
    print("\n📋 Data preview:")
    print("-" * 70)
    print(df.head(10).to_string())
    print("-" * 70)
    
    # Step 2: Calculate statistics
    logger.info("📊 Step 2/4: Calculating statistics...")
    stats = calculate_statistics(df)
    num_duplicates, duplicates_df = check_duplicates(df)
    
    # Show statistics
    print_statistics(stats, num_duplicates)
    
    # If there are duplicates, show which ones
    if num_duplicates > 0:
        logger.warning(f"⚠️ Duplicates found:")
        print("\n🔍 Duplicate IDs:")
        print(duplicates_df[["id", "description", "organism"]].to_string())
    
    # Step 3: Generate visualizations
    logger.info("🎨 Step 3/4: Generating visualizations...")
    
    # Length histogram
    plot_length_distribution(df, OUTPUT_HIST)
    
    # Organism bar chart
    plot_top_organisms(df, OUTPUT_BAR, top_n=10)
    
    # Step 4: Final summary
    logger.info("📝 Step 4/4: Generating final summary...")
    
    print("\n" + "="*70)
    print("✅ EXPLORATORY ANALYSIS COMPLETED")
    print("="*70)
    
    print("\n📁 Generated files:")
    print(f"   • {OUTPUT_HIST}")
    print(f"   • {OUTPUT_BAR}")
    
    print("\n📋 PHASE 1 CHECKLIST:")
    print(f"   [{'✓' if (BASE_DIR / 'data' / 'raw').exists() else '✗'}] Folders created correctly")
    print(f"   [{'✓' if (BASE_DIR / 'data' / 'raw' / 'raw_candidates.fasta').exists() else '✗'}] raw_candidates.fasta file")
    print(f"   [{'✓' if INPUT_CSV.exists() else '✗'}] candidates_info.csv file")
    print(f"   [{'✓' if OUTPUT_HIST.exists() else '✗'}] distribucion_longitud.png chart")
    print(f"   [{'✓' if OUTPUT_BAR.exists() else '✗'}] top_organismos.png chart")
    
    # Check if the data makes biological sense
    print("\n🔬 BIOLOGICAL VALIDATION:")
    
    if 100 <= stats['mean_length'] <= 800:
        print(f"   ✅ Average length ({stats['mean_length']:.0f} aa) is in expected range")
    else:
        print(f"   ⚠️ Average length ({stats['mean_length']:.0f} aa) outside typical range")
    
    if stats['total_sequences'] >= 50:
        print(f"   ✅ Enough sequences ({stats['total_sequences']}) for analysis")
    else:
        print(f"   ⚠️ Few sequences ({stats['total_sequences']}), consider expanding search")
    
    if num_duplicates == 0:
        print(f"   ✅ No duplicates - clean data")
    else:
        print(f"   ⚠️ {num_duplicates} duplicates detected - review before Phase 2")
    
    print("\n" + "="*70)
    print("\n🎉 Phase 1 completed! You are ready for Phase 2: Sequence Analysis")
    print("\n💡 Suggestions for Phase 2:")
    print("   • Perform multiple alignment with Clustal Omega")
    print("   • Search for conserved domains with InterProScan")
    print("   • Compare with known PETases using BLAST")
    print("   • Structure prediction with AlphaFold/ESMFold")
    
    print("\n" + "🧬"*30 + "\n")


# ============================================================================
# ENTRY POINT
# ============================================================================
if __name__ == "__main__":
    # Configure logging
    logger = setup_logging()
    
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n⚠️ Execution cancelled by user.")
        logger.info("Execution cancelled by user")
    except Exception as e:
        logger.error(f"❌ Fatal error: {e}")
        raise
