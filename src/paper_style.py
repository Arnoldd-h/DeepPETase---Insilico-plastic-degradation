#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from pathlib import Path
from typing import Optional, Tuple, List

# ============================================================================
# PROFESSIONAL COLOR PALETTE (Colorblind-friendly)
# ============================================================================

# Main palette - Inspired by Nature/Science
COLORS = {
    # Primary colors
    'primary': '#1A5276',        # Deep blue (professional)
    'secondary': '#148F77',      # Teal green
    'accent': '#E74C3C',         # Coral red (for emphasis)
    'highlight': '#F39C12',      # Golden orange
    
    # Data colors
    'data_blue': '#2E86AB',      # Data blue
    'data_green': '#28A745',     # Success green
    'data_orange': '#E67E22',    # Warning orange  
    'data_red': '#C0392B',       # Strong red
    'data_purple': '#8E44AD',    # Purple
    'data_teal': '#17A2B8',      # Teal
    
    # Professional grays
    'dark': '#2C3E50',           # Near black
    'gray_dark': '#566573',      # Dark gray
    'gray_medium': '#7F8C8D',    # Medium gray
    'gray_light': '#BDC3C7',     # Light gray
    'gray_pale': '#ECF0F1',      # Very light gray
    
    # Backgrounds
    'background': '#FFFFFF',     # Pure white
    'background_alt': '#FAFAFA', # Off-white
    'grid': '#E5E8E8',           # Grid gray
    
    # Quality categories (pLDDT, etc.)
    'quality_high': '#27AE60',   # Green - High quality
    'quality_medium': '#F1C40F', # Yellow - Medium
    'quality_low': '#E67E22',    # Orange - Low
    'quality_very_low': '#E74C3C', # Red - Very low
    
    # Hits vs Non-hits
    'hit': '#E74C3C',            # Red for positive hits
    'non_hit': '#95A5A6',        # Gray for non-hits
    'best_hit': '#C0392B',       # Dark red for best hit
}

# Sequential palette for heatmaps
CMAP_SEQUENTIAL = LinearSegmentedColormap.from_list(
    'paper_sequential',
    ['#EBF5FB', '#3498DB', '#1A5276']
)

# Divergent palette 
CMAP_DIVERGENT = LinearSegmentedColormap.from_list(
    'paper_divergent',
    ['#E74C3C', '#FDFEFE', '#27AE60']
)

# Categorical palette (maximum contrast, colorblind-safe)
PALETTE_CATEGORICAL = [
    '#2E86AB',  # Blue
    '#E74C3C',  # Red
    '#28A745',  # Green
    '#F39C12',  # Orange
    '#8E44AD',  # Purple
    '#17A2B8',  # Teal
    '#E67E22',  # Dark orange
    '#1ABC9C',  # Turquoise
]

# Palette for quality gradients
PALETTE_QUALITY = ['#C0392B', '#E74C3C', '#F39C12', '#F1C40F', '#27AE60', '#1E8449']


# ============================================================================
# TYPOGRAPHY CONFIGURATION
# ============================================================================

FONT_CONFIG = {
    'family': 'sans-serif',
    'sans_serif': ['Arial', 'Helvetica Neue', 'Helvetica', 'DejaVu Sans', 'Liberation Sans'],
    'size_title': 14,
    'size_subtitle': 12,
    'size_label': 11,
    'size_tick': 10,
    'size_legend': 9,
    'size_annotation': 9,
    'weight_title': 'bold',
    'weight_label': 'medium',
}


# ============================================================================
# FIGURE CONFIGURATION
# ============================================================================

FIGURE_CONFIG = {
    # Standard sizes (inches) - journal compatible
    'single_column': (3.5, 2.8),      # Single column (89mm)
    'one_half_column': (5.5, 4.0),    # 1.5 columns
    'double_column': (7.2, 5.0),      # Double column (183mm)
    'full_page': (7.2, 9.0),          # Full page
    
    # For presentations
    'presentation': (12, 7),
    'poster': (16, 10),
    
    # DPI
    'dpi_screen': 150,
    'dpi_print': 300,
    'dpi_poster': 600,
    
    # Margins and spacing
    'padding': 0.15,
    'subplot_spacing': 0.25,
}


# ============================================================================
# MAIN STYLE FUNCTION
# ============================================================================

def apply_publication_style(
    style: str = 'nature',
    font_scale: float = 1.0,
    context: str = 'paper'
) -> None:
    """
    Apply professional style configuration for Q1 publications.
    
    Args:
        style: 'nature', 'science', 'cell', 'minimal'
        font_scale: Scale factor for fonts (1.0 = default)
        context: 'paper', 'poster', 'presentation'
    
    Example:
        >>> apply_publication_style(style='nature')
    """
    
    # Reset first
    plt.rcdefaults()
    
    # Base configuration according to context
    if context == 'poster':
        base_font = 14
        base_line = 2.0
    elif context == 'presentation':
        base_font = 12
        base_line = 1.8
    else:  # paper
        base_font = 10
        base_line = 1.2
    
    # Apply scale
    base_font = int(base_font * font_scale)
    
    # Complete matplotlib configuration
    publication_params = {
        # === FONTS ===
        'font.family': 'sans-serif',
        'font.sans-serif': FONT_CONFIG['sans_serif'],
        'font.size': base_font,
        'font.weight': 'normal',
        
        # === TEXT SIZES ===
        'axes.titlesize': base_font + 2,
        'axes.labelsize': base_font,
        'xtick.labelsize': base_font - 1,
        'ytick.labelsize': base_font - 1,
        'legend.fontsize': base_font - 1,
        'legend.title_fontsize': base_font,
        'figure.titlesize': base_font + 4,
        
        # === FONT WEIGHTS ===
        'axes.titleweight': 'bold',
        'axes.labelweight': 'medium',
        'figure.titleweight': 'bold',
        
        # === COLORS ===
        'text.color': COLORS['dark'],
        'axes.labelcolor': COLORS['dark'],
        'axes.edgecolor': COLORS['gray_dark'],
        'xtick.color': COLORS['gray_dark'],
        'ytick.color': COLORS['gray_dark'],
        
        # === AXES ===
        'axes.linewidth': base_line,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.facecolor': COLORS['background'],
        'axes.axisbelow': True,
        
        # === TICKS ===
        'xtick.major.size': 5,
        'xtick.major.width': base_line * 0.8,
        'xtick.minor.size': 3,
        'xtick.minor.width': base_line * 0.5,
        'xtick.direction': 'out',
        'xtick.major.pad': 5,
        
        'ytick.major.size': 5,
        'ytick.major.width': base_line * 0.8,
        'ytick.minor.size': 3,
        'ytick.minor.width': base_line * 0.5,
        'ytick.direction': 'out',
        'ytick.major.pad': 5,
        
        # === GRID ===
        'axes.grid': True,
        'grid.color': COLORS['grid'],
        'grid.linewidth': 0.5,
        'grid.alpha': 0.7,
        'grid.linestyle': '-',
        
        # === LEGEND ===
        'legend.frameon': True,
        'legend.framealpha': 0.95,
        'legend.edgecolor': COLORS['gray_light'],
        'legend.fancybox': True,
        'legend.shadow': False,
        'legend.borderpad': 0.5,
        'legend.labelspacing': 0.4,
        'legend.handlelength': 1.5,
        'legend.handleheight': 0.7,
        'legend.handletextpad': 0.5,
        'legend.columnspacing': 1.0,
        'legend.markerscale': 1.0,
        
        # === FIGURE ===
        'figure.facecolor': COLORS['background'],
        'figure.edgecolor': 'none',
        'figure.dpi': FIGURE_CONFIG['dpi_screen'],
        'figure.constrained_layout.use': True,
        'figure.autolayout': False,
        
        # === SAVING ===
        'savefig.dpi': FIGURE_CONFIG['dpi_print'],
        'savefig.facecolor': COLORS['background'],
        'savefig.edgecolor': 'none',
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.1,
        'savefig.transparent': False,
        
        # === HISTOGRAMS ===
        'hist.bins': 'auto',
        
        # === SCATTER ===
        'scatter.edgecolors': 'white',
        
        # === LINES ===
        'lines.linewidth': 1.5,
        'lines.markersize': 6,
        'lines.markeredgewidth': 0.8,
        'lines.markeredgecolor': 'white',
        
        # === PATCHES (bars, etc.) ===
        'patch.linewidth': 0.8,
        'patch.edgecolor': 'white',
        'patch.force_edgecolor': True,
        
        # === IMAGES ===
        'image.cmap': 'viridis',
        'image.interpolation': 'nearest',
        
        # === MATH ===
        'mathtext.fontset': 'dejavusans',
    }
    
    # Apply journal-specific styles
    if style == 'nature':
        # Nature uses more whitespace
        publication_params.update({
            'axes.spines.top': False,
            'axes.spines.right': False,
        })
    elif style == 'science':
        # Science is more compact
        publication_params['figure.constrained_layout.h_pad'] = 0.02
        publication_params['figure.constrained_layout.w_pad'] = 0.02
    elif style == 'cell':
        # Cell uses full borders
        publication_params.update({
            'axes.spines.top': True,
            'axes.spines.right': True,
        })
    elif style == 'minimal':
        # Minimalist
        publication_params.update({
            'axes.grid': False,
            'axes.spines.top': False,
            'axes.spines.right': False,
            'axes.spines.left': True,
            'axes.spines.bottom': True,
        })
    
    # Apply configuration
    plt.rcParams.update(publication_params)
    
    # Register custom colormaps
    try:
        mpl.colormaps.register(CMAP_SEQUENTIAL, name='paper_seq', force=True)
        mpl.colormaps.register(CMAP_DIVERGENT, name='paper_div', force=True)
    except Exception:
        pass  # Already registered


def get_figure_size(
    layout: str = 'single_column',
    aspect_ratio: Optional[float] = None
) -> Tuple[float, float]:
    """
    Get figure dimensions according to journal layout.
    
    Args:
        layout: 'single_column', 'one_half_column', 'double_column', 'full_page'
        aspect_ratio: Custom ratio (width/height)
    
    Returns:
        Tuple (width, height) in inches
    """
    base_size = FIGURE_CONFIG.get(layout, FIGURE_CONFIG['double_column'])
    
    if aspect_ratio:
        return (base_size[0], base_size[0] / aspect_ratio)
    
    return base_size


def save_figure(
    fig: plt.Figure,
    filename: str,
    output_dir: Path,
    formats: List[str] = ['png', 'pdf'],
    dpi: Optional[int] = None
) -> List[Path]:
    """
    Save figure in multiple formats with publication quality.
    
    Args:
        fig: Matplotlib figure
        filename: Base filename (without extension)
        output_dir: Output directory
        formats: List of formats ('png', 'pdf', 'svg', 'eps')
        dpi: Custom DPI (None = use default)
    
    Returns:
        List of saved file paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    saved_files = []
    
    for fmt in formats:
        output_path = output_dir / f"{filename}.{fmt}"
        
        # DPI according to format
        if dpi is None:
            if fmt == 'png':
                save_dpi = FIGURE_CONFIG['dpi_print']
            elif fmt in ['pdf', 'svg', 'eps']:
                save_dpi = None  # Vector format
            else:
                save_dpi = FIGURE_CONFIG['dpi_print']
        else:
            save_dpi = dpi
        
        # Save
        fig.savefig(
            output_path,
            format=fmt,
            dpi=save_dpi,
            bbox_inches='tight',
            facecolor=COLORS['background'],
            edgecolor='none',
            pad_inches=0.1
        )
        
        saved_files.append(output_path)
    
    return saved_files


def add_panel_label(
    ax: plt.Axes,
    label: str,
    loc: str = 'upper left',
    fontsize: int = 14,
    fontweight: str = 'bold'
) -> None:
    """
    Add panel label (A, B, C, etc.) in journal style.
    
    Args:
        ax: Matplotlib axes
        label: Label text ('A', 'B', etc.)
        loc: Location ('upper left', 'upper right', etc.)
        fontsize: Font size
        fontweight: Font weight
    """
    # Predefined positions
    positions = {
        'upper left': (-0.12, 1.08),
        'upper right': (1.02, 1.08),
        'lower left': (-0.12, -0.08),
        'lower right': (1.02, -0.08),
    }
    
    x, y = positions.get(loc, positions['upper left'])
    
    ax.text(
        x, y, label,
        transform=ax.transAxes,
        fontsize=fontsize,
        fontweight=fontweight,
        color=COLORS['dark'],
        va='top',
        ha='left'
    )


def style_axis(
    ax: plt.Axes,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    title: Optional[str] = None,
    grid: bool = True,
    spines: List[str] = ['bottom', 'left']
) -> None:
    """
    Apply consistent style to an axis.
    
    Args:
        ax: Matplotlib axes
        xlabel: X-axis label
        ylabel: Y-axis label  
        title: Subplot title
        grid: Show grid
        spines: List of spines to show
    """
    # Configure spines
    for spine in ['top', 'right', 'bottom', 'left']:
        ax.spines[spine].set_visible(spine in spines)
        if spine in spines:
            ax.spines[spine].set_linewidth(1.2)
            ax.spines[spine].set_color(COLORS['gray_dark'])
    
    # Labels
    if xlabel:
        ax.set_xlabel(xlabel, fontweight='medium', color=COLORS['dark'])
    if ylabel:
        ax.set_ylabel(ylabel, fontweight='medium', color=COLORS['dark'])
    if title:
        ax.set_title(title, fontweight='bold', color=COLORS['dark'], pad=10)
    
    # Grid
    if grid:
        ax.grid(True, alpha=0.5, linewidth=0.5, color=COLORS['grid'])
        ax.set_axisbelow(True)
    else:
        ax.grid(False)


def create_colorbar(
    mappable,
    ax: plt.Axes,
    label: str = '',
    orientation: str = 'vertical',
    shrink: float = 0.8
) -> plt.colorbar:
    """
    Create colorbar with professional style.
    
    Args:
        mappable: Mappable object (scatter, imshow, etc.)
        ax: Associated axes
        label: Colorbar label
        orientation: 'vertical' or 'horizontal'
        shrink: Shrink factor
    
    Returns:
        Colorbar object
    """
    cbar = plt.colorbar(
        mappable, 
        ax=ax,
        orientation=orientation,
        shrink=shrink,
        aspect=20,
        pad=0.02
    )
    
    cbar.set_label(label, fontweight='medium', color=COLORS['dark'])
    cbar.ax.tick_params(labelsize=FONT_CONFIG['size_tick'] - 1)
    cbar.outline.set_linewidth(0.8)
    cbar.outline.set_edgecolor(COLORS['gray_light'])
    
    return cbar


def annotate_stats(
    ax: plt.Axes,
    stats_dict: dict,
    loc: str = 'upper right',
    fontsize: int = 9
) -> None:
    """
    Add statistics box to the plot.
    
    Args:
        ax: Matplotlib axes
        stats_dict: Dictionary with statistics {'label': value}
        loc: Text location
        fontsize: Font size
    """
    # Build text
    lines = [f"{k}: {v}" for k, v in stats_dict.items()]
    text = '\n'.join(lines)
    
    # Positions
    positions = {
        'upper right': (0.97, 0.97),
        'upper left': (0.03, 0.97),
        'lower right': (0.97, 0.03),
        'lower left': (0.03, 0.03),
    }
    
    x, y = positions.get(loc, positions['upper right'])
    ha = 'right' if 'right' in loc else 'left'
    va = 'top' if 'upper' in loc else 'bottom'
    
    ax.text(
        x, y, text,
        transform=ax.transAxes,
        fontsize=fontsize,
        fontfamily='monospace',
        verticalalignment=va,
        horizontalalignment=ha,
        bbox=dict(
            boxstyle='round,pad=0.5',
            facecolor=COLORS['background_alt'],
            edgecolor=COLORS['gray_light'],
            alpha=0.95,
            linewidth=0.8
        )
    )


# ============================================================================
# UTILITY FUNCTIONS FOR SPECIFIC PLOTS
# ============================================================================

def create_quality_cmap():
    """Create colormap for quality data (pLDDT, etc.)."""
    return LinearSegmentedColormap.from_list(
        'quality',
        [COLORS['quality_very_low'], COLORS['quality_low'], 
         COLORS['quality_medium'], COLORS['quality_high']],
        N=256
    )


def get_category_colors(n_categories: int) -> List[str]:
    """
    Get list of colors for n categories.
    
    Args:
        n_categories: Number of categories
    
    Returns:
        List of hexadecimal colors
    """
    if n_categories <= len(PALETTE_CATEGORICAL):
        return PALETTE_CATEGORICAL[:n_categories]
    else:
        # Interpolate more colors
        cmap = plt.cm.get_cmap('tab20')
        return [mpl.colors.rgb2hex(cmap(i / n_categories)) for i in range(n_categories)]


# ============================================================================
# INITIALIZATION
# ============================================================================

# Apply default style on import
# apply_publication_style()

if __name__ == "__main__":
    # Style demo
    print("=" * 60)
    print("PAPER STYLE - Q1 Scopus Configuration Demo")
    print("=" * 60)
    
    apply_publication_style(style='nature')
    
    # Demo figure
    fig, axes = plt.subplots(1, 2, figsize=get_figure_size('double_column'))
    
    # Panel A: Histogram
    np.random.seed(42)
    data = np.random.normal(0, 1, 1000)
    axes[0].hist(data, bins=30, color=COLORS['primary'], edgecolor='white', alpha=0.8)
    style_axis(axes[0], xlabel='Value', ylabel='Frequency', title='Distribution')
    add_panel_label(axes[0], 'A')
    annotate_stats(axes[0], {'n': 1000, 'μ': f'{data.mean():.2f}', 'σ': f'{data.std():.2f}'})
    
    # Panel B: Scatter
    x = np.random.randn(100)
    y = x * 0.5 + np.random.randn(100) * 0.5
    scatter = axes[1].scatter(x, y, c=y, cmap='paper_seq', s=50, alpha=0.8)
    style_axis(axes[1], xlabel='X Variable', ylabel='Y Variable', title='Correlation')
    add_panel_label(axes[1], 'B')
    create_colorbar(scatter, axes[1], 'Value')
    
    plt.suptitle('Publication Quality Figure Demo', fontweight='bold')
    
    demo_path = Path(__file__).parent.parent / 'results' / 'figures' / 'style_demo.png'
    demo_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(demo_path, dpi=300)
    plt.close()
    
    print(f"\n✓ Demo saved at: {demo_path}")
    print("\nAvailable colors:")
    for name, color in list(COLORS.items())[:10]:
        print(f"  • {name}: {color}")
    print("  ...")
    print(f"\nTotal: {len(COLORS)} colors defined")
