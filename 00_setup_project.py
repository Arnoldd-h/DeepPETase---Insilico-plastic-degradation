#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
00_setup_project.py - Initial Project Configuration for Deep-PETase-Mining
================================================================================

Description:
    This script sets up the directory structure for the genomic data mining
    project focused on searching for hypothetical PETases in metagenomes.
    
    The created structure follows best practices for organizing bioinformatics
    projects, separating raw data, processed data, metadata, and results.

Date: December 2025
Project: Deep-PETase-Mining - Phase 1

Usage:
    python 00_setup_project.py

Dependencies:
    - Python 3.9+
    - pathlib (included in standard library)
================================================================================
"""

from pathlib import Path
from typing import List
import logging
from datetime import datetime


def setup_logging() -> logging.Logger:
    """
    Configure the logging system to record script operations.
    
    Returns:
        logging.Logger: Configured logger object for displaying console messages.
    
    Note for students:
        Logging is fundamental in bioinformatics to track what operations
        were performed and when, especially when working with long pipelines.
    """
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger(__name__)


def create_directory_structure(base_path: Path) -> List[Path]:
    """
    Create the necessary directory structure for the project.
    
    Args:
        base_path (Path): Base path where the project structure will be created.
    
    Returns:
        List[Path]: List of successfully created directories.
    
    The created structure is:
        📁 project/
        ├── 📁 data/
        │   ├── 📁 raw/          # Raw FASTA files from NCBI
        │   ├── 📁 processed/    # Filtered and clean sequences
        │   └── 📁 metadata/     # CSVs with organism information
        ├── 📁 src/              # Analysis scripts
        ├── 📁 results/
        │   ├── 📁 figures/      # Exploratory plots
        │   └── 📁 logs/         # Log files
        ├── 📄 README.md
        └── 📄 .gitignore
    """
    # Define all directories that need to be created
    directories: List[Path] = [
        base_path / "data" / "raw",           # Raw FASTA from NCBI
        base_path / "data" / "processed",     # Filtered sequences
        base_path / "data" / "metadata",      # CSVs with metadata
        base_path / "src",                     # Scripts de análisis
        base_path / "results" / "figures",    # Gráficos
        base_path / "results" / "logs",       # Logs de ejecución
    ]
    
    created_dirs: List[Path] = []
    
    for directory in directories:
        try:
            # mkdir with parents=True creates all necessary parent directories
            # exist_ok=True prevents errors if the directory already exists
            directory.mkdir(parents=True, exist_ok=True)
            created_dirs.append(directory)
            logger.info(f"✅ Directory created/verified: {directory}")
        except PermissionError as e:
            logger.error(f"❌ Permission error creating {directory}: {e}")
        except Exception as e:
            logger.error(f"❌ Unexpected error creating {directory}: {e}")
    
    return created_dirs


def create_readme(base_path: Path) -> bool:
    """
    Create a README.md file with basic project information.
    
    Args:
        base_path (Path): Base path of the project.
    
    Returns:
        bool: True if created successfully, False otherwise.
    
    Note for students:
        A good README is your project's cover letter.
        It should explain what it does, how to install it, and how to use it.
    """
    readme_path: Path = base_path / "README.md"
    
    readme_content: str = """# 🧬 Deep-PETase-Mining

## Descripción del Proyecto

Pipeline bioinformático para la minería de datos genómicos enfocado en la 
identificación de enzimas PETasa hipotéticas provenientes de metagenomas de 
ambientes con exposición a plásticos (vertederos, compost, suelos contaminados).

## Estructura del Proyecto

```
Deep-PETase-Mining/
├── data/
│   ├── raw/          # Secuencias FASTA crudas descargadas de NCBI
│   ├── processed/    # Secuencias filtradas y curadas
│   └── metadata/     # Información de organismos y anotaciones
├── src/              # Scripts de análisis
├── results/
│   ├── figures/      # Visualizaciones y gráficos
│   └── logs/         # Registros de ejecución
└── README.md
```

## Requisitos

- Python 3.9+
- Biopython
- Pandas
- Matplotlib/Seaborn
- tqdm

## Instalación

```bash
pip install biopython pandas matplotlib seaborn tqdm
```

## Uso

### Fase 1: Adquisición de Datos

1. **Configurar el proyecto:**
   ```bash
   python 00_setup_project.py
   ```

2. **Minar secuencias de NCBI:**
   ```bash
   python 01_mining_ncbi.py
   ```

3. **Análisis exploratorio:**
   ```bash
   python 02_exploratory_analysis.py
   ```

## Flujo de Trabajo

```mermaid
graph TD
    A[Definir Keywords] --> B[NCBI Entrez API]
    B --> C[Fetch GenBank]
    C --> D[Es Hypothetical Protein?]
    D -->|Si| E[100-2000 aa?]
    D -->|No| F[Descartar]
    E -->|Si| G[candidates.fasta]
    E -->|No| F
    G --> H[metadata.csv]
```

## Autor

[Tu nombre] - Proyecto de Maestría en Bioinformática

## Licencia

MIT License - Ver archivo LICENSE para mas detalles.

---
*Generado automaticamente el """ + datetime.now().strftime("%Y-%m-%d") + "*\n"
    
    try:
        readme_path.write_text(readme_content, encoding='utf-8')
        logger.info(f"✅ README.md created at: {readme_path}")
        return True
    except Exception as e:
        logger.error(f"❌ Error creating README.md: {e}")
        return False


def create_gitignore(base_path: Path) -> bool:
    """
    Create a .gitignore file to exclude unnecessary files from version control.
    
    Args:
        base_path (Path): Base path of the project.
    
    Returns:
        bool: True if created successfully, False otherwise.
    
    Note for students:
        The .gitignore is essential for:
        1. Not uploading large data to GitHub (data/)
        2. Avoiding conflicts with temporary files (__pycache__)
        3. Keeping the repository clean and professional
    """
    gitignore_path: Path = base_path / ".gitignore"
    
    # Content of the .gitignore optimized for bioinformatics projects
    gitignore_content: str = """# ========================================
# .gitignore for Deep-PETase-Mining
# ========================================

# --- Data (very important to exclude) ---
# Raw data can be very large and should be downloaded locally
data/

# --- Python ---
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
build/
develop-eggs/
dist/
downloads/
eggs/
.eggs/
lib/
lib64/
parts/
sdist/
var/
wheels/
*.egg-info/
.installed.cfg
*.egg

# --- Virtual environments ---
.env
.venv
env/
venv/
ENV/
env.bak/
venv.bak/

# --- Jupyter Notebooks ---
.ipynb_checkpoints/
*.ipynb_checkpoints

# --- IDEs ---
.idea/
.vscode/
*.swp
*.swo
*~

# --- Operating system ---
.DS_Store
Thumbs.db
desktop.ini

# --- Temporary files ---
*.tmp
*.temp
*.log
*.bak

# --- Large results (optional, uncomment if needed) ---
# results/figures/*.png
# results/logs/*.log

# --- Sensitive configuration files ---
config.ini
secrets.yaml
*.key
"""
    
    try:
        gitignore_path.write_text(gitignore_content, encoding='utf-8')
        logger.info(f"✅ .gitignore created at: {gitignore_path}")
        return True
    except Exception as e:
        logger.error(f"❌ Error creating .gitignore: {e}")
        return False


def print_project_tree(base_path: Path) -> None:
    """
    Print a visual representation of the created directory tree.
    
    Args:
        base_path (Path): Base path of the project.
    """
    print("\n" + "="*60)
    print("🌳 PROJECT STRUCTURE CREATED:")
    print("="*60)
    
    tree: str = f"""
    {base_path.name}/
    ├── 📁 data/
    │   ├── 📁 raw/          → Raw FASTA sequences
    │   ├── 📁 processed/    → Filtered sequences
    │   └── 📁 metadata/     → CSV files
    ├── 📁 src/              → Analysis scripts
    ├── 📁 results/
    │   ├── 📁 figures/      → PNG plots
    │   └── 📁 logs/         → Log files
    ├── 📄 README.md         → Documentation
    └── 📄 .gitignore        → Git exclusions
    """
    print(tree)
    print("="*60)


def main() -> None:
    """
    Main function that orchestrates project creation.
    
    Executes sequentially:
    1. Directory creation
    2. README.md creation
    3. .gitignore creation
    4. Directory tree visualization
    """
    print("\n" + "🧬"*30)
    print("   DEEP-PETASE-MINING - Project Configuration")
    print("🧬"*30 + "\n")
    
    # Get the path of the directory where this script is located
    # This ensures reproducibility regardless of where it's executed
    base_path: Path = Path(__file__).parent.resolve()
    
    logger.info(f"📂 Project base directory: {base_path}")
    
    # Step 1: Create directory structure
    logger.info("Step 1/3: Creating directory structure...")
    created_dirs = create_directory_structure(base_path)
    
    # Step 2: Create README.md
    logger.info("Step 2/3: Creating README.md...")
    create_readme(base_path)
    
    # Step 3: Create .gitignore
    logger.info("Step 3/3: Creating .gitignore...")
    create_gitignore(base_path)
    
    # Show final result
    print_project_tree(base_path)
    
    # Final summary
    print(f"\n✨ Configuration completed successfully!")
    print(f"📊 Directories created: {len(created_dirs)}")
    print(f"📄 Files created: README.md, .gitignore")
    print(f"\n👉 Next step: Run 'python 01_mining_ncbi.py' to start mining")
    print("\n" + "🧬"*30 + "\n")


# ============================================================================
# ENTRY POINT
# ============================================================================
if __name__ == "__main__":
    # Initialize the logger as a global variable to use in all functions
    logger = setup_logging()
    
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n⚠️ Execution cancelled by user.")
    except Exception as e:
        logging.error(f"❌ Fatal error: {e}")
        raise
