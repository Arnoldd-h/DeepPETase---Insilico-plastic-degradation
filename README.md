# 🧬 Deep-PETase-Mining

> **AI-powered discovery of novel plastic-degrading enzymes from environmental metagenomes**

![Graphical Abstract](results/figures/Graphical%20Abstract.png)

---

## 📋 Table of Contents

- [Overview](#overview)
- [Pipeline](#pipeline)
- [Results](#results)
- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Technologies](#technologies)
- [License](#license)

---

## 🔬 Overview

**Deep-PETase-Mining** is a comprehensive bioinformatics pipeline designed to discover novel PET-degrading enzymes (PETases) from environmental metagenomic data using a combination of:

- 🔍 **Data mining** from NCBI protein databases
- 🧬 **Homology filtering** with Smith-Waterman alignment
- 🎯 **Motif scanning** for catalytic sites (G-x-S-x-G)
- 🧠 **AI-powered structure prediction** with ESMFold
- 💊 **Virtual screening** with AutoDock Vina

### 🎯 Key Findings

| Metric | Value |
|--------|-------|
| Initial candidates | 39,000 sequences |
| After homology filtering | 8,500 sequences |
| High-quality structures | 1,718 models |
| Top hits (< -7.0 kcal/mol) | **27 candidates** |
| Best binding affinity | **-7.89 kcal/mol** |

---

## 🔄 Pipeline

The pipeline consists of **5 phases**:

### Phase 1: Data Acquisition
- Mining NCBI protein database for hypothetical proteins
- Keywords: `plastic`, `landfill`, `metagenome`
- Length filter: 100-2000 amino acids

### Phase 2: Homology Filtering
- Reference: IsPETase (*Ideonella sakaiensis*)
- Smith-Waterman local alignment
- Identity threshold: 30-90%
- Catalytic motif scan: G-x-S-x-G

### Phase 3: Structure Prediction
- CD-HIT clustering (80% identity)
- ESMFold structure prediction
- pLDDT quality filtering (≥70)

### Phase 4: Virtual Screening
- Ligand: MHET (PET degradation intermediate)
- Receptor preparation with OpenBabel
- Blind docking with AutoDock Vina

### Phase 5: Analysis & Reporting
- Statistical analysis
- Publication-quality figures
- Top candidates selection

---

## 📊 Results

### Top 10 PETase Candidates

| Rank | Accession | Binding Affinity | pLDDT |
|:----:|-----------|:----------------:|:-----:|
| 1 | GKS65910.1 | -7.89 kcal/mol | 85.2 |
| 2 | PXZ46872.1 | -7.75 kcal/mol | 82.1 |
| 3 | QTV19888.1 | -7.68 kcal/mol | 79.8 |
| 4 | OCN00452.1 | -7.62 kcal/mol | 81.5 |
| 5 | GKS58565.1 | -7.58 kcal/mol | 83.7 |
| 6 | GAB6179983.1 | -7.52 kcal/mol | 80.3 |
| 7 | GAA4408724.1 | -7.48 kcal/mol | 78.9 |
| 8 | KGP88451.1 | -7.45 kcal/mol | 82.4 |
| 9 | QQE08272.1 | -7.41 kcal/mol | 81.1 |
| 10 | EIW90505.1 | -7.38 kcal/mol | 79.6 |

---

## 💻 Installation

### Prerequisites

- Python 3.9+
- AutoDock Vina 1.2.7
- OpenBabel 3.1+
- CD-HIT

### Python Dependencies

```bash
pip install biopython pandas matplotlib seaborn tqdm requests
```

---

## 🚀 Usage

Run scripts sequentially:

```bash
# Phase 1: Data Acquisition
python 00_setup_project.py
python 01_mining_ncbi.py
python 02_exploratory_analysis.py

# Phase 2: Homology Filtering
python 03_get_reference.py
python 04_filter_homology.py
python 05_motif_scan.py
python 06_plot_filtering.py

# Phase 3: Structure Prediction
python 07_cluster_sequences.py
python 08_predict_structures_esm.py
python 09_filter_structures_plddt.py
python 10_plot_plddt.py

# Phase 4: Virtual Screening
python 11_prep_ligand.py
python 12_prep_receptors.py
python 13_run_virtual_screening.py
python 14_select_top_hits.py

# Phase 5: Analysis
python 15_generate_paper_plots.py
```

---

## 📁 Project Structure

```
Deep-PETase-Mining/
├── 00_setup_project.py          # Initialize directories
├── 01_mining_ncbi.py            # NCBI data mining
├── 02_exploratory_analysis.py   # EDA and visualization
├── 03_get_reference.py          # Download IsPETase reference
├── 04_filter_homology.py        # Homology filtering
├── 05_motif_scan.py             # Catalytic motif detection
├── 06_plot_filtering.py         # Filtering visualization
├── 07_cluster_sequences.py      # CD-HIT clustering
├── 08_predict_structures_esm.py # ESMFold prediction
├── 09_filter_structures_plddt.py# Quality filtering
├── 10_plot_plddt.py             # pLDDT visualization
├── 11_prep_ligand.py            # MHET preparation
├── 12_prep_receptors.py         # Receptor conversion
├── 13_run_virtual_screening.py  # AutoDock Vina docking
├── 14_select_top_hits.py        # Top candidates selection
├── 15_generate_paper_plots.py   # Publication figures
├── data/
│   ├── ligands/                 # MHET ligand files
│   ├── metadata/                # Sequence metadata
│   ├── processed/               # Filtered sequences
│   ├── raw/                     # Raw FASTA files
│   └── references/              # IsPETase reference
├── results/
│   ├── docking_results/         # Vina output files
│   ├── figures/                 # Generated plots
│   ├── structures/              # PDB/PDBQT files
│   └── top_hits_pdb/            # Best candidates
└── src/
    └── paper_style.py           # Publication styling
```

---

## 🛠️ Technologies

| Category | Tools |
|----------|-------|
| **Data Mining** | NCBI Entrez API, Biopython |
| **Alignment** | Smith-Waterman (Biopython) |
| **Clustering** | CD-HIT |
| **Structure Prediction** | ESMFold (Meta AI) |
| **Molecular Docking** | AutoDock Vina 1.2.7 |
| **Visualization** | Matplotlib, Seaborn |
| **Data Processing** | Pandas, NumPy |

---

## 📄 License

This project is licensed under the MIT License.

---

## 📧 Contact

For questions or collaborations, please open an issue on this repository.

---

<p align="center">
  <b>🌍 Contributing to the fight against plastic pollution through AI-driven enzyme discovery</b>
</p>