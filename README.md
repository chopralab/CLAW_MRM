# CLAW-MRM

**CLAW-MRM** (Comprehensive Lipidomics Automation Workflow) - a pipeline for processing and analyzing MRM (Multiple Reaction Monitoring) lipidomics data.

**Created by:** Sanjay Iyer

## Overview

CLAW is designed to streamline the analysis of MRM-based lipidomics experiments. In this tutorial we will parse mzML data and format in organized pandas DataFrames and export as csv files. We will also reformat this data for various tools such as EdgeR (differential expression) and pathway analysis.

## Repository Structure

```
CLAW_tutorial/
├── CLAW_tutorial.ipynb          # Tutorial Jupyter notebook
├── README.md
├── scripts/                     # Core analysis scripts
│   ├── CLAW.py                  # Main CLAW processing script
│   ├── CLAW_preedgeR.py         # Pre-EdgeR formatting script
│   └── CLAW_pathwayanalysis.py  # Pathway analysis formatting
├── lipid_database/              # Reference databases
│   ├── Custom_MRM.csv           # Custom MRM transitions
│   └── Lipid_Database.xlsx      # Lipid reference database
├── projects/                    # Project directories
│   └── lipid_load/              # Example project
│       ├── mzml/                # mzML data files
│       ├── labels/              # Sample metadata & labels
│       ├── results/             # CLAW output results
│       ├── pre_edger/           # EdgeR-formatted data
│       └── pathway_analysis/    # Pathway analysis formatted data
│     
├── CLAW_MRM/                    # CLAW MRM source code
```

## Getting Started

### Prerequisites

- Python 3.9+ (recommend using Python 3.11+ for best performance)
- Required packages: `pandas`, `numpy`, `pymzml`, `matplotlib`, `openpyxl`

### Environment Setup

Create a conda environment with all required dependencies:

```bash
# Using default name "CLAW"
conda env create -f requirements/CLAW.yml

# Or specify a custom environment name (make any name you want in place of my_custom_name)
conda env create -f requirements/CLAW.yml -n my_custom_name 

# Activate the environment
conda activate CLAW  # or your custom name
```

### Quick Start

1. Open `CLAW_tutorial.ipynb` in Jupyter
2. Follow the step-by-step instructions for:
   - Loading and processing mzML files
   - Matching lipids against the custom MRM database
   - Generating pre-EdgeR formatted files
   - Creating Pathway Analysis input files


## Scripts

| Script | Description |
|--------|-------------|
| `CLAW.py` | Main processing script for mzML files and lipid matching |
| `CLAW_preedgeR.py` | Formats intensity data for EdgeR differential analysis |
| `CLAW_pathwayanalysis.py` | Generates pathway analysis formatted output |

## Contact

For questions or collaboration:

- **Sanjay Iyer** - iyer95@purdue.edu
- **Chopra Lab** - choprait@purdue.edu