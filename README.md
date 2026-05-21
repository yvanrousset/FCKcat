# FCKcat

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20325541.svg)](https://doi.org/10.5281/zenodo.20325541)

This repository accompanies the manuscript:

**Overcoming systematic data biases enables accurate prediction of enzyme k_cat fold-changes for computational protein design**

Yvan Rousset, Alexander Kroll, and Martin J. Lercher  
*bioRxiv* 2026 — https://doi.org/10.64898/2026.01.23.701068

## Status

Code is available in this repository. The data and model files required to reproduce the analyses are available on Zenodo:

https://doi.org/10.5281/zenodo.20325541

To run the notebooks, download `data.zip` from Zenodo and extract it into the repository root. This should create a `data/` folder with the required datasets, models, mappings, baseline predictions, and results.

## Repository structure

```text
FCKcat/
├── notebooks/
│   ├── nb1_bias_analysis.ipynb              # Variance partitioning (Figures 1, S1)
│   ├── nb2_model_predictions.ipynb          # FCKcat predictions on test pairs (Figures 3, S5, S6)
│   ├── nb3_model_comparison.ipynb           # Comparison with DLKcat and EITLEM (Figure S2)
│   ├── nb4_model_selection.ipynb            # Model selection and ablations
│   └── utils/
│       ├── helpers.py                       # Shared plotting and preprocessing utilities
│       └── CCB_plot_style_0v4.mplstyle      # Matplotlib style sheet
├── data/                                    # Created after extracting data.zip from Zenodo
│   ├── datasets/                            # Train/test pair dataframes
│   ├── models/                              # Trained FCKcat model (XGBoost)
│   ├── mapping/                             # ESM2 embeddings and sequence index maps
│   ├── EITLEM_data/                         # EITLEM baseline predictions and metadata
│   ├── DLKcat_data/                         # DLKcat baseline predictions and metadata
│   └── results/                             # Aggregated evaluation results
├── requirements.txt
└── README.md
```

## Reproducing the analyses

### 1. Install dependencies

```bash
pip install -r requirements.txt
```

### 2. Download the data

Download `data.zip` from Zenodo:

https://doi.org/10.5281/zenodo.20325541

Extract it into the root of this repository:

```bash
unzip data.zip
```

After extraction, the repository should contain:

```text
FCKcat/
├── data/
│   ├── datasets/
│   ├── models/
│   ├── mapping/
│   ├── EITLEM_data/
│   ├── DLKcat_data/
│   └── results/
├── notebooks/
├── requirements.txt
└── README.md
```

### 3. Run the notebooks

Notebooks are numbered in logical order:

| Notebook | Description | Figures |
|---|---|---|
| `nb1_bias_analysis.ipynb` | Variance partitioning of kcat across enzyme groups using linear mixed-effects models | 1, S1 |
| `nb2_model_predictions.ipynb` | FCKcat Δkcat predictions on all test pairs; collective absolute kcat estimator | 3, S4, S5, S6 |
| `nb3_model_comparison.ipynb` | Head-to-head comparison with DLKcat and EITLEM on unseen sequences | 4, S2 |
| `nb4_model_selection.ipynb` | Model selection and feature ablation experiments | S3 |

## Data availability

The data and model files required to reproduce the analyses are available on Zenodo:

> Rousset, Y. (2026). *Data and model files for FCKcat: prediction of enzyme kcat fold-changes*. Zenodo. https://doi.org/10.5281/zenodo.20325541

The Zenodo archive contains the `data/` directory expected by the notebooks, including train/test dataframes, the trained FCKcat XGBoost model, ESM2 embeddings and mapping files, DLKcat and EITLEM baseline data, and aggregated evaluation results.

## Citation

If you use this repository, please cite:

> Rousset Y., Kroll A., Lercher M.J. (2026). *Overcoming systematic data biases enables accurate prediction of enzyme k_cat fold-changes for computational protein design.* bioRxiv. https://doi.org/10.64898/2026.01.23.701068

## Contact

For questions, please open an issue in this repository.
