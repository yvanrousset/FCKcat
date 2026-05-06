# FCKcat

This repository accompanies the manuscript:

**Overcoming systematic data biases enables accurate prediction of enzyme k_cat fold-changes for computational protein design**

Yvan Rousset, Alexander Kroll, and Martin J. Lercher  
*bioRxiv* 2026 — https://doi.org/10.64898/2026.01.23.701068

## Status

Code is available here. Data are **not yet publicly released** and will be deposited on **Zenodo upon submission** of the manuscript (link will be added here).

> To run the notebooks, download the data archive from Zenodo and place its contents in the `data/` folder.

## Repository structure

```
FCKcat/
├── notebooks/
│   ├── nb1_bias_analysis.ipynb              # Variance partitioning (Figures 1, S1)
│   ├── nb2_model_predictions.ipynb          # FCKcat predictions on test pairs (Figures 3, S5, S6)
│   ├── nb3_model_comparison.ipynb           # Comparison with DLKcat and EITLEM (Figure S2)
│   ├── nb4_model_selection.ipynb            # Model selection and ablations
│   └── utils/
│       ├── helpers.py                       # Shared plotting and preprocessing utilities
│       └── CCB_plot_style_0v4.mplstyle      # Matplotlib style sheet
├── data/                                    # Empty — download from Zenodo
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

Download the data archive from Zenodo *(link to be added)* and extract it into the `data/` folder.

### 3. Run the notebooks

Notebooks are numbered in logical order:

| Notebook | Description | Figures |
|---|---|---|
| `nb1_bias_analysis.ipynb` | Variance partitioning of kcat across enzyme groups using linear mixed-effects models | 1, S1 |
| `nb2_model_predictions.ipynb` | FCKcat Δkcat predictions on all test pairs; collective absolute kcat estimator | 3, S5, S6 |
| `nb3_model_comparison.ipynb` | Head-to-head comparison with DLKcat and EITLEM on unseen sequences | S2 |
| `nb4_model_selection.ipynb` | Model selection and feature ablation experiments | — |

## Citation

If you use this repository, please cite:

> Rousset Y., Kroll A., Lercher M.J. (2026). *Overcoming systematic data biases enables accurate prediction of enzyme k_cat fold-changes for computational protein design.* bioRxiv. https://doi.org/10.64898/2026.01.23.701068

The Zenodo data record will be linked here upon release.

## Contact

For questions, please open an issue in this repository.
