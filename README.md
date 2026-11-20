# SPRSNMF

MATLAB code associated with **Self-Paced Robust Sparse Non-Negative Matrix Factorization with Optimal Graph Learning**.

Jianyu Miao, Huibo Liu, Liping Yu, Ning Li, and Huiying Hou. *Big Data and Cognitive Computing*, 2026, 10, 379.  
Paper: [https://doi.org/10.3390/bdcc10110379](https://doi.org/10.3390/bdcc10110379)

## Overview

SPRSNMF is an unsupervised non-negative matrix factorization method for clustering. It combines robust reconstruction, self-paced sample weighting, adaptive graph learning, and a non-convex sparsity constraint on the coefficient matrix. The graph and low-dimensional representation are learned together.

For a data matrix `X` with features in rows and samples in columns, the model learns a non-negative basis and a non-negative sample representation. Cluster assignments are obtained from the learned representation. Ground-truth labels are used to evaluate clustering results.

## Repository structure

```text
SPRSNMF_main/
├── core/                         Core factorization and update functions
├── dataset/                      MATLAB benchmark datasets
├── utils/                        Data loading, normalization, and evaluation
├── SPRSNMF_results_dn.m           Parameter search and clustering evaluation
├── SPRSNMF_ablation.m             Ablation experiment driver
├── SPRSNMF_convergence.m          Convergence plots from saved results
├── SPRSNMF_sensitivity.m          Parameter sensitivity plots from saved tables
├── SPRSNMF_sample_weight.m        Self-paced growth-step comparison plots
├── clustering_visualization.m    Original-data and representation visualization
├── robustness_gaussian.m          Gaussian-noise robustness plots from saved tables
├── obj_g.m                       Objective-curve experiment script
└── obj_f.m                       Plots saved objective curves
```

Some development identifiers and output paths retain the name `SPLAGNMF` or `splagnmf`.

## Requirements

- MATLAB.
- Statistics and Machine Learning Toolbox for the core initialization with `kmeans` and the visualization scripts using `tsne`.
- Saved experiment outputs for scripts that read result tables or objective curves. These outputs are not included in this repository.

## Datasets

The `dataset/` directory contains:

`att`, `BA`, `COIL20`, `colon`, `GLIOMA`, `JAFFE`, `lung`, `lymphoma`, `oh15`, `TOX_171`, `UMIST`, and `warpPIE10P`.

The data loader expects a variable named `X` and a label vector named `Y`, `y`, or `gnd`. It returns the data as a features-by-samples matrix, makes the data non-negative, and normalizes each sample to unit Euclidean norm.

From the repository root, a data-loading example is:

```matlab
rootDir = pwd;
addpath(fullfile(rootDir, 'core'));
addpath(fullfile(rootDir, 'utils'));

[X, Y, info] = load_dataset(fullfile(rootDir, 'dataset'), 'att');
```

## Preparing the current code snapshot

The uploaded files retain their original development configuration. Before running the experiment scripts, address the following:

1. **Function name:** `core/SPRSNMF.m` currently declares the function `splagnmf`, and the experiment scripts call `splagnmf`. Align the file name and callers; for example, rename this file to `core/splagnmf.m` in your local checkout.
2. **Core helpers:** `pairwise_sqdist.m` and `project_simplex.m` are referenced by the core functions but are not included in the current snapshot. Add these dependencies to the MATLAB path before factorization.
3. **Dataset paths:** several scripts pass the repository root to `load_dataset`, while the bundled data are in `dataset/`. Pass `fullfile(rootDir, 'dataset')` as the loader's first argument.
4. **Dataset selection:** `SPRSNMF_results_dn.m` currently selects `YaleB_32x32`, which is not bundled. Select an included dataset, such as `att`, or supply the selected dataset separately.
5. **Experiment-specific inputs:** update the local path in `SPRSNMF_convergence.m`. The ablation driver additionally references `splagnmf_without_graph`, `splagnmf_without_sparse`, and `splagnmf_without_spl`, which are not included. Plotting scripts require the result files named in their headers.

## Model interface

After the preparation above, the core function is called as:

```matlab
[F, obj] = splagnmf(X, c, alpha, lambda, beta, k);
pred = cluster_embedding(F, c);
metrics = clustering_metrics(Y, pred);
```

| Input/output | Meaning |
| --- | --- |
| `X` | Non-negative data matrix, features by samples |
| `c` | Number of clusters / factorization rank |
| `alpha` | Weight of the representation's graph regularization term |
| `lambda` | Weight of the spectral graph-learning term |
| `beta` | Weight of the sparsity term |
| `k` | Neighborhood-size parameter used in graph learning |
| `F` | Non-negative sample representation, samples by clusters |
| `obj` | Objective values recorded during iterations |

The evaluation utility reports ACC, NMI, and Purity. Hyperparameters are dataset dependent. The parameter-search scripts contain large grids; reduce the parameter lists and repeat counts for an initial check.

## Citation

If you use this work, please cite:

```bibtex
@article{miao2026sprsnmf,
  title   = {Self-Paced Robust Sparse Non-Negative Matrix Factorization with Optimal Graph Learning},
  author  = {Miao, Jianyu and Liu, Huibo and Yu, Liping and Li, Ning and Hou, Huiying},
  journal = {Big Data and Cognitive Computing},
  year    = {2026},
  volume  = {10},
  pages   = {379},
  doi     = {10.3390/bdcc10110379}
}
```
