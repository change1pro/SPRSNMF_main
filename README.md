# SPRSNMF

## Self-Paced Robust Sparse Non-Negative Matrix Factorization with Optimal Graph Learning

Jianyu Miao, Huibo Liu, Liping Yu, Ning Li, and Huiying Hou

*Big Data and Cognitive Computing*, 2026, 10, 379.

[DOI: 10.3390/bdcc10110379](https://doi.org/10.3390/bdcc10110379)

## Abstract

Non-negative matrix factorization (NMF) is a well-established parts-based representation learning paradigm for high-dimensional data. However, existing NMF methods often struggle to capture the intrinsic structure of data and lack sufficient robustness when dealing with corrupted data. To address these limitations, we propose a self-paced robust sparse NMF framework with optimal graph learning, termed SPRSNMF. Specifically, we integrate $\ell_{2,1}$ norm minimization with non-convex sparse regularization into a unified framework. To be specific, the $\ell_{2,1}$ norm-based reconstruction loss is utilized to enhance the robustness against noise and outliers, while $\ell_{2,1/2}$ is imposed on the representation matrix to promote sparsity. Furthermore, the self-paced learning progressively selects samples from “easy” to “hard” so that the influence of noise and outliers is suppressed during optimization. Moreover, by combining adaptive graph learning with the low-rank constraint, we derive an optimal graph learning scheme to effectively uncover the intrinsic structure of data. We develop an effective block coordinate descent (BCD) algorithm to solve the resulting non-convex optimization problem, and theoretically establish its convergence and computational complexity analysis. Extensive experiments on nine benchmark datasets from various domains demonstrate the superiority of SPRSNMF over state-of-the-art approaches.

## Citation

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
