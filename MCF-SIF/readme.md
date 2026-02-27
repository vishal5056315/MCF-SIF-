# MCF-SIF: Mean Curvature Filter-based Structure-aware Image Fusion (MATLAB)

This repository contains the MATLAB implementation of the MCF-SIF method proposed in our manuscript:

> “Photoacoustic Biomedical Signal Fusion for Multi-Patch Based Vascular Details: A Curvature-Based Approach”

The code implements a **two-stage image fusion pipeline** for multi-patch photoacoustic vascular images:

1. **Structure-aware fusion (baseline):**  
   Gradient-based saliency maps and spatially smoothed decision maps are used to select structurally dominant pixels across input patches.

2. **Mean curvature-based weight refinement (MCF):**  
   Binary decision maps are iteratively refined using a mean curvature motion PDE to obtain continuous, geometry-aware weight maps that preserve vessel boundaries while suppressing noise and artifacts.

The implementation targets **photoacoustic vascular images of rat ear** acquired in multiple patches and wavelengths, but the method can be applied to other multimodal or multi-patch imaging scenarios.

---

## Repository Structure

- `main_demo.m`  
  Main entry point.  
  - Demonstrates how to run the MCF-SIF fusion on example photoacoustic image patches.  
  - Calls the fusion routine (`MCF_SIF.m`) and then uses `matrices_new.m` and related helpers to compute quantitative image quality metrics.

- `MCF_SIF.m`  
  Core implementation of the **Mean Curvature Filter-based Structure-aware Image Fusion** method.  
  - Performs normalization, adaptive smoothing, gradient-based saliency computation, spatial consistency (box filtering), binary weight initialization, and iterative mean curvature refinement of the weight maps.

- `mean_curvature_filter.m`  
  Implements the iterative mean curvature motion PDE used to refine the binary weight maps into continuous, geometry-aware weight fields.

- `boxfilter.m`  
  Fast box-filter implementation used for enforcing spatial consistency of decision maps and weight maps.

- `smoothing.m`  
  Local smoothing / adaptive filtering routine (e.g., Wiener-like behavior) used in the preprocessing stage to reduce noise while preserving edges.

- `matrices_new.m`  
  Computes all quantitative image quality metrics used in the paper (e.g., SF, QABF, NABF, entropy and other fusion/quality indices) for the fused outputs and baselines.

- `chi_Sq.m`  
  Performs the chi-square analysis on metric distributions (e.g., to evaluate stability of fusion quality across images or methods).

- `stat.m`  
  Helper functions for summarizing statistics across runs/methods (means, standard deviations, etc.).

- `imhist_fn.m`, `joint_hist_fn.m`  
  Utility functions for computing (joint) histograms, used by some of the information-theoretic and statistical metrics in `matrices_new.m`.

- `readme.md` (this file)  
  High-level description, usage instructions, and citation information for the repository.

---

## Requirements

- MATLAB
- Image Processing Toolbox (for basic image operations and visualization)
- Statistics and Machine Learning Toolbox (recommended for some metrics/statistics, depending on your local modifications)

---

## How to Run

1. Clone or download this repository into a folder on your machine.
2. Add the folder (and its subfolders, if any) to the MATLAB path:
   ```matlab
   addpath(genpath('path_to_this_repo'));
