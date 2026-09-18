# Bayesian Optimisation Framework for UT&C-BEM

This repository contains a MATLAB framework for Bayesian multi-objective optimisation using the Urban Tethys-Chloris Building Energy Model (UT&C-BEM).

The framework combines Gaussian process surrogate modelling with Bayesian optimisation and augmented Tchebycheff scalarisation to efficiently explore a computationally expensive multi-objective parameter space.

## Requirements

- MATLAB R2025b
- Statistics and Machine Learning Toolbox

The repository contains the modified UT&C-BEM implementation required to run the framework, so no separate installation of UT&C-BEM is required.

## Getting Started

Clone or download this repository and open it in MATLAB.

Add the repository and all subfolders to the MATLAB path:

```matlab
addpath(genpath(pwd));
```

Before running the framework, open `setup.m` and check the paths and model/optimisation settings.

## Main Workflow

The main workflow consists of three consecutive steps.

### 1. Generate the initial LHS dataset

Run:

```matlab
generate_initial_LHS_dataset
```

This script generates the initial Latin Hypercube Sampling (LHS) design and evaluates the corresponding parameter combinations using UT&C-BEM.

The resulting dataset provides the initial observations for the subsequent Bayesian optimisation.

### 2. Run ParEGO Bayesian optimisation

Run:

```matlab
run_parego_BO_10runs
```

This script performs 10 independent ParEGO Bayesian optimisation runs.

The optimisation uses Gaussian process surrogate models and augmented Tchebycheff scalarisation to explore different regions of the multi-objective solution space.

The results of the individual optimisation runs are stored for subsequent analysis.

### 3. Analyse the Pareto front

Run:

```matlab
dense_GP_pareto_all_runs
```

This script combines the results from the Bayesian optimisation runs and uses Gaussian process surrogate models to densely explore the objective space and identify the resulting Pareto front.

This step is intended to be performed after the ParEGO optimisation runs.

## LHS Baseline and Chebyshev Densification

As a separate workflow, the repository also contains:

```matlab
lhs_baseline_and_chebyshev_densification
```

This script provides a baseline and comparison approach based on the initial LHS design and Chebyshev-based densification.

It can be run independently if the required input data are available. Its main purpose is to provide a comparison with the solutions obtained from the Bayesian optimisation workflow described above.

## Optimisation Objectives

The framework is designed for multi-objective optimisation of urban climate and thermal-comfort metrics.

For the application presented in the associated study, the objectives are:

- Peak air temperature at 2 m (`T_peak`)
- Peak Universal Thermal Climate Index at 1.1 m (`UTCI_peak`)
- Minimum nocturnal air temperature following the hottest day (`T_night`)

All objectives are formulated as minimisation objectives.

## Bayesian Optimisation Approach

The framework uses Gaussian process surrogate models to approximate the computationally expensive UT&C-BEM simulations.

ParEGO is used to transform the multi-objective optimisation problem into a sequence of scalar optimisation problems using augmented Tchebycheff scalarisation. Different scalarisation weights allow different regions of the Pareto front to be explored across optimisation runs.

## Repository Structure

```text
Bayesian-Optimization-Framework-for-UT-C/
│
├── generate_initial_LHS_dataset.m
├── run_parego_BO_10runs.m
├── dense_GP_pareto_all_runs.m
│
├── lhs_baseline_and_chebyshev_densification.m
│
├── setup.m
├── getFixedRefPoint.m
├── UTC.m
│
├── ...
│
├── README.md
└── LICENSE
```

## Reproducibility

The code provided in this repository corresponds to the implementation used for the associated research study.

For reproducibility, use the specific repository version associated with the corresponding Zenodo release.

## UT&C-BEM

This repository contains a modified version of the Urban Tethys-Chloris Building Energy Model (UT&C-BEM), originally developed by Naika Meili, Simone Fatichi and contributors.

The UT&C-BEM implementation included in this repository has been modified for integration with the Bayesian optimisation framework.

The original UT&C-BEM software release is available at:

https://doi.org/10.5281/zenodo.14788637

Please cite the original UT&C-BEM publication and software release when using the model.

## Citation

If you use this framework in your research, please cite the associated publication and the corresponding Zenodo release.

publication citation will be added upon acceptance of the paper

[(https://doi.org/10.5281/zenodo.22706138)]

## Licence

The Bayesian optimisation framework developed for this repository is released under the MIT License.

The repository also contains modified UT&C-BEM code. The UT&C-BEM code remains subject to the usage conditions specified by its original authors. Please refer to the UT&C-BEM section above and the original software release for the applicable terms and citation requirements.
