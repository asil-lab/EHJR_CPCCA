# Cooperative Gaussian Process-based Model Predictive Control for Safe Multi-Agent Navigation

This repository contains the code for reproducing the results in "Cooperative Gaussian Process-based Model Predictive Control for Safe Multi-Agent Navigation".

## Description

The code implements a cooperative control approach for multi-agent systems using Gaussian Process-based Model Predictive Control (GP-MPC) to ensure safe navigation in environments with dynamic obstacles.

- `./DDAAModule`: Contains the core module implementation with algorithms for distributed decision-making and agent coordination.

## Key Points

- The paths specified in the code are relative to the location of the main code file. Adjust paths according to your system setup if needed.
- The main simulation can be run with different parameters. Be sure to update save folders and configuration parameters as necessary.
- This project is implemented in Julia for high-performance numerical computation.

## Support and Questions

For questions and issues, please use the [Issues](../../issues) section of this repository.

## Supported Platforms

- Julia 1.0 and higher

## Citation

```bibtex
@article{Riemens2025,
  author = {Riemens, Ellen H. J. and van der Veen, Alle-Jan and Rajan, Raj T.},
  title = {Cooperative Gaussian Process-based Model Predictive Control for Safe Multi-Agent Navigation},
  journal = {},
  year = {2025}
}
```

## Authors

- **Ellen H. J. Riemens** - E.H.J.Riemens@tudelft.nl
- **Alle-Jan van der Veen** - A.J.vanderVeen@tudelft.nl
- **Raj T. Rajan** - R.T.Rajan@tudelft.nl

Delft University of Technology, Faculty of EEMCS, Mekelweg 4, 2628 CD Delft, The Netherlands

## Getting Started

The code is written in Julia.

### Required Julia Packages

Install the following packages:

- `Statistics`
- `LinearAlgebra`
- `Plots`
- `Distributions`
- `OptimizationProblems`
- `JuMP`
- `MOSEK` (or other optimizer)


### Running Simulations

To run the main simulation:

```julia
julia Test_general.jl
```

