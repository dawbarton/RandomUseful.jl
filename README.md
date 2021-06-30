# RandomUseful

[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

A collection of functions that I don't know where to put anywhere else.

* Fourier differentiation matrices - `fourier_diff` for creating differentiation matrices (first and second order) with arbitrary number types on an equispaced grid (assumes periodicity).
* Chebyshev differentiation matrices - `cheb_diff` for creating differentiation matrices (first order) with arbitrary number types on a grid given by the Chebyshev nodes.
* Stability calculations for delay differentiation equations - `dde_stability` using the algorithm of  *Pseudospectral Differencing Methods for Characteristic Roots of Delay Differential Equations*; D. Breda, S. Maset, and R. Vermiglio, SIAM Journal on Scientific Computing 2005 27:2, 482-495.
