# Origin
The codes presented in this repository originate from @mafaldaramoa.
The orignial repository is: https://github.com/mafaldaramoa/ceo-adapt-vqe.git

# ADAPT-VQE Simulation Code

This repository contains code to simulate the OO-OP, OO-RP, RO-OP and RO-RP variants of the Adaptive Derivative-Assembled Problem-Tailored (ADAPT) - Variational Quantum Eigensolver (VQE). In particular, this code was used in the following paper:
* [Strategies for Overcoming Gradient Troughs in the ADAPT-VQE Algorithm](https://arxiv.org/abs/2512.25004)

The original code by @mafaldaramoa with a implementation for a wide array of ADAPT-VQE variants and pools was used in the following papers:
* [Reducing the Resources Required by ADAPT-VQE Using Coupled Exchange Operators and Improved Subroutines](https://arxiv.org/abs/2407.08696)
* [Reducing measurement costs by recycling the Hessian in adaptive variational quantum algorithms](https://arxiv.org/abs/2401.05172)

## Installation

You can install `arXiv:2512.25004_CODES` as follows:

Use anaconda to install the necessary dependencies and a proper environment.
```
conda env create -f xO-yP_env.yml
```
The xO-yP_env.yml file can be found in the folder [arXiv:2512.25004_CODES](https://github.com/jonas-stadelmann/xO-yP_adapt-vqe/tree/b990b1282186af6f918c89ff0b5092e2019abce9/arXiv%3A2512.25004_CODES) within the repository fork.

This ensures the installation of a proper environment with all required packages and installations. This conda environment uses Python 3.11.11 in our case.

Note that PySCF does not support Windows. You can use Windows Subsystem for Linux (WSL) to install a Linux distribution (e.g. Ubuntu), then install Anaconda.

The results within the paper [Strategies for Overcoming Gradient Troughs in the ADAPT-VQE Algorithm](https://arxiv.org/abs/2512.25004) were obtained with the following package versions:

```
qiskit 2.0.0
pyscf 2.8.0
openfermion 1.7.0
openfermionpyscf 0.5
scipy 1.15.2
numpy 1.26.4
```

!!! While the new implementations for adding operators in different positions as stated in [Strategies for Overcoming Gradient Troughs in the ADAPT-VQE Algorithm](https://arxiv.org/abs/2512.25004) were in theory implemented for all versions of ADAPT-VQE within the original code of @mafaldaramoa, I only tested and made sure for the scheme to work for select_via_gradient as this was the main focus of the paper. I cannot guarantee the scheme to work for other implementations of ADAPT-VQE within the code. It is incumbent upon the user to verify this. !!!

## Test Systems

All example scripts use the $H_6$ molecule with $4 \text{ Å}$ as investigated in the paper. Other molecules can be implemented in the submodule molecules.

## Simulation Time

For larger molecules, such as $H_6$, simulations might take several hours to complete. To speed up simulations, you may create an eigendecomposition of the pool (see method `create_eig_decomp()` in submodule `pools` and method `load()` in `algorithms.adapt_vqe`). While the eigendecomposition itself takes hours to compute, once it is created it can be used for the simulation of any system with the same number of qubits (with the same pool).

## Supported Variants

** While the new implementations for adding operators in different positions as stated in [1] were in theory implemented for all versions of ADAPT-VQE within the original code of @mafaldaramoa, I only tested and made sure for the scheme to work for select_via_gradient as this was the main focus of the paper. I cannot guarantee the scheme to work for other implementations of ADAPT-VQE within the code. It is incumbent upon the user to verify this.

For all options regarding the ADAPT-VQE implementation, see `AdaptVQE` class constructor in `algorithms.adaptvqe`. The current implemention supports (**) Hessian recycling [3], TETRIS [5] and orbital optimization [8], as well as a variety of selection and convergence criteria.

A variety of pool options are also supported, namely all CEO variants (OVP, MVP, DVG, DVE) [2], the qubit pool [6], the QE pool [7], and fermionic pools - GSD, SD, Spin-Adapted GSD, etc [4]. For details, see submodule `pools`.

## References

[1] [Strategies for Overcoming Gradient Troughs in the ADAPT-VQE Algorithm](https://arxiv.org/abs/2512.25004)

[2] [Reducing the Resources Required by ADAPT-VQE Using Coupled Exchange Operators and Improved Subroutines](https://arxiv.org/abs/2407.08696)

[3] [Reducing measurement costs by recycling the Hessian in adaptive variational quantum algorithms](https://arxiv.org/abs/2401.05172)

[4] [An adaptive variational algorithm for exact molecular simulations on a quantum compute](https://www.nature.com/articles/s41467-019-10988-2)

[5] [TETRIS-ADAPT-VQE: An adaptive algorithm that yields shallower, denser circuit ansätze](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.6.013254)

[6] [Qubit-ADAPT-VQE: An Adaptive Algorithm for Constructing Hardware-Efficient Ansätze on a Quantum Processor](https://journals.aps.org/prxquantum/abstract/10.1103/PRXQuantum.2.020310)

[7] [Qubit-excitation-based adaptive variational quantum eigensolver](https://www.nature.com/articles/s42005-021-00730-0)

[8] [Self-Consistent Field Approach for the Variational Quantum Eigensolver: Orbital Optimization Goes Adaptive](https://pubs.acs.org/doi/10.1021/acs.jpca.3c05882)
