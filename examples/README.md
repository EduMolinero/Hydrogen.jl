# Examples

The `examples/` folder contains different cases demonstrating the usage of the code.

## Simple High Harmonic Generation (HHG)

This folder contains an input file for a basic High Harmonic Generation (HHG) simulation of the Hydrogen atom. 
The example illustrates the fundamental setup for HHG calculations.

## Optimization

This folder contains code that performs two parallel calculations:

1. One without any Rydberg states between the ground state and the continuum.
2. One with a single Rydberg state at energy \\(E_1\\).

### Purpose

The goal is to optimize the parameters \\(a\\) and \\(\\eta\\) of the soft-core potential to ensure the correct value of the energy levels:

- The ground state energy should match the \\\Delta\\ value provided in both cases.
- The first excited state should either be pushed into the continuum or set to \\(E_1\\).

### Process

The optimization loop follows these steps:

1. Build the Hamiltonian of the system.
2. Compute the eigenvalues for the ground state and first excited state.
3. Check if they match the desired values.
4. If not, adjust the parameters and repeat until convergence is achieved.

This example demonstrates how to control the energy spectrum of the system by tuning the potential parameters.


