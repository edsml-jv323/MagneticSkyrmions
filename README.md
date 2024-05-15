# 2D Spin-Lattice System on HoneyComb and Square/Rectangular Lattices

This package provides tools to simulate and visualise the interactions and behaviors of spins within 2D lattice systems, with a focus on the Metropolis Monte Carlo algorithm for rectangular lattices.

## Foundations: The Heisenberg's Spin Model 👨‍🔬

The core of this package is inspired by the Heisenberg's spin model, adapted and simplified for this work. The primary assumptions are:

- **Two dimensional lattice:** Atoms are arranged in a two-dimensional lattice.
- **Fixed Atomic Position:** The atoms maintain a fixed position in space, meaning their distance from their nearest neigbours remains constant.
- **Magnetic Atoms:** Every atom is magnetic, so they have a spin vector `s=(Sx, Sy, Sz)`.
- **Constant Spin Magintude:** The magnitude of all spins is conserved \(|s| = 1\).

## Overview

- **Square/Rectangular Lattice Dynamics:**
    - Simulates how spins interact within a rectangular lattice, defined by dimensions `nx` and `ny`.
    - **Energy Calculations:** Includes methods to compute various energy contributions such as Zeeman, Anisotropy, Exchange, and DMI.
    - **Metropolis Algorithm:** Implements the Metropolis Monte Carlo technique to achieve spin configurations that represent equilibrium states. The implementation provided accepts changes that lower the energy but does not include probabilistic acceptance based on the Boltzmann factor. For simplicity T=0 was considered.
    - **Visualisation:** Offers visualisation tools to showcase the spin configurations, useful to understand their orientation and magnitude after applying Metropolis Algorithm.
    Notably, the visualisation capabilities extend to visualising **skyrmions** - topological spin textures renowned for their potential in spintronic applications. Observing these entities in the lattice systems offers valuable insights into their behavior, formation, and interactions.

- **Honeycomb Lattice Dynamics** 🐝
    - **Dual Atom Representation**: Recognises that every unit cell contains two different atoms, labelled A and B, each with its own spin vector.
    - **Mean Spin:** Calculates the average spin value across the entire lattice, useful for understanding general orientations.
    - **Magnetisation:** Evaluates the net magnetization of the lattice, providing insights into collective magnetic behavior.
    - **Visualisation:** Visualise random or manually set spin configurations in the unique honeycomb lattice.

## Getting started guide

- **Step 1**: Create a fresh new conda environment with Python 3.11 as it follows:
`conda create -n mcsim python=3.11`
- **Step 2**: Activate this new conda environment: `conda activate mcsim`
- **Step 3**: Clone the repository by running: `git clone url_of_the_git_repository`
- **Step 4**: Navigate to where the repository is cloned: `cd repository_name`
- **Step 5**: Install the `mcsim` package by running: `pip install .`

Now you can run all the methods 😊

## Usage

Here's a basic guide on how to use the package:

**Initialization**

```python
import mcsim

#Create a 2D square/rectangular lattice with random spins

spins = mcsim.Spins(n=(10,10))
s.randomise()

#Tune the parameters for calculating the energies
system = mcsim.System(s=s, B=(0, 0, 0.1), K=0.01, u=(0, 0, 1), J=0.5, D=0.5)
```
**Monte Carlo simulation**

```python
driver = mcsim.Driver()
driver.drive(system, n=10_000) #Select number of iterations
```

**Analysis**

```python
system.s.mean #Displays the mean spin of the final state
system.s.plot() #plots the spin configuration in the final state
```

**For the Honeycomb lattice**
```python
s_honey = mcsim.HoneycombSpins((5, 5))
s_honey.randomise()
print("Mean Spin: ", s_honey.mean())
print("Net Magnetization: ", s_honey.magnetization())
s_honey.plot_honeycomb()
```

## Understanding the Code 💡
For those who are new to `NumPy` or prefer a more step-by-step approach, each function's docstring contains an "unoptimised" version of the algorithm. This version uses basic `for` loops to detail the operations.

This is intended to provide a clearer picture of the logic behind each function. Once you're comfortable with the basics, you can then delve into the optimised `NumPy` version to appreciate the efficiencies gained.

We recommend new users or those unfamiliar with `NumPy` to check these unoptimised codes first.

## ✅ Testing

This project uses `pytest` for unit testing to ensure code quality and functionality. Regular testing helps in identifying potential issues early and ensures robustness of the codebase. 

We aim to maintain high test coverage to ensure the tool's reliability. If you're contributing to the project, please ensure your contirbutions are well-tested, and maintain or improve the current test coverage 😄




## 🚧Development Status🚧
Please note that this code is still in active development and might undergo significant changes.

**Honeycomb Lattice:** As of the current version, the honeycomb lattice primarily offers visualisation for random spins. We have plans in our roadmap to extend functionalities, including energy computations and Monte Carlo simulations, for the honeycomb lattice in future releases.

**Probabilistic Monte Carlo:** In future updates, we plan to enhance the Monte Carlo simulations to incorporate probabilistic acceptance of spin configurations. This approach will be based on the Boltzmann factor, allowing for more realistic simulations of systems at finite temperature.

**Future additions:** We are also considering the implementation of the Ising model on a honeycomb lattice. Stay tuned for more features and advancements in our subsequent releases!

## License
Distributed under the Apache License. See `license.md` for more information.

## Contact
Jorge Veiras Yanes - jv323@ic.ac.uk
