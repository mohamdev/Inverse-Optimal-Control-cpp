# Inverse Optimal Control (IOC) with Bayesian Optimization

## Overview

This repository contains a C++ implementation of **Inverse Optimal Control (IOC)** for box lifting task, leveraging **Bayesian Optimization** using the `bayesopt` implementation. **Direct Optimal Control (DOC)** is implemented through **Pinocchio** implementation of rigid bodies algorithms, **Ipopt** solver, and **CppADCodegen** for gradient computation.

## Features

- **IOC:** Automatically determine the cost function weights for a system by observing optimal behavior.
- **Bayesian Optimization:** Solves IOC by exploring the DOC cost space, in function of costs weights.
- **Direct Optimal Control (DOC):** Formulates the motion generation through the minimization of a weighted sum of cost functions

## Requirements

### Libraries

1. **[Pinocchio](https://github.com/stack-of-tasks/pinocchio):** state-of-the-art Rigid Body Algorithms for poly-articulated systems.
2. **[Ipopt](https://coin-or.github.io/Ipopt/):** An open-source software package for large-scale nonlinear optimization.
4. **[CppADCodegen](https://github.com/joaoleal/CppADCodeGen):** Automatic differentiation with support for code generation.
5. **[BayesOpt](https://github.com/rmcantin/bayesopt):** A library for Bayesian optimization.

### Tools

- **CMake:** For building and managing the project.
- **gcc/g++:** C++ compiler.

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/mohamdev/Inverse-Optimal-Control-cpp.git
cd Inverse-Optimal-Control-cpp
```

### 2. Install Dependencies

Make sure to install all required libraries (Ipopt, CppADCodegen, and BayesOpt) and ensure they are available in your system's library paths.

### 3. Build the Project

```bash
mkdir build && cd build
cmake ..
make
```

## Directory Structure

```plaintext
.
├── apps          # Main executables for testing and demonstrating IOC & DOC
├── cg_libs       # CppADCodegen support files
├── data          # Example data files (inputs, outputs, reference trajectories)
├── include       # Header files
├── output_ioc    # Outputs generated by IOC
├── src           # Source files for IOC, DOC, and supporting functionality
└── CMakeLists.txt # Build configuration
```

## How to Run

### Example: Perform IOC with Bayesian Optimization

1. Navigate to the build directory:
   ```bash
   cd build
   ```
2. Run the IOC executable (replace `example_executable` with the relevant binary name):
   ```bash
   ./apps/main_IOC_bayesopt
   ```

## Acknowledgments

- **Ipopt:** For enabling robust DOC solutions.
- **CppADCodegen:** For making gradient computations fast and efficient.
- **BayesOpt:** For providing a flexible Bayesian optimization framework.

## Contact

For questions or suggestions, feel free to open an issue[ or cont](https://github.com/mohamdev)act [mohamdev](https://github.com/mohamdev).

