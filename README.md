# Matrix Operations with AVX, OpenMP, and MPI

## Overview

This repository contains a set of matrix operations implemented using AVX (Advanced Vector Extensions), OpenMP (Open Multi-Processing), and MPI (Message Passing Interface) to harness the power of High-Performance Computing (HPC) for efficient parallel computation.

## Directory Structure

- **include:** This directory contains header files for each type of matrix.
  
- **src:** The source code for the matrix operations is located in this directory. Different implementations for AVX and OpenMP can be found here.

- **bin:** The compiled binaries for `main.cpp`.

- **MPI:** This directory contains the MPI code.

## Dependencies

- AVX: Ensure your CPU supports AVX instructions.
- OpenMP: Compiler with OpenMP support (G++, Clang).
- MPI: MPI implementation (OpenMPI, MPICH).

To recompile the code:
```bash
    make
    ```
