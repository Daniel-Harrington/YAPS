# Galaxy Collision Simulation

![Galaxy Simulation](https://img.shields.io/badge/Simulation-Galaxy%20Collision-blue)
![NVIDIA CUDA](https://img.shields.io/badge/CUDA-NVIDIA-green)
![License](https://img.shields.io/badge/License-MIT-lightgrey)

A high-performance astrophysics simulation modeling black hole migration resulting from the merger of a large host galaxy with a smaller galaxy. The simulation uses **CUDA acceleration** and **Fourier transforms** to compute gravitational interactions between 100 million particles, including two supermassive black holes (SMBHs). The project aims to study the **final parsec problem** in galaxy mergers.

---

## Table of Contents
- [About](#about)
- [Features](#features)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)

---

## About

This project is part of **UTSC PHYD57**, focusing on galaxy evolution simulations. The model includes:

- **100 million particles** representing stars, dark matter, and gas.
- **2 Supermassive Black Holes (SMBHs)** resulted from recently merged galaxes
- **Fourier-based gravitational computations** to resolve forces in O(NLog(N)) computations
- **CUDA-accelerated** particle-grid interpolations and integration.

The simulation aims to observe the behavior of galaxies as they merge, with a particular focus on the **final parsec problem**, where the SMBHs slow down due to dynamical friction.

---

## Features

- **CUDA Acceleration**: All computations are GPU-accelerated
- **FFT-based Gravitational Solvers**: Uses CUFFT to compute accelerations in Fourier space.
- **Integration Method**:Leap Frog 2nd Order Symplectic Integrator
---

## Dependencies

The following tools and libraries are required to build and run the simulation:

- **NVIDIA CUDA Toolkit (>=9.2)**
- **PGI Fortran Compiler**: For compiling CUDA Fortran code.
- **CUFFT Library**: For GPU-accelerated Fourier transforms.

---

## Installation

### Clone the Repository

```bash
git clone https://github.com/Daniel-Harrington/YAPS.git
cd src
```

### Compilation

Use the following command to compile the code:

```bash
pgf95 -mp -Mcuda=fastmath,cc35,cc50,cc60,fma,flushz,lineinfo -ta=nvidia -tp=haswell -O2 -fast -Minfo=all -mcmodel=medium Galaxy_Collison.f95 -o Galaxy_Collision -L/usr/local/cuda-9.2/lib64 -lcufft -lcupti
```

---

## Usage

### Running the Simulation

After successful compilation, run the executable to start the simulation:

```bash
./Galaxy_Collision
```

### Input Parameters

- **Number of particles (N)**: `100,000,000` (adjustable in the source code).
- **Grid dimensions**: Default `nx=512`, `ny=512`, `nz=256`.
- **Timestep (\( \Delta t \))**: Default `1e-6`.
