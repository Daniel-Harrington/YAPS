Below is a professional and comprehensive `README.md` file for your project, designed to be copy-paste friendly and formatted for GitHub:

```markdown
# Galaxy Collision Simulation

![Galaxy Simulation](https://img.shields.io/badge/Simulation-Galaxy%20Collision-blue)
![NVIDIA CUDA](https://img.shields.io/badge/CUDA-NVIDIA-green)
![License](https://img.shields.io/badge/License-MIT-lightgrey)

A high-performance astrophysics simulation modeling the collision of a large galaxy with a smaller galaxy. The simulation uses **CUDA acceleration** and **Fourier transforms** to compute gravitational interactions between 100 million particles, including two supermassive black holes (SMBHs). The project aims to study the **final parsec problem** in galaxy mergers.

---

## Table of Contents
- [About](#about)
- [Features](#features)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)
- [Results](#results)
- [Project Structure](#project-structure)
- [Contributing](#contributing)
- [License](#license)

---

## About

This project is part of **UTSC PHYD57**, focusing on galaxy evolution simulations. The model includes:

- **100 million particles** representing stars, dark matter, and gas.
- **2 Supermassive Black Holes (SMBHs)** modeled with enhanced resolution.
- **Fourier-based gravitational computations** to resolve forces with high accuracy.
- CUDA-accelerated particle-grid interpolations and integration.

The simulation aims to observe the behavior of galaxies as they merge, with a particular focus on the **final parsec problem**, where the SMBHs slow down due to dynamical friction.

---

## Features

- **CUDA Acceleration**: All computations are GPU-accelerated for high performance.
- **FFT-based Gravitational Solvers**: Uses CUFFT to compute accelerations in Fourier space.
- **Adaptive Resolution**: Fine resolution around SMBHs and merging regions.
- **Realistic Galaxy Initialization**: Generates logarithmic spirals and realistic orbital parameters for merging galaxies.
- **Energy Conservation**: Implements leapfrog integration with careful energy tracking.

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
git clone https://github.com/your-repo-name/Galaxy-Collision-Simulation.git
cd Galaxy-Collision-Simulation
```

### Compilation

Use the following command to compile the code:

```bash
pgf95 -mp -Mcuda=fastmath,cc35,cc50,cc60,fma,flushz,lineinfo -ta=nvidia -tp=haswell -O2 -fast -Minfo=all -mcmodel=medium Galaxy_Collison.f95 -L/usr/local/cuda-9.2/lib64 -lcufft -lcupti
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
- **Grid dimensions**: Default `nx=128`, `ny=128`, `nz=64`.
- **Timestep (\( \Delta t \))**: Default `1e-6`.

### Output Files

- `gravity_grid_output.csv`: Gravitational accelerations for each timestep.
- `particledata.csv`: Particle positions, velocities, and accelerations.

---

## Results

The simulation visualizes the following:

- **SMBH trajectories** as the galaxies merge.
- **Density profiles** of the galaxies over time.
- **Energy conservation** metrics to verify the accuracy of the simulation.
- The **final parsec** dynamics of SMBHs in the merger.

---

## Project Structure

```
Galaxy-Collision-Simulation/
├── Galaxy_Collision.f95        # Main simulation code
├── particledata.csv            # Initial particle data
├── gravity_grid_output.csv     # Gravitational accelerations
├── README.md                   # Project documentation
└── LICENSE                     # License file
```

---

## Contributing

Contributions are welcome! Please follow these steps:

1. Fork the repository.
2. Create a new branch (`git checkout -b feature-branch`).
3. Commit your changes (`git commit -m "Add new feature"`).
4. Push to the branch (`git push origin feature-branch`).
5. Open a Pull Request.

---

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.
```

### Key Notes:
- Replace `your-repo-name` in the `git clone` command with the actual repository name.
- Update the file paths and results section as you produce outputs.

This README should make your project look professional and accessible to collaborators and evaluators. Let me know if you need further customizations!
