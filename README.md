# 🌊 Landlab to FEniCSx : Finite Element Modeling of Shallow Water Equations 🌊

[![Docker](https://img.shields.io/badge/Docker-Supported-blue?logo=docker)](https://www.docker.com/)
[![Python](https://img.shields.io/badge/Python-3.10-blue?logo=python)](https://www.python.org/)
<a href="https://www.iitm.ac.in/">
  <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/images_logos/iitm_logo.png" alt="FEniCS Logo" width="20">
</a>
<a href="https://flaxandteal.co.uk/">
  <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/images_logos/flaxandteal_logo.jpeg" alt="FEniCS Logo" width="20">
</a>
<a href="https://fenicsproject.org/">
  <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/images_logos/fenics_logo.png" alt="FEniCS Logo" width="20">
</a>

⚡ **Numerical simulation of 2D Shallow Water Equations (SWEs)** using **FEniCS and DOLFINx** on bathymetry data.  
✔️ This repository contains code for solving the **1D & 2D Shallow Water Equations** (SWEs) using `DOLFINx`. The setup involves running the code inside a **Docker container**.  
🚀 **Run simulations in Docker** with both **Continuous Galerkin (CG)** and **Discontinuous Galerkin (DG)** methods.  

![Project Banner](https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/images_logos/projectBanner_gif.gif)

## ✍️ Requirements:
    Docker Desktop
    Dolfinx 0.10.0.0
    numpy
    matplotlib
    imageio


## 🛠️ Installation and Setup

### 👉 Install Docker Desktop
Download and install **Docker Desktop** from the official site:  
🔗 [Docker Installation Guide](https://docs.docker.com/desktop/setup/install/windows-install/)

### 👉 Run DOLFINx in Docker
Open **PowerShell** in **Administrator mode** and run:

```sh
docker run -ti dolfinx/dolfinx:stable   # Stable version
docker run -ti dolfinx/dolfinx:nightly  # Nightly (latest) version
```
This installs DOLFINx in two different modes: stable and nightly.

### 👉 Run Jupyter Notebook in Nightly Mode
To run a Jupyter Notebook using the nightly image in Docker, run:

```sh
docker run --init -ti -p 8888:8888 dolfinx/lab:nightly
```
Then, open the address start with <http://localhost:8888> in the browser.

### 👉 Enter the Running Docker Container Before Cloning the Repo
First, find the container ID using:

```sh
docker ps
```
Then, enter the container (replace <container_id> with the actual ID):

```sh
docker exec -it <container_id> bash
```
### 👉 Install Dependencies
Run the following inside the current root of the container:

```sh
pip install imageio
pip install imageio[ffmpeg]
```
### 👉 Clone and Run the Repository
To clone the repository, run in the current root:

```sh
git clone https://github.com/partha-sakha-paul/FEniCS_Landlab.git
```

## 🐬 Running the Code
The repository contains three main files in different working directories. Each file corresponds to a specific Shallow Water Equations (SWEs in 1D and 2D(CG and DG)) test case.

### ✅ Test 1 : SWEs in 1D -- Nonbreaking Wave Propagation Over a Horizontal Plane

#### Case 1 : In a slightly higher frictional area(`rural case`)
`Directory:` FEniCS_Landlab/SWE_1D_FVMtoFEM

`Run the script:`

```sh
cd FEniCS_Landlab/SWE_1D_FVMtoFEM
python swe1D_Test1_Case1.py
```
#### Case 2 : In a slightly lower frictional area(`urbun case`)
`Directory:` FEniCS_Landlab/SWE_1D_FVMtoFEM

`Run the script:`

```sh
cd FEniCS_Landlab/SWE_1D_FVMtoFEM
python swe1D_Test1_Case2.py
```

### ✅ Test 2 : SWEs 2D - Without Bathymetry experiment(`wetdry = True`) with tidal wave(as in `Test 3 & 4`)
`Directory:` FEniCS_Landlab/SWEs2D_exp_without_bathymetry

`Run the script:`

```sh
cd FEniCS_Landlab/SWEs2D_exp_without_bathymetry
python main.py
```

<h2 align="center">🌊 Simulation Results</h2>

<table align="center">
  <tr>
    <td align="center">
      <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/exp_without_bathymtery_CG/station_(9000%2C%203650)_timeseries_exp_without_bathymetry_CG.png" alt="Output 1" width="250">
      <br><b>Figure 1: At station (9000, 3650)</b>
    </td>
    <td align="center">
      <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/exp_without_bathymtery_CG/station_(11000%2C%203650)_timeseries_exp_without_bathymetry_CG.png" alt="Output 2" width="250">
      <br><b>Figure 2: At station (11000, 3650)</b>
    </td>
    <td align="center">
      <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/exp_without_bathymtery_CG/station_(13500%2C%203650)_timeseries_exp_without_bathymetry_CG.png" alt="Output 3" width="250">
      <br><b>Figure 3: At station (13800, 3650)</b>
    </td>
    <td align="center">
      <video width="250" controls>
        <source src="https://github.com/YOUR_USERNAME/YOUR_REPO/blob/main/videos/simulation.mp4" type="video/mp4">
        Your browser does not support the video tag.
      </video>
      <br><b>Video: Full Simulation</b>
    </td>
  </tr>
</table>


### ✅ Test 3 : SWEs on a 2D Sloped Beach (`Continuous Galerkin`)
`Directory:` FEniCS_Landlab/SWEs_2D_SlopedBeachProblem_SWEmniCS_CG

`Run the script:`

```sh
cd FEniCS_Landlab/SWEs_2D_SlopedBeachProblem_SWEmniCS_CG
python swe2d_sloped_beach_cg_main.py
```

### ✅ Test 4 : SWEs on a 2D Sloped Beach (`Discontinuous Galerkin`)
`Directory:` FEniCS_Landlab/SWEs_2D_SlopedBeachProblem_SWEmniCS_DG

`Run the script:`

```sh
cd FEniCS_Landlab/SWEs_2D_SlopedBeachProblem_SWEmniCS_DG
python swe2d_sloped_beach_dg_main.py
```
