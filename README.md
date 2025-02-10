# FEniCSx_Landlab
This repository contains code for solving the **1D & 2D Shallow Water Equations** (SWEs) using `DOLFINx`. The setup involves running the code inside a **Docker container**.


## Requirements:
    Docker Desktop
    Dolfinx 0.10.0.0
    Numpy
    Matplotlib


## Installation and Setup

### 1. Install Docker Desktop
Download and install **Docker Desktop** from the official site:  
🔗 [Docker Installation Guide](https://docs.docker.com/desktop/setup/install/windows-install/)

### 2. Run DOLFINx in Docker
Open **PowerShell** in **Administrator mode** and run:

```sh
docker run -ti dolfinx/dolfinx:stable   # Stable version
docker run -ti dolfinx/dolfinx:nightly  # Nightly (latest) version
```
This installs DOLFINx in two different modes: stable and nightly.

### 3. Run Jupyter Notebook in Nightly Mode
To run a Jupyter Notebook using the nightly image in Docker, run:

```sh
docker run --init -ti -p 8888:8888 dolfinx/lab:nightly
```
Then, open the address start with <http://localhost:8888> in the browser.

### 4. Enter the Running Docker Container Before Cloning the Repo
First, find the container ID using:

```sh
docker ps
```
Then, enter the container (replace <container_id> with the actual ID):

```sh
docker exec -it <container_id> bash
```
### 5. Install Dependencies
Run the following inside the current root of the container:

```sh
pip install imageio
pip install imageio[ffmpeg]
```
### 6. Clone and Run the Repository
To clone the repository, run in the current root:

```sh
git clone https://github.com/partha-sakha-paul/FEniCS_Landlab.git
```

## Running the Code
The repository contains three main files in different working directories. Each file corresponds to a specific Shallow Water Equations (SWEs in 1D and 2D(CG and DG)) test case.

### Test 1 : SWEs in 1D -- Nonbreaking Wave Propagation Over a Horizontal Plane

#### Case 1 : In a slightly higher frictional area(rural case)
`Directory:` FEniCS_Landlab/SWE_1D_FVMtoFEM

`Run the script:`

```sh
cd FEniCS_Landlab/SWE_1D_FVMtoFEM
python swe1D_Test1_Case1.py
```
#### Case 2 : In a slightly lower frictional area(urbun case)
`Directory:` FEniCS_Landlab/SWE_1D_FVMtoFEM

`Run the script:`

```sh
cd FEniCS_Landlab/SWE_1D_FVMtoFEM
python swe1D_Test1_Case2.py
```

### Test 2 : SWEs on a 2D Sloped Beach (Continuous Galerkin)
`Directory:` FEniCS_Landlab/SWEs_2D_SlopedBeachProblem_SWEmniCS_CG

`Run the script:`

```sh
cd FEniCS_Landlab/SWEs_2D_SlopedBeachProblem_SWEmniCS_CG
python swe2d_sloped_beach_cg_main.py
```

### Test 3 : SWEs on a 2D Sloped Beach (Discontinuous Galerkin)
`Directory:` FEniCS_Landlab/SWEs_2D_SlopedBeachProblem_SWEmniCS_DG

`Run the script:`

```sh
cd FEniCS_Landlab/SWEs_2D_SlopedBeachProblem_SWEmniCS_DG
python swe2d_sloped_beach_dg_main.py
```
