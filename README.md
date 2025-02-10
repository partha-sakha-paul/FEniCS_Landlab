# FEniCSx_Landlab
This repository contains code for solving the **1D & 2D Shallow Water Equations** (SWEs) using `DOLFINx`. The setup involves running the code inside a **Docker container**.


## Requirements:
    Docker Desktop
    Dolfinx 0.10.0.0
    Numpy
    Matplotlib# FEniCS_Landlab - Solving 2D Shallow Water Equations


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
Then, open http://localhost:8888 in your browser.

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
Change to the working directory:

```sh
cd FEniCS_Landlab/SWEs_2D_SlopedBeachProblem_SWEmniCS_CG
```
Run the main script:

```sh
python swe2d_sloped_beach_cg_main.py
```
