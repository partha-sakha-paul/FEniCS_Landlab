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
📒The repository contains five main files in different working directories. Each file corresponds to a specific Shallow Water Equations (SWEs in 1D and 2D(CG and DG)) test case.

### ✅ Test 1 : SWEs in 1D -- Nonbreaking Wave Propagation Over a Horizontal Plane

#### Case 1 : In a slightly higher frictional area(`rural case`)
`Directory:` FEniCS_Landlab/SWE_1D_FVMtoFEM

`Details:` This experiment consists of a wave propagation over a horizontal plane that is initially dry by assuming a constant velocity distribution over the direction of wave propagation with the moving boundary condition h(ut, t) = 0.

    Parameters: n = 0.01(higher friction), u_const = 0.4(constant velocity in x direction) -- Figure 1

`Run the script:`

```sh
cd FEniCS_Landlab/SWE_1D_FVMtoFEM
python swe1D_Test1_Case1.py
```
#### Case 2 : In a slightly lower frictional area(`urbun case`)
`Directory:` FEniCS_Landlab/SWE_1D_FVMtoFEM

`Details:` Same as Case 1 but with different parameters.

    Parameters: n = 0.005(lower friction), u_const = 0.635(constant velocity in x direction) -- Figure 2

`Run the script:`

```sh
cd FEniCS_Landlab/SWE_1D_FVMtoFEM
python swe1D_Test1_Case2.py
```

<h2 align="center">🌊 Simulation Results</h2>

<table align="center">
  <tr>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs1D/water_depth_evolution_1D_Case1.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs1D/water_depth_evolution_1D_Case1.png" alt="Output 1" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs1D/water_depth_evolution_1D_Case1.png">Figure 1: Case-1 - Higher Frictional Area(Rural case)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs1D/water_depth_evolution_1D_Case2.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs1D/water_depth_evolution_1D_Case2.png" alt="Output 2" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs1D/water_depth_evolution_1D_Case2.png">Figure 2: Case-2 - Lower Frictional Area(Urban case)</a></b>
  </tr>
</table>

### ✅ Test 2 : SWEs 2D - Without Bathymetry experiment(`wetdry = True`) with tidal wave(as in `Test 3 & 4`)
`Directory:` FEniCS_Landlab/SWEs2D_exp_without_bathymetry

`Details:` This experiment evaluates SWEMniCS' wetting and drying capability using Kärna’s α scheme with a `CG` approach on a `13,800 m × 7,200 m` water surface without bathymetry. The `initial` water surface is `flat`, and a `harmonic tide (2m amplitude, 12h period)` is applied at the `open left` boundary, while `other` boundaries are `walls`. The simulation runs for `7 days`, using Manning’s friction `0.02` s/m<sup>1/3</sup>, and records elevations and velocities at `x = 9000 m, 11,000 m, and 13,500 m.`  --  `(Figure 3, 4, 5)`

`Run the script:`

```sh
cd FEniCS_Landlab/SWEs2D_exp_without_bathymetry
python main.py
```

<h2 align="center">🌊 Simulation Results</h2>

<table align="center">
  <tr>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/exp_without_bathymtery_CG/station_(9000%2C%203650)_timeseries_exp_without_bathymetry_CG.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/exp_without_bathymtery_CG/station_(9000%2C%203650)_timeseries_exp_without_bathymetry_CG.png" alt="Output 1" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/exp_without_bathymtery_CG/station_(9000%2C%203650)_timeseries_exp_without_bathymetry_CG.png">Figure 3: At station (9000, 3650)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/exp_without_bathymtery_CG/station_(11000%2C%203650)_timeseries_exp_without_bathymetry_CG.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/exp_without_bathymtery_CG/station_(11000%2C%203650)_timeseries_exp_without_bathymetry_CG.png" alt="Output 2" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/exp_without_bathymtery_CG/station_(11000%2C%203650)_timeseries_exp_without_bathymetry_CG.png">Figure 4: At station (11000, 3650)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/exp_without_bathymtery_CG/station_(13500%2C%203650)_timeseries_exp_without_bathymetry_CG.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/exp_without_bathymtery_CG/station_(13500%2C%203650)_timeseries_exp_without_bathymetry_CG.png" alt="Output 3" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/exp_without_bathymtery_CG/station_(13500%2C%203650)_timeseries_exp_without_bathymetry_CG.png">Figure 5: At station (13500, 3650)</a></b>
    </td>
    <td align="center">
      <video width="250" controls>
        <source src="https://drive.google.com/file/d/1sSBj288MSJ1-PuQOx2mpAKI8TERTMlPe/view?usp=sharing" type="video/mp4">
        Your browser does not support the video tag.
      </video>
      <br>
      <b>🎥 <a href="https://drive.google.com/file/d/1sSBj288MSJ1-PuQOx2mpAKI8TERTMlPe/view?usp=sharing">Video: Full Simulation</a></b>
    </td>
  </tr>
</table>



### ✅ Test Case Study 1 : SWEs on a 2D Sloped Beach (`Continuous Galerkin`)
`Directory:` FEniCS_Landlab/SWEs_2D_SlopedBeachProblem_SWEmniCS_CG

`Details:` This test case evaluates SWEMniCS' wetting and drying capability using Kärna’s α scheme with a `CG` approach on a `13,800 m × 7,200 m` sloping beach. The `initial` water surface is `flat`, and a `harmonic tide (2m amplitude, 12h period)` is applied at the `open left` boundary, while `other` boundaries are `walls`. The simulation runs for `7 days`, using Manning’s friction `0.02` s/m<sup>1/3</sup>, and records elevations and velocities at `x = 9000 m, 11,000 m, and 13,500 m.`

`Run the script:`

```sh
cd FEniCS_Landlab/SWEs_2D_SlopedBeachProblem_SWEmniCS_CG
python swe2d_sloped_beach_cg_main.py
```

<h2 align="center">🌊 Simulation Results</h2>

<table align="center">
  <tr>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_1_water_surface_elevation.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_1_water_surface_elevation.png" alt="Output 1" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_1_water_surface_elevation.png">Figure 6: At station (9000, 3650)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_1_x_velocity.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_1_x_velocity.png" alt="Output 2" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_1_x_velocity.png">Figure 7: At station (9000, 3650)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_1_y_velocity.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_1_y_velocity.png" alt="Output 3" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_1_y_velocity.png">Figure 8: At station (9000, 3650)</a></b>
    </td>
  </tr>
  <tr>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_2_water_surface_elevation.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_2_water_surface_elevation.png" alt="Output 4" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_2_water_surface_elevation.png">Figure 9: At station (11000, 3650)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_2_x_velocity.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_2_x_velocity.png" alt="Output 5" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_2_x_velocity.png">Figure 10: At station (11000, 3650)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_2_y_velocity.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_2_y_velocity.png" alt="Output 6" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_2_y_velocity.png">Figure 11: At station (11000, 3650)</a></b>
    </td>
  </tr>
  <tr>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_3_water_surface_elevation.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_3_water_surface_elevation.png" alt="Output 7" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_3_water_surface_elevation.png">Figure 12: At station (13500, 3650)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_3_x_velocity.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_3_x_velocity.png" alt="Output 8" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_3_x_velocity.png">Figure 13: At station (13500, 3650)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_3_y_velocity.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_3_y_velocity.png" alt="Output 9" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_CG/station_3_y_velocity.png">Figure 14: At station (13500, 3650)</a></b>
    </td>
  </tr>
</table>

<h3 align="center">🎥 Full Simulation Video</h3>
<p align="center">
  <a href="https://drive.google.com/file/d/1eGYL2J2Qrgfv_kozJB67OOYfqGQVYyKK/view?usp=sharing">
    <img src="https://drive.google.com/file/d/1eGYL2J2Qrgfv_kozJB67OOYfqGQVYyKK/view?usp=sharing" alt="Simulation Video" width="400">
  </a>
  <br>
  <b>🎥 <a href="https://drive.google.com/file/d/1eGYL2J2Qrgfv_kozJB67OOYfqGQVYyKK/view?usp=sharing">Video: Full Simulation</a></b>
</p>

### ✅ Test Case Study 2 : SWEs on a 2D Sloped Beach (`Discontinuous Galerkin`)
`Directory:` FEniCS_Landlab/SWEs_2D_SlopedBeachProblem_SWEmniCS_DG

`Details:` This test case evaluates SWEMniCS' wetting and drying capability using Kärna’s α scheme with a `DG` approach on a `13,800 m × 7,200 m` sloping beach. The `initial` water surface is `flat`, and a `harmonic tide (2m amplitude, 12h period)` is applied at the `open left` boundary, while `other` boundaries are `walls`. The simulation runs for `7 days`, using Manning’s friction `0.02` s/m<sup>1/3</sup>, and records elevations and velocities at `x = 9000 m, 11,000 m, and 13,500 m.`

`Run the script:`

```sh
cd FEniCS_Landlab/SWEs_2D_SlopedBeachProblem_SWEmniCS_DG
python swe2d_sloped_beach_dg_main.py
```

<h2 align="center">🌊 Simulation Results</h2>

<table align="center">
  <tr>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_1_water_surface_elevation.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_1_water_surface_elevation.png" alt="Output 1" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_1_water_surface_elevation.png">Figure 15: At station (9000, 3650)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_1_x_velocity.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_1_x_velocity.png" alt="Output 2" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_1_x_velocity.png">Figure 16: At station (9000, 3650)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_1_y_velocity.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_1_y_velocity.png" alt="Output 3" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_1_y_velocity.png">Figure 17: At station (9000, 3650)</a></b>
    </td>
  </tr>
  <tr>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_2_water_surface_elevation.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_2_water_surface_elevation.png" alt="Output 4" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_2_water_surface_elevation.png">Figure 18: At station (11000, 3650)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_2_x_velocity.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_2_x_velocity.png" alt="Output 5" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_2_x_velocity.png">Figure 19: At station (11000, 3650)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_2_y_velocity.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_2_y_velocity.png" alt="Output 6" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_2_y_velocity.png">Figure 20: At station (11000, 3650)</a></b>
    </td>
  </tr>
  <tr>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_3_water_surface_elevation.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_3_water_surface_elevation.png" alt="Output 7" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_3_water_surface_elevation.png">Figure 21: At station (13500, 3650)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_3_x_velocity.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_3_x_velocity.png" alt="Output 8" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_3_x_velocity.png">Figure 22: At station (13500, 3650)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_3_y_velocity.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_3_y_velocity.png" alt="Output 9" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/SWEs2D_slopeBeach_DG/station_3_y_velocity.png">Figure 23: At station (13500, 3650)</a></b>
    </td>
  </tr>
</table>

<h3 align="center">🎥 Full Simulation Video</h3>
<p align="center">
  <a href="https://drive.google.com/file/d/1kMRGNzSVy5nVHfbbiv38oP04piiggBKQ/view?usp=sharing">
    <img src="https://drive.google.com/file/d/1kMRGNzSVy5nVHfbbiv38oP04piiggBKQ/view?usp=sharing" alt="Simulation Video" width="400">
  </a>
  <br>
  <b>🎥 <a href="https://drive.google.com/file/d/1kMRGNzSVy5nVHfbbiv38oP04piiggBKQ/view?usp=sharing">Video: Full Simulation</a></b>
</p>


### ✅ Test 3 : SWEs on Real DEM Bathymetry
`Directory:` FEniCS_Landlab/Real_bathy_cg_exp

`Details:` This test evaluates the finite element implementation of the shallow water equations using real-world bathymetry data. The bathymetry is obtained from a DEM file `dem_swe.tif` downloaded from `OpenTopography`. The experiment follows the same parameters as Test Case Study 1, with continuous Galerkin (CG) elements.
The domain size from the dem data: `Lx = 7739.01 m, Ly = 14811.74 m`. The mesh resolution is set to: `nx = 6, ny = 12`

`Run the script:`

```sh
cd FEniCS_Landlab/Real_bathy_cg_exp
python main.py
```

<h2 align="center">🌊 Simulation Results</h2>

<table align="center">
  <tr>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_1_water_surface_elevation.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_1_water_surface_elevation.png" alt="Output 1" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_1_water_surface_elevation.png">Figure 24: At station (6450, 2470)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_1_x_velocity.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_1_x_velocity.png" alt="Output 2" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_1_x_velocity.png">Figure 25: At station (6450, 2470)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_1_y_velocity.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_1_y_velocity.png" alt="Output 3" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_1_y_velocity.png">Figure 26: At station (6450, 2470)</a></b>
    </td>
  </tr>
  <tr>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_2_water_surface_elevation.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_2_water_surface_elevation.png" alt="Output 4" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_2_water_surface_elevation.png">Figure 27: At station (6450, 6180)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_2_x_velocity.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_2_x_velocity.png" alt="Output 5" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_2_x_velocity.png">Figure 28: At station (6450, 6180)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_2_y_velocity.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_2_y_velocity.png" alt="Output 6" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_2_y_velocity.png">Figure 29: At station (6450, 6180)</a></b>
    </td>
  </tr>
  <tr>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_3_water_surface_elevation.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_3_water_surface_elevation.png" alt="Output 7" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_3_water_surface_elevation.png">Figure 30: At station (6450, 11110)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_3_x_velocity.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_3_x_velocity.png" alt="Output 8" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_3_x_velocity.png">Figure 31: At station (6450, 11110)</a></b>
    </td>
    <td align="center">
      <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_3_y_velocity.png">
        <img src="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_3_y_velocity.png" alt="Output 9" width="250">
      </a>
      <br>
      <b>🖼️ <a href="https://github.com/partha-sakha-paul/FEniCS_Landlab/blob/main/results/real_bathy_cg/station_3_y_velocity.png">Figure 32: At station (6450, 11110)</a></b>
    </td>
  </tr>
</table>

<h3 align="center">🎥 Full Simulation Video</h3>
<p align="center">
  <a href="https://drive.google.com/file/d/1kMRGNzSVy5nVHfbbiv38oP04piiggBKQ/view?usp=sharing">
    <img src="https://drive.google.com/file/d/1kMRGNzSVy5nVHfbbiv38oP04piiggBKQ/view?usp=sharing" alt="Simulation Video" width="400">
  </a>
  <br>
  <b>🎥 <a href="https://drive.google.com/file/d/1kMRGNzSVy5nVHfbbiv38oP04piiggBKQ/view?usp=sharing">Video: Full Simulation</a></b>
</p>
