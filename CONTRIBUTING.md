# Contributing to FEniCS_Landlab

Thank you for considering contributing to FEniCS_Landlab! Your contributions help improve the simulation of 2D Shallow Water Equations using FEniCSx and DOLFINx on bathymetry data.

## How to Contribute

1. **Fork the Repository**: Click the "Fork" button at the top right of the repository page to create your own copy.

2. **Clone Your Fork**:
   ```bash
   git clone https://github.com/your-username/FEniCS_Landlab.git
   ```

3. **Set Up the Development Environment**:
   - **Install Docker Desktop**: Follow the [Docker Installation Guide](https://docs.docker.com/desktop/setup/install/windows-install/).
   - **Run DOLFINx in Docker**:
     - For the stable version:
       ```bash
       docker run -ti dolfinx/dolfinx:stable
       ```
     - For the nightly (latest) version:
       ```bash
       docker run -ti dolfinx/dolfinx:nightly
       ```
   - **Run Jupyter Notebook in Nightly Mode**:
     ```bash
     docker run --init -ti -p 8888:8888 dolfinx/lab:nightly
     ```
     Access the notebook at `http://localhost:8888` using the token provided in the terminal.

5. **Create a New Branch**:
   ```bash
   git checkout -b feature/your-feature-name
   ```

6. **Make Your Changes**:
   - Ensure the code adheres to the project's coding standards.
   - Update documentation as needed.

7. **Test Your Changes**:
   - Run existing tests to ensure your changes don't break the codebase.
   - Add new tests for your contributions.

8. **Commit Your Changes**:
   ```bash
   git add .
   git commit -m "Add feature: [brief description]"
   ```

9. **Push to Your Fork**:
   ```bash
   git push origin feature/your-feature-name
   ```

10. **Submit a Pull Request**:
    - Navigate to the original repository.
    - Click on "New Pull Request".
    - Select your branch and provide a clear description of your changes.

## Additional Resources

- [FEniCS Project](https://fenicsproject.org/)
- [Landlab Documentation](https://landlab.readthedocs.io/en/latest/tutorials/index.html)

We appreciate your contributions and look forward to collaborating with you!
