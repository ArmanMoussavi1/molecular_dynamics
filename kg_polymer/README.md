# Polymer Simulation Setup Tutorial
Written by Arman Moussavi
Last updated: 10/16/2025

This tutorial guides you through setting up, running, and analyzing polymer simulations using the provided codebase.

## Prerequisites
- Git installed on your system
- Python installed (for running the build and analysis scripts)
- LAMMPS installed (for running simulations)
- Access to an HPC cluster (if running simulations on one)

## Step-by-Step Instructions

1. **Clone the Repository**  
   Clone the repository to your desired destination directory:
   ```bash
   git clone <repository_url> <destination_folder>
   ```

2. **Navigate to the Initial Configuration Folder**  
   Move to the folder where the initial polymer structure will be generated:
   ```bash
   cd <destination_folder>/write_initial_config
   ```

3. **Generate the Initial Polymer Structure**  
   Run the Python script to build the initial configuration:
   ```bash
   python build_melt.py
   ```

4. **Navigate to the Generation Folder**  
   Move to the folder containing the simulation scripts:
   ```bash
   cd ../generation
   ```

5. **Run the Generation Script**  
   - If running locally, execute the LAMMPS input script:
     ```bash
     lmp -i generation.inp
     ```
   - If running on an HPC cluster, modify the submission script (`submit_generation.sub`) as needed for your cluster's configuration, then submit the job:
     ```bash
     sbatch submit_generation.sub
     ```

6. **Analyze the Final Polymer Structure**  
   Navigate to the analysis folder and run the Python script to analyze the final structure of the polymer system:
   ```bash
   cd ../analysis
   python conformations.py
   ```

## Notes
- Ensure all dependencies (e.g., LAMMPS, Python libraries) are installed before running the scripts.
- For HPC cluster usage, check your cluster's documentation for specific submission script requirements.
- Verify file paths in the scripts match your system setup.
