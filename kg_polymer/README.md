Polymer Simulation Setup Tutorial
This tutorial guides you through setting up and running polymer simulations using the provided codebase.
Written by Arman Moussavi
Last updated: 10/16/2025

Prerequisites

Git installed on your system
Python installed (for running the build script)
LAMMPS installed (for running simulations)
Access to an HPC cluster (if running simulations on one)

Step-by-Step Instructions

Clone the RepositoryClone the repository to your desired destination directory:
git clone <repository_url> <destination_folder>


Navigate to the Initial Configuration FolderMove to the folder where the initial polymer structure will be generated:
cd <destination_folder>/write_initial_config


Generate the Initial Polymer StructureRun the Python script to build the initial configuration:
python build_melt.py


Navigate to the Generation FolderMove to the folder containing the simulation scripts:
cd ../generation


Run the Generation Script  

If running locally, execute the LAMMPS input script:lmp -i generation.inp


If running on an HPC cluster, modify the submission script (submit_generation.sub) as needed for your cluster's configuration, then submit the job:sbatch submit_generation.sub





Notes

Ensure all dependencies (e.g., LAMMPS, Python libraries) are installed before running the scripts.
For HPC cluster usage, check your cluster's documentation for specific submission script requirements.
Verify file paths in the scripts match your system setup.
