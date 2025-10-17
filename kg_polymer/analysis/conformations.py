# ===============================================================
#  POLYMER ANALYSIS SCRIPT — EDUCATIONAL VERSION
#  Author: Arman Moussavi
#  Last updated: 2025-04-01
#
#  PURPOSE:
#    This script reads a LAMMPS polymer data file and calculates 
#    structural metrics such as:
#        - Radius of gyration (Rg)
#        - End-to-end distance (Ree)
#        - Asphericity (shape anisotropy)
#        - Chain size (number of atoms per polymer)
#
#  AUDIENCE:
#    Beginners in polymer physics and molecular simulation.
#    Comments explain both the physical meaning and code logic.
# ===============================================================

# --------------------------
# Import necessary libraries
# --------------------------
import numpy as np                  # For numerical operations (arrays, math)
import re                           # For pattern matching (reading file structure)
import networkx as nx               # For graph analysis (polymer connectivity)
import matplotlib.pyplot as plt     # (Optional) For visualization
import os                           # For file handling

# ===============================================================
# 1. PARSING A LAMMPS DATA FILE
# ===============================================================
# This function reads a LAMMPS ".data" file that contains information about
# atoms, bonds, and angles — all needed to reconstruct the polymer topology.
#
# PHYSICS BACKGROUND:
# - "Atoms" represent monomers or beads.
# - "Bonds" connect pairs of atoms (like covalent bonds).
# - "Angles" connect triplets of atoms (used for stiffness or bending).
#
# GOAL:
# Extract this information into Python lists for later analysis.
# ===============================================================

def parse_lammps_data(file_path):
    # Prepare empty lists to store atom, bond, and angle information
    atoms = []
    bonds = []
    angles = []

    # Read the entire file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Initialize flags to track which section we’re currently reading
    atom_section = False
    bond_section = False
    angle_section = False

    # Dictionary to hold box dimensions (simulation box size)
    box_bounds = {}

    # Regular expressions to detect box limits in x, y, and z directions
    box_regex = {
        'x': re.compile(r'^([-\d\.eE]+)\s+([-\d\.eE]+)\s+xlo\s+xhi'),
        'y': re.compile(r'^([-\d\.eE]+)\s+([-\d\.eE]+)\s+ylo\s+yhi'),
        'z': re.compile(r'^([-\d\.eE]+)\s+([-\d\.eE]+)\s+zlo\s+zhi')
    }

    # ----------------------------------------------
    # STEP 1: Extract the box size from the header
    # ----------------------------------------------
    for line in lines:
        line = line.strip()
        for axis in ['x', 'y', 'z']:
            match = box_regex[axis].match(line)
            if match:
                lo, hi = map(float, match.groups())
                box_bounds[axis] = hi - lo
                break

    # ----------------------------------------------
    # STEP 2: Extract atoms, bonds, and angles
    # ----------------------------------------------
    for line in lines:
        stripped = line.strip()

        # Detect when each section starts
        if stripped == "Atoms # full":
            atom_section = True
            continue
        elif stripped == "Bonds":
            bond_section, atom_section = True, False
            continue
        elif stripped == "Angles":
            angle_section, bond_section, atom_section = True, False, False
            continue
        elif stripped == "Velocities":
            # Stop reading structural data when velocities begin
            bond_section = atom_section = angle_section = False

        # --------------------------
        # ATOMS SECTION
        # --------------------------
        if atom_section and stripped and not stripped.startswith("#"):
            parts = stripped.split()
            atom_id = int(parts[0])
            mol_id = int(parts[1])
            atom_type = int(parts[2])
            x, y, z = map(float, parts[4:7])          # positions
            ix, iy, iz = map(int, parts[7:10])        # image flags (periodicity)
            atoms.append((atom_id, mol_id, atom_type, (x, y, z), (ix, iy, iz)))

        # --------------------------
        # BONDS SECTION
        # --------------------------
        if bond_section and stripped and not stripped.startswith("#"):
            parts = stripped.split()
            bond_id = int(parts[0])
            bond_type = int(parts[1])
            atom1, atom2 = int(parts[2]), int(parts[3])
            bonds.append((bond_id, bond_type, atom1, atom2))

        # --------------------------
        # ANGLES SECTION
        # --------------------------
        if angle_section and stripped and not stripped.startswith("#"):
            parts = stripped.split()
            angle_id = int(parts[0])
            angle_type = int(parts[1])
            atom1, atom2, atom3 = map(int, parts[2:5])
            angles.append((angle_id, angle_type, atom1, atom2, atom3))

    # Return all structural data and box dimensions
    return atoms, bonds, angles, box_bounds


# ===============================================================
# 2. CALCULATING POLYMER STRUCTURAL METRICS
# ===============================================================
# PHYSICS IDEA:
# - A polymer can be treated as a chain (graph) of connected monomers.
# - We can calculate shape-related quantities to understand its conformation.
#
# QUANTITIES CALCULATED:
#   • Rg (radius of gyration): spread of the polymer around its center of mass.
#   • Ree (end-to-end distance): distance between the two chain ends.
#   • Asphericity: how elongated or spherical the shape is.
#   • Chain size: number of monomers in each chain.
# ===============================================================

def calculate_metrics(atoms, bonds, box_bounds):
    # ----------------------------------------------
    # STEP 1: Map atom IDs → coordinates
    # ----------------------------------------------
    atom_coords = {atom[0]: np.array(atom[3]) for atom in atoms}

    # ----------------------------------------------
    # STEP 2: Unwrap periodic boundary conditions
    # ----------------------------------------------
    # LAMMPS "image flags" tell how many times an atom crosses the boundary.
    # Multiply them by the box length to get unwrapped (true) coordinates.
    Lx, Ly, Lz = box_bounds['x'], box_bounds['y'], box_bounds['z']
    atom_coords = {
        atom[0]: np.array(atom[3]) + np.array(atom[4]) * np.array([Lx, Ly, Lz])
        for atom in atoms
    }

    # ----------------------------------------------
    # STEP 3: Build the polymer network as a graph
    # ----------------------------------------------
    # Each atom = node; each bond = edge connecting two nodes.
    G = nx.Graph()
    for _, _, atom1, atom2 in bonds:
        G.add_edge(atom1, atom2)

    # Identify all connected components (individual chains or clusters)
    clusters = list(nx.connected_components(G))

    # Storage for results
    rg_results, ree_results, asph_results, sizes = [], [], [], []

    # ----------------------------------------------
    # STEP 4: Loop over each polymer chain
    # ----------------------------------------------
    for i, cluster in enumerate(clusters, start=1):
        if len(cluster) < 2:
            print(f"Cluster {i}: only one atom, skipping.")
            continue

        # Get coordinates of all atoms in this chain
        coords = np.array([atom_coords[a] for a in cluster])

        # (a) Compute center of mass
        com = np.mean(coords, axis=0)

        # (b) Radius of gyration: average spread around COM
        squared_distances = np.sum((coords - com) ** 2, axis=1)
        rg = np.sqrt(np.mean(squared_distances))
        rg_results.append(rg)
        sizes.append(len(cluster))

        # Find chain ends (atoms with degree = 1)
        subgraph = G.subgraph(cluster)
        ends = [node for node, deg in subgraph.degree() if deg == 1]

        # (c) End-to-end distance only defined for linear chains
        if len(ends) == 2:
            ree = np.linalg.norm(atom_coords[ends[0]] - atom_coords[ends[1]])
            ree_results.append(ree)

            # (d) Asphericity from eigenvalues of gyration tensor
            centered = coords - com
            S = np.matmul(centered.T, centered) / len(coords)
            eigvals = np.sort(np.linalg.eigvalsh(S))[::-1]
            trace = np.sum(eigvals)
            if trace == 0:
                asph = 0.0
            else:
                products = (eigvals[0]*eigvals[1] + eigvals[0]*eigvals[2] + eigvals[1]*eigvals[2])
                asph = 1 - 3 * products / trace**2
            asph_results.append(asph)
        else:
            print(f"Cluster {i}: not linear (has {len(ends)} ends), skipping Ree/Asphericity.")

    # ----------------------------------------------
    # STEP 5: Summarize results for all chains
    # ----------------------------------------------
    if rg_results:
        metrics = {
            'Radius of Gyration (Rg)': {
                'mean': np.mean(rg_results),
                'std': np.std(rg_results)
            },
            'End-to-End Distance (Ree)': {
                'mean': np.mean(ree_results) if ree_results else np.nan,
                'std': np.std(ree_results) if ree_results else np.nan
            },
            'Asphericity': {
                'mean': np.mean(asph_results) if asph_results else np.nan,
                'std': np.std(asph_results) if asph_results else np.nan
            },
            'Chain Size': {
                'mean': np.mean(sizes),
                'std': np.std(sizes)
            },
        }
    else:
        metrics = {}
        print("No multi-atom clusters found. Check your data file.")

    return metrics


# ===============================================================
# 3. RUNNING THE SCRIPT
# ===============================================================
# You can replace 'generation.data' with your own LAMMPS data file.
# ===============================================================

file_path = "../generation/generation.data"

atoms, bonds, angles, box_bounds = parse_lammps_data(file_path)
metrics = calculate_metrics(atoms, bonds, box_bounds)

# ---------------------------------------------------------------
# Print a summary table of polymer structural metrics
# ---------------------------------------------------------------
if metrics:
    print("\nPolymer Conformation Metrics:")
    print("{:<25} {:<12} {:<12}".format("Metric", "Mean", "Std"))
    print("-" * 50)
    for key, val in metrics.items():
        mean_val = val['mean'] if not np.isnan(val['mean']) else "N/A"
        std_val = val['std'] if not np.isnan(val['std']) else "N/A"
        print("{:<25} {:<12} {:<12}".format(key, mean_val, std_val))
else:
    print("\nNo valid polymer data found.")
