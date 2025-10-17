# ===============================================================
#  POLYMER ANALYSIS SCRIPT — EDUCATIONAL VERSION
#  Author: Arman Moussavi
#  Last updated: 2025-10-16
#
#  PURPOSE:
#    Reads a LAMMPS polymer data file and calculates:
#        - Radius of gyration (Rg)
#        - End-to-end distance (Ree)
#        - Asphericity
#        - Chain size
#    Then visualizes what each metric *means* physically.
# ===============================================================

import numpy as np
import re
import networkx as nx
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D

# ---------------------------------------------------------------
# 1. PARSING FUNCTION (same as before)
# ---------------------------------------------------------------
def parse_lammps_data(file_path):
    atoms, bonds, angles = [], [], []
    with open(file_path, 'r') as file:
        lines = file.readlines()

    atom_section = bond_section = angle_section = False
    box_bounds = {}
    box_regex = {
        'x': re.compile(r'^([-\d\.eE]+)\s+([-\d\.eE]+)\s+xlo\s+xhi'),
        'y': re.compile(r'^([-\d\.eE]+)\s+([-\d\.eE]+)\s+ylo\s+yhi'),
        'z': re.compile(r'^([-\d\.eE]+)\s+([-\d\.eE]+)\s+zlo\s+zhi')
    }

    for line in lines:
        line = line.strip()
        for axis in ['x', 'y', 'z']:
            match = box_regex[axis].match(line)
            if match:
                lo, hi = map(float, match.groups())
                box_bounds[axis] = hi - lo
                break

    for line in lines:
        stripped = line.strip()
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
            bond_section = atom_section = angle_section = False

        if atom_section and stripped and not stripped.startswith("#"):
            parts = stripped.split()
            atom_id = int(parts[0])
            mol_id = int(parts[1])
            atom_type = int(parts[2])
            x, y, z = map(float, parts[4:7])
            ix, iy, iz = map(int, parts[7:10])
            atoms.append((atom_id, mol_id, atom_type, (x, y, z), (ix, iy, iz)))

        if bond_section and stripped and not stripped.startswith("#"):
            parts = stripped.split()
            bond_id = int(parts[0])
            bond_type = int(parts[1])
            atom1, atom2 = int(parts[2]), int(parts[3])
            bonds.append((bond_id, bond_type, atom1, atom2))

        if angle_section and stripped and not stripped.startswith("#"):
            parts = stripped.split()
            angle_id = int(parts[0])
            angle_type = int(parts[1])
            atom1, atom2, atom3 = map(int, parts[2:5])
            angles.append((angle_id, angle_type, atom1, atom2, atom3))

    return atoms, bonds, angles, box_bounds

# ---------------------------------------------------------------
# 2. METRIC CALCULATION
# ---------------------------------------------------------------
def calculate_metrics(atoms, bonds, box_bounds):
    atom_coords = {atom[0]: np.array(atom[3]) for atom in atoms}
    Lx, Ly, Lz = box_bounds['x'], box_bounds['y'], box_bounds['z']
    atom_coords = {
        atom[0]: np.array(atom[3]) + np.array(atom[4]) * np.array([Lx, Ly, Lz])
        for atom in atoms
    }

    G = nx.Graph()
    for _, _, a1, a2 in bonds:
        G.add_edge(a1, a2)

    clusters = list(nx.connected_components(G))
    rg_results, ree_results, asph_results, sizes = [], [], [], []
    example_data = {}

    for i, cluster in enumerate(clusters, start=1):
        if len(cluster) < 2:
            continue

        coords = np.array([atom_coords[a] for a in cluster])
        com = np.mean(coords, axis=0)
        squared_distances = np.sum((coords - com) ** 2, axis=1)
        rg = np.sqrt(np.mean(squared_distances))
        rg_results.append(rg)
        sizes.append(len(cluster))

        subgraph = G.subgraph(cluster)
        ends = [n for n, d in subgraph.degree() if d == 1]
        if len(ends) == 2:
            ree = np.linalg.norm(atom_coords[ends[0]] - atom_coords[ends[1]])
            ree_results.append(ree)

            centered = coords - com
            S = np.matmul(centered.T, centered) / len(coords)
            eigvals = np.sort(np.linalg.eigvalsh(S))[::-1]
            trace = np.sum(eigvals)
            asph = 1 - 3 * (eigvals[0]*eigvals[1] + eigvals[0]*eigvals[2] + eigvals[1]*eigvals[2]) / trace**2 if trace != 0 else 0
            asph_results.append(asph)

            # Store an example chain for visualization
            if not example_data:
                example_data = {
                    'coords': coords,
                    'com': com,
                    'ends': [atom_coords[ends[0]], atom_coords[ends[1]]],
                    'rg': rg,
                    'ree': ree,
                    'asph': asph
                }

    metrics = {
        'Rg': {'mean': np.mean(rg_results), 'std': np.std(rg_results)},
        'Ree': {'mean': np.mean(ree_results), 'std': np.std(ree_results)},
        'Asphericity': {'mean': np.mean(asph_results), 'std': np.std(asph_results)},
        'Chain Size': {'mean': np.mean(sizes), 'std': np.std(sizes)},
    }

    return metrics, example_data, sizes

# ---------------------------------------------------------------
# 3. PLOTTING FUNCTIONS
# ---------------------------------------------------------------
def plot_metrics(example_data, sizes):
    coords = example_data['coords']
    com = example_data['com']
    end1, end2 = example_data['ends']

    fig = plt.figure(figsize=(14, 10))

    # (1) Radius of Gyration
    ax1 = fig.add_subplot(221, projection='3d')
    ax1.scatter(coords[:,0], coords[:,1], coords[:,2], c='b', label='Monomers')
    ax1.scatter(*com, c='r', s=100, label='Center of Mass')
    sphere = plt.Circle((0,0), example_data['rg'], fill=False)
    ax1.set_title(f"Radius of Gyration (Rg ≈ {example_data['rg']:.2f})")
    ax1.legend()

    # (2) End-to-End Distance
    ax2 = fig.add_subplot(222, projection='3d')
    ax2.plot(coords[:,0], coords[:,1], coords[:,2], 'o-', color='gray', alpha=0.7)
    ax2.plot([end1[0], end2[0]], [end1[1], end2[1]], [end1[2], end2[2]], 'r-', lw=3, label='End-to-End')
    ax2.set_title(f"End-to-End Distance (Ree ≈ {example_data['ree']:.2f})")
    ax2.legend()

    # (3) Asphericity (qualitative)
    ax3 = fig.add_subplot(223, projection='3d')
    centered = coords - com
    ax3.scatter(centered[:,0], centered[:,1], centered[:,2], c='g')
    ax3.set_title(f"Asphericity (A ≈ {example_data['asph']:.2f})")

    # (4) Chain Size Distribution
    ax4 = fig.add_subplot(224)
    ax4.hist(sizes, bins=10, color='purple', alpha=0.7)
    ax4.set_xlabel("Chain Size (# of atoms)")
    ax4.set_ylabel("Frequency")
    ax4.set_title("Chain Size Distribution")

    plt.tight_layout()
    plt.savefig("polymer_metrics.png")

# ---------------------------------------------------------------
# 4. RUN
# ---------------------------------------------------------------
file_path = "../generation/generation.data"
atoms, bonds, angles, box_bounds = parse_lammps_data(file_path)
metrics, example_data, sizes = calculate_metrics(atoms, bonds, box_bounds)

print("\nPolymer Conformation Metrics:")
print("{:<20} {:<10} {:<10}".format("Metric", "Mean", "Std"))
print("-" * 40)
for key, val in metrics.items():
    print("{:<20} {:<10.3f} {:<10.3f}".format(key, val['mean'], val['std']))

# Plot metric illustrations
if example_data:
    plot_metrics(example_data, sizes)
