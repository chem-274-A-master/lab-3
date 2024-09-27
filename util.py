import MDAnalysis as mda
import numpy as np
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import networkx as nx
import matplotlib.pyplot as plt

def build_supercell(u, frame_number=-1):

    # Move to the last frame of the trajectory
    u.trajectory[-1]
    
    # Get the atom positions and box dimensions
    positions = u.atoms.positions
    box = u.dimensions[:3]  # Only box lengths
    
    # Initialize an empty list to hold new replicated atoms and positions
    new_positions = []
    new_atom_types = []
    
    # Replicate the system (2x2x2) by shifting the coordinates
    for x_shift in range(2):
        for y_shift in range(2):
            for z_shift in range(2):
                shift_vector = np.array([x_shift, y_shift, z_shift]) * box
                shifted_positions = positions + shift_vector
                new_positions.append(shifted_positions)
    
    # Stack all the replicated coordinates
    replicated_positions = np.vstack(new_positions)
    
    # Create a new Universe with the replicated positions and the original topology
    # Use mda.Merge to create a new system with the same atom topology replicated across the new positions
    # We'll replicate the atom group multiple times in a new merged universe
    replicated_universe = mda.Merge(*[u.atoms for _ in range(8)])  # 2x2x2 = 8 total replicas
    replicated_universe.atoms.positions = replicated_positions
    
    return replicated_universe
    


def measure_cluster(universe, selection='resname OCT', min_distance=6.0):
    """
    Function to cluster octane molecules in the given MDAnalysis Universe
    based on a minimum distance cutoff, and return cluster sizes over time.

    Parameters:
    universe : MDAnalysis.Universe
        The MDAnalysis Universe object for the trajectory.
    selection : str
        Atom selection string to select octane molecules (default 'resname OCT').
    min_distance : float
        Distance cutoff (in Ångstroms) to consider molecules as part of the same cluster.

    Returns:
    clusters_per_frame : list of lists
        A list where each entry corresponds to the clusters found in that frame.
    cluster_sizes_per_frame : list of lists
        A list of lists where each entry corresponds to the sizes of clusters in each frame.
    """
    # Select octane molecules based on the selection string
    octane = universe.select_atoms(selection)

    # Number of octane molecules (residues)
    n_octane = len(octane.residues)

    clusters_per_frame = []
    cluster_sizes_per_frame = []

    # Loop over trajectory frames
    for ts in universe.trajectory:
        mySet = set()

        # Calculate distances between octane molecules
        coms = octane.residues.center_of_mass(compound='residues')  # Get COM for all residues at once

        for i in range(n_octane - 1):
            for j in range(i + 1, n_octane):
                # Calculate the minimum distance between COMs
                dist = np.linalg.norm(coms[i] - coms[j])

                # If within the cutoff distance, add to the set of connected components
                if dist <= min_distance:
                    mySet.add((i, j))

        # Create a graph to represent clusters
        graph = nx.from_edgelist(mySet)

        # Track clusters and their sizes
        frame_clusters = []
        frame_cluster_sizes = []
        for i in range(n_octane):
            if i not in graph:
                frame_clusters.append([i])
                frame_cluster_sizes.append(1)
            else:
                cluster = list(nx.node_connected_component(graph, i))
                frame_clusters.append(cluster)
                frame_cluster_sizes.append(len(cluster))

        clusters_per_frame.append(frame_clusters)
        cluster_sizes_per_frame.append(frame_cluster_sizes)

    return clusters_per_frame, cluster_sizes_per_frame



# Example usage: Plot the largest cluster size over time
def plot_largest_cluster_size(universe, selection='resname OCT', min_distance=6.0):
    # Get the clusters and cluster sizes per frame
    _, cluster_sizes_per_frame = cluster_octane(universe, selection, min_distance)

    # Extract the largest cluster size in each frame
    largest_cluster_sizes = [max(sizes) for sizes in cluster_sizes_per_frame]

    # Plot the largest cluster size over time
    plt.figure(figsize=(8, 6))
    plt.plot(largest_cluster_sizes, label="Largest Cluster Size")
    plt.xlabel("Frame")
    plt.ylabel("Cluster Size")
    plt.title("Largest Cluster Size Over Time")
    plt.legend()
    plt.show()

