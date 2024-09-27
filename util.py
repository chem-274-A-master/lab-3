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