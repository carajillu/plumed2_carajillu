import mdtraj
import numpy as np
from scipy.spatial import cKDTree
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Process some floats.")
    parser.add_argument('--input', type=str, default="system.pdb", help='Protein-ligand input system')
    parser.add_argument('--ligresname', type=str, default=None, help='residue name of the ligand')
    parser.add_argument('--r_max',type=float, default=0.15, help="maximum distance between placed probe and closest ligand atom")
    parser.add_argument('--r_min',type=float, default=0.30, help="minimum distance between placed probes")
    parser.add_argument('--output', type=str, default="pseudo.pdb", help='PDB file of probes overlapping ligand')

    args = parser.parse_args()
    return args

def get_ligand(traj,lig_name):
    """
    Returns the coordinates of the specifide ligand from an MDTraj trajectory.
    
    Parameters:
    traj : mdtraj.Trajectory
        MDTraj trajectory object.
    lig_name:
        str residue name of the ligand
    
    Returns:
    heavy_atoms_traj : mdtraj.Trajectory
        MDTraj trajectory object containing only the ligand.
    """
    ligand_atom_indices = [atom.index for atom in traj.topology.atoms if atom.residue.name == lig_name]
    ligand_atoms_traj = traj.atom_slice(ligand_atom_indices)
    return ligand_atoms_traj

def get_heavy_atoms(traj):
    """
    Returns the coordinates of all heavy atoms (non-hydrogen) from an MDTraj trajectory.
    
    Parameters:
    traj : mdtraj.Trajectory
        MDTraj trajectory object.
    
    Returns:
    heavy_atoms_traj : mdtraj.Trajectory
        MDTraj trajectory object containing only the heavy atoms.
    """
    heavy_atom_indices = [atom.index for atom in traj.topology.atoms if atom.element.symbol != 'H']
    heavy_atoms_traj = traj.atom_slice(heavy_atom_indices)
    return heavy_atoms_traj


def generate_maximal_set_B(A, r_max, r_min, num_candidates=1000):
    """
    Generates a maximal set B of points given an array of coordinates A,
    ensuring each point in B is at most r_max from the closest point in A
    and at least r_min from any other point in B.
    
    Parameters:
    A : np.ndarray
        Array of shape (n, 3) containing the coordinates of points in A.
    r_max : float
        Maximum distance constraint for points in B from A.
    r_min : float
        Minimum distance constraint between points in B.
    num_candidates : int
        Number of candidate points to sample per atom.
    
    Returns:
    B_traj : mdtrajtraj.Trajectory
        MDTraj trajectory object containing the coordinates of points in B.
    """
    # Generate candidate points around A
    candidates = []
    for point in A:
        random_offsets = np.random.uniform(-1, 1, (num_candidates, 3))
        random_offsets /= np.linalg.norm(random_offsets, axis=1, keepdims=True)  # Normalize to unit vectors
        random_offsets *= np.random.uniform(0, r_max, (num_candidates, 1))  # Scale distances up to r_max
        new_points = point + random_offsets
        candidates.append(new_points)
    
    candidates = np.vstack(candidates)  # Shape (num_candidates * n, 3)
    
    # Filter candidates to ensure at most r_max distance to A
    tree_A = cKDTree(A)
    valid_mask = tree_A.query(candidates, distance_upper_bound=r_max)[0] != np.inf
    candidates = candidates[valid_mask]
    
    # Iteratively select points for B while ensuring r_min separation
    B = []
    tree_B = cKDTree(np.empty((0, 3)))  # Empty tree initially
    
    for candidate in candidates:
        if len(B) == 0 or tree_B.query(candidate, distance_upper_bound=r_min)[0] == np.inf:
            B.append(candidate)
            tree_B = cKDTree(np.array(B))  # Update tree with new point
    
    B = np.array(B)
    
    # Create MDTraj trajectory object for B
    topology = mdtraj.Topology()
    chain = topology.add_chain()
    for _ in range(B.shape[0]):
        residue = topology.add_residue("B", chain)
        topology.add_atom("B", mdtraj.element.carbon, residue)  # Arbitrary element choice
    
    B_traj = mdtraj.Trajectory(B[None, :, :], topology)  # Add batch dimension
    
    return B_traj

if __name__=="__main__":
   args=parse_args()
   z=mdtraj.load(args.input)
   if args.ligresname is not None:
      ligand=get_ligand(z,args.ligresname)
   else:
       ligand=z
   ligand_heavy=get_heavy_atoms(ligand)
   pseudo=generate_maximal_set_B(ligand_heavy.xyz[0],r_max=args.r_max,r_min=args.r_min)
   pseudo.save_pdb(args.output)