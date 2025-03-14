import json
import pypdb
import pandas as pd
import subprocess
import os
import argparse
import mdtraj
import numpy as np
import sys
import shutil

excluded_residues = [
    # WATER
    "WAT", "HOH", "SOL",

    #HEME GROUPS
    "HEM",  # Heme B (most common, hemoglobin, cytochromes)
    "HEC",  # Heme C (cytochrome C, covalently bound)
    "HEA",  # Heme A (cytochrome c oxidase)
    "HEM3", # Iron(III) protoporphyrin IX (oxidized heme)
    "HEO",  # Heme O (precursor to Heme A, bacterial cytochromes)
    "DHE",  # Deuterated Heme (experimental use)
    "HDD",  # Heme D (bacterial oxidases)
    "HAS"  # Heme AS (bacterial variant of Heme A)
]


def parse():
    parser = argparse.ArgumentParser(description="Process some floats.")
    parser.add_argument('--json', type=str, default="cryptobench.json", help='cryptobench dataset (json)')
    parser.add_argument("--r_max",type=float,default=0.6,help="maximum distance to consider that protein and ligand atoms are interacting.")
    args = parser.parse_args()
    return args

def json_to_pd(json_path: str):
    json_db=json.load(open(json_path,"r"))
    holo_pdb_id=[]
    holo_chain=[]
    ligand=[]
    ligand_index=[]
    ligand_chain=[]
    holo_pymol_selection=[]
    ligand_pymol_selection=[] # not present in the original json
    holo_mdtraj_selection=[]  # not present in the original json
    ligand_mdtraj_selection=[] # not present in the original json
    for key in json_db.keys():
        holo_pdb_id.append(json_db[key]["holo_pdb_id"])
        holo_chain.append(json_db[key]["holo_chain"])
        ligand.append(json_db[key]["ligand"])
        ligand_index.append(json_db[key]["ligand_index"])
        ligand_chain.append(json_db[key]["ligand_chain"])
        holo_pymol_selection.append(json_db[key]["holo_pymol_selection"])
        holo_mdtraj_selection.append(sel_pymol2mdtraj(json_db[key]["holo_pymol_selection"]))
        ligand_pymol_selection.append(f"{json_db[key]["holo_pdb_id"]} and chain {json_db[key]["ligand_chain"]} and resn {json_db[key]["ligand"]} and resi {json_db[key]["ligand_index"]}")
        ligand_mdtraj_selection.append(f"resname \'{json_db[key]["ligand"]}\' and resSeq {json_db[key]["ligand_index"]}")
    df=pd.DataFrame({
         "holo_pdb_id": holo_pdb_id,
         "holo_chain": holo_chain,
         "ligand": ligand,
         "ligand_index": ligand_index,
         "ligand_chain": ligand_chain,
         "holo_pymol_selection": holo_pymol_selection,
         "holo_mdtraj_selection":holo_mdtraj_selection,
         "ligand_pymol_selection":ligand_pymol_selection,
         "ligand_mdtraj_selection":ligand_mdtraj_selection
          })
    
    return df

def sel_pymol2mdtraj(pymol_selection_str: str):
    pymol_selection_str=pymol_selection_str.split()
    pymol_sel_lst=[]
    items_to_exclude=["","resi","and","chain"]
    for item in pymol_selection_str:
        item=item.strip("(").strip(")")
        if item in items_to_exclude:
            continue
        if "+" in item:
            item=item.split("+")
        pymol_sel_lst.append(item)
    #mdtraj_sel_str=f"chain {pymol_sel_lst[1]} and (resid {" or resid ".join(pymol_sel_lst[2])})" # chain id not available in mdtraj
    mdtraj_sel_str=f"resSeq {" or resSeq ".join(pymol_sel_lst[2])}"
    return mdtraj_sel_str

def pdb_get_chain(pdb_str: str, holo_chain: str, ligand_chain: str):
    #return pdb_str
    pdb_lst=pdb_str.split("\n")
    chain_lst=[]
    for line in pdb_lst:
        spline=line.split()
        if (len(spline)==0):
            continue
        if "ATOM" not in spline[0] and "HETATM" not in spline[0]:
            continue
        if line[21]!=holo_chain and line[21]!=ligand_chain:
            continue
        chain_lst.append(line)
    return("\n".join(chain_lst))

def filter_atoms_by_distance(A: mdtraj.Trajectory, B: mdtraj.Trajectory, frame: int=0, r_max: float=0.6):
    """
    Returns a new MDTraj object containing only the atoms in A that are within r_max distance
    from the closest atom in B.
    
    Parameters:
        A (md.Trajectory): The first MDTraj object (atoms to filter).
        B (md.Trajectory): The second MDTraj object (reference atoms).
        r_max (float): The maximum distance threshold.
    
    Returns:
        md.Trajectory: A new MDTraj object containing only the selected atoms from A.
    """
    # Compute pairwise distances between all atoms in A and all atoms in B
    xyz_A=A.xyz[frame]
    xyz_B=B.xyz[frame]
    distances=[]
    # Efficient computation of pairwise distances using broadcasting
    distances = np.linalg.norm(xyz_A[:, np.newaxis, :] - xyz_B[np.newaxis, :, :], axis=-1)

    
    # Reshape to get the minimum distance to any atom in B for each atom in A
    min_distances = distances.reshape(A.n_atoms, B.n_atoms).min(axis=1)
    
    # Select atoms in A that have at least one neighbor in B within r_max
    selected_indices = np.where(min_distances < r_max)[0]
    
    # Slice the trajectory to keep only the selected atoms
    filtered_traj = A.atom_slice(selected_indices)
    
    return filtered_traj

def build_benchmark(df:pd.DataFrame,r_max:float=0.6):
    rootdir=os.getcwd()
    n_success=[]
    n_fail_noligand=[]
    n_fail_noreceptor=[]
    n_fail_mdtraj=[]
    n_fail_download=[]
    for i in range(len(df.index)):
        
        holo_pdb_id=df.holo_pdb_id[i]

        #Ligand might have been in excluded gorups
        if df.ligand[i] in excluded_residues:
            print(f"************************* {holo_pdb_id} ligand not allowed: {df.ligand[i]} *************************")
            n_fail_noligand.append(holo_pdb_id)
            continue

        #get pdb structure
        pdb_str=pypdb.get_pdb_file(holo_pdb_id)
        if pdb_str is None: #if get_pdb fails, it returns None
           n_fail_download.append(holo_pdb_id)
           print(f"************************* failed to download file {holo_pdb_id}.pdb  *************************")
           if os.path.isfile(f"{holo_pdb_id}.pdb"):
                shutil.rmtree(f"{holo_pdb_id}.pdb")
           continue

        with open(f"{holo_pdb_id}.pdb","w") as f:
             f.write(pdb_str)

        try:
            structure_obj=mdtraj.load(f"{holo_pdb_id}.pdb")
        except:
            n_fail_mdtraj.append(holo_pdb_id)
            print(f"************************* {holo_pdb_id} cannot be loaded with mdtraj *************************")
            continue
        
        #filter unwanted things
        selection=[]
        for atom in structure_obj.topology.atoms:
            if atom.residue.name not in excluded_residues and atom.element!=mdtraj.element.hydrogen:
                selection.append(atom.index)
        structure_obj=structure_obj.atom_slice(selection)
        

        #Get get ligand and receptor objects
        rec_sel=[]
        lig_sel=structure_obj.topology.select(df.ligand_mdtraj_selection[i])
        for atom in structure_obj.topology.atoms:
            if atom.index not in lig_sel:
               rec_sel.append(atom.index)

        
        ligand=structure_obj.atom_slice(lig_sel)
        receptor=structure_obj.atom_slice(rec_sel)
        
        
        #Filter atoms that are too far
        ligand_slice=filter_atoms_by_distance(ligand,receptor,r_max=r_max)
        receptor_slice=filter_atoms_by_distance(receptor,ligand_slice,r_max=r_max)
        if (ligand_slice.n_atoms==0) or (receptor_slice.n_atoms==0):
            n_fail_noreceptor.append(holo_pdb_id)
            print(f"************************* ligand in {holo_pdb_id} seems too far from the protein to be correct *************************")
            continue

        #export structures
        rootdir=os.getcwd()
        os.makedirs(holo_pdb_id,exist_ok=True)
        os.chdir(holo_pdb_id)
        ligand_slice.save_pdb("ligand.pdb")
        receptor_slice.save_gro("pocket.gro")
        shutil.move(f"{rootdir}/{holo_pdb_id}.pdb",f"{rootdir}/{holo_pdb_id}/{holo_pdb_id}.pdb")
        os.chdir(rootdir)
        n_success.append(holo_pdb_id)
        print(f"************************* successfully exported structures for {holo_pdb_id} *************************")

    print(f"Number of successfully processed structres: {len(n_success)}")

    print(f"Number of structures that failed due to illegal ligand: {len(n_fail_noligand)}")
    if len(n_fail_noligand)>0:
       print(n_fail_noligand)
       os.makedirs("illegal_ligands",exist_ok=True)
       for item in n_fail_noligand:
           shutil.move(f"{rootdir}/{item}.pdb",f"{rootdir}/illegal_ligands/{item}.pdb")

    print(f"Number of structures that failed due to wrong ligand (far from receoptor): {len(n_fail_noreceptor)}")
    if len(n_fail_noreceptor)>0:
       print(n_fail_noreceptor)
       os.makedirs("wrong_ligands",exist_ok=True)
       for item in n_fail_noreceptor:
           shutil.move(f"{rootdir}/{item}.pdb",f"{rootdir}/illegal_ligands/{item}.pdb")

    print(f"Number of structures that couuld not be loaded on mdtraj: {len(n_fail_mdtraj)}")
    if len(n_fail_mdtraj)>0:
       print(n_fail_mdtraj)
       os.makedirs("mdtraj_error",exist_ok=True)
       for item in n_fail_mdtraj:
           shutil.move(f"{rootdir}/{item}.pdb",f"{rootdir}/illegal_ligands/{item}.pdb")

    print(f"Number of structures that could not be downloaded: {len(n_fail_download)}")
    if len(n_fail_download)>0:
       print(n_fail_download)
       os.makedirs("download_error",exist_ok=True)
       for item in n_fail_download:
           shutil.move(f"{rootdir}/{item}.pdb",f"{rootdir}/illegal_ligands/{item}.pdb")

    return



if __name__=="__main__":
    args=parse()
    df=json_to_pd(args.json)
    print(df)
    df.to_csv("cryptobench_main.csv", index=False, sep=" ")
    build_benchmark(df,r_max=args.r_max)