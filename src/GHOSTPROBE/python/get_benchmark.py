import json
import pypdb
import pandas as pd
import subprocess
import os
import argparse
import mdtraj

def parse():
    parser = argparse.ArgumentParser(description="Process some floats.")
    parser.add_argument('--json', type=str, default="cryptobench.json", help='cryptobench dataset (json)')
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


def build_benchmark(df:pd.DataFrame):
    rootdir=os.getcwd()
    for i in range(len(df.index)):
        holo_pdb_id=df.holo_pdb_id[i]
        holo_chain=df.holo_chain[i]
        ligand_chain=df.ligand_chain[i]
        os.makedirs(holo_pdb_id,exist_ok=True)
        os.chdir(holo_pdb_id)
        filtered_pdb=pdb_get_chain(pypdb.get_pdb_file(holo_pdb_id),holo_chain,ligand_chain)

        with open(f"{holo_pdb_id}.pdb","w") as f:
            f.write(filtered_pdb)

        structure=mdtraj.load(f"{holo_pdb_id}.pdb")
        
        pocket_sel=structure.topology.select(df.holo_mdtraj_selection[i])
        pocket_obj=structure.atom_slice(pocket_sel)
        pocket_obj.save_gro("pocket.gro") # so that we can use make_ndx

        ligand_sel=structure.topology.select(df.ligand_mdtraj_selection[i])
        print(df.ligand_mdtraj_selection[i])
        print(ligand_sel)
        ligand_obj=structure.atom_slice(ligand_sel)
        ligand_obj.save_pdb("ligand.pdb")
        
        os.chdir(rootdir)
    
    return



if __name__=="__main__":
    args=parse()
    df=json_to_pd(args.json)
    print(df)
    build_benchmark(df)