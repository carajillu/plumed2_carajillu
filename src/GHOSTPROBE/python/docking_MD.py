"""
This program would take as input a dry MD trajectory of the protein (-f), a MD trajectory of the corresponding Ghost Probes (-p), 
relevant probe-$i-stats.csv (-s) and a benzene.pdbqt file (could expand for other ligands?).

1) parse arguments

2) Load protein and probes trajectory on mdtraj

3) 

For each probe at each snapshot, we will dock a benzene in a small docking box centered on the probe.


The output would be a consolidated probe-stats.csv file which will also include the best docking score
"""

import argparse
import mdtraj



def parse():
    return

def mk_prepare_receptor(protein_frame: mdtraj.Trajectory):
    return

def prepare_vina_config(probe_frame: mdtraj.Trajectory, probe_id: int):
    return

def run_vina_docking(vina_config_file: str):
    return

if __name__=="__main__":
    print("Hello world")
