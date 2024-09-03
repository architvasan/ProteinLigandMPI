"""This module contains functions docking a molecule to a receptor using Openeye.
The code is adapted from this repository: https://github.com/inspiremd/Model-generation
"""

from pose_generation.docking_openeye import *
import MDAnalysis as mda
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from functools import cache, partial
from pathlib import Path
from typing import List, Optional
import numpy as np
from openeye import oechem, oedocking, oeomega
import pandas as pd
from tqdm import tqdm
from docking_utils import smi_to_structure
from utils import exception_handler
import os
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from argparse import ArgumentParser, SUPPRESS

'''
Running Code
'''

### Parameters setting
parser = ArgumentParser()#add_help=False)

parser.add_argument(
    "-d", "--data", type=str, required=True, help="input data"
)

parser.add_argument(
    "-s", "--smicol", type=str, required=True, help="smiles_column for input data"
)

parser.add_argument(
    "-r", "--receptor", type=str, required=True, help="receptor file"
)

parser.add_argument(
    "-t", "--tempdir", type=str, required=True, help="temporary directory"
)

parser.add_argument(
    "-o", "--outtraj", type=str, required=True, help="output trajectory direcotyr"
)

parser.add_argument(
    "-S", "--scorepatt", type=str, required=True, help="score pattern"
)

parser.add_argument(
    "-P", "--protpdb", type=str, required=True, help="protein pdb"
)

args = parser.parse_args()
df_smiles = pd.read_csv(args.data)#'../passed_files/pains_pass.csv')#input/RTCB_TopHits_22B_100k.csv')[:1000] # smilesdatafile. smiles column is 'smiles'
recep_file = args.receptor#'input/rtcb-7p3b-receptor-5GP-A-DU-601.oedu' # receptor oedu file for openeye 
max_confs = 10 # confs to generate
score_cutoff = -10 # below this score generate pose
temp_dir = args.tempdir #'lig_confs' # store ligand poses temporarily
output_traj = args.outtraj #'output_combined_trajectories' # store pdb + dcd files for complex poses
score_pattern = args.scorepatt#'scores/rtcb_scores' #store scores like this (will have #ranks files)
protein_pdb = args.protpdb#'prot.pdb' #protein pdb file to use to store. Will save everything in this file to complex


### Running docking in parallel
run_list_docking(df_smiles, recep_file, max_confs, score_cutoff, protein_pdb, temp_dir, output_traj, score_pattern, 'temp')#args.data, args.receptor, max_confs, score_cutoff, args.protpdb, args.tempdir, args.outtraj, args.scorepatt, 'temp')#

