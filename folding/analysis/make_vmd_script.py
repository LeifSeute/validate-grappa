import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
from pathlib import Path
import os

import argparse

if not __name__ == "__main__":
    raise ImportError("This script should be run as a standalone script.")

# Parse the arguments
parser = argparse.ArgumentParser(description="Create a VMD script for visualizing the aligned trajectory and the structure with the lowest RMSD.")
parser.add_argument("name", type=str, help="The name of the system. E.g. chignolin_grappa")
args = parser.parse_args()

NAME = args.name

name = '_'.join(NAME.split('_')[:-1])

# Paths
pdb_path = Path.cwd().parent / "results" / NAME / "equilibration" / "pep.pdb"
trj_paths = list((Path.cwd().parent / "results" / NAME / "mds").iterdir())
trj_path = trj_paths[0] / "md_centered.trr" if trj_paths else None
ref_path = Path.cwd().parent / "references" / f"{name}.pdb"
ref_centered_path = Path.cwd().parent / "references" / f"{name}_centered.pdb"
aligned_trj_path = str(trj_path.parent / "aligned_md.trr")
min_rmsd_pdb_path = str(trj_path.parent / "min_rmsd_frame.pdb")

# Load the structures
u = mda.Universe(pdb_path, trj_path)
ref = mda.Universe(ref_path)

# Center the reference structure
ref.atoms.translate(-ref.atoms.center_of_mass())

# write the ref:
ref.atoms.write(str(ref_centered_path))

# Align the trajectory to the reference structure and find the frame with minimal RMSD
rmsd_values = []
with mda.Writer(aligned_trj_path, u.atoms.n_atoms) as W:
    for ts in u.trajectory:
        align.alignto(u, ref, select="protein and name CA", weights="mass")

        rmsd_val = rmsd(u.select_atoms("protein and name CA").positions, ref.select_atoms("protein and name CA").positions, superposition=False)
        rmsd_values.append(rmsd_val)
        W.write(u.atoms)

min_rmsd_index = rmsd_values.index(min(rmsd_values))
u.trajectory[min_rmsd_index]
align.alignto(u, ref, select="protein and name CA", weights="mass")
u.atoms.write(min_rmsd_pdb_path)

# VMD script content
vmd_script_content = f"""
# VMD script for loading and visualizing the aligned trajectory and the structure with lowest RMSD

# Load the reference structure
mol new {ref_centered_path} type pdb
mol modstyle 0 top NewCartoon
mol modcolor 0 top ColorID 0

# Load the molecule and the aligned trajectory
mol new {min_rmsd_pdb_path} type pdb
mol modstyle 0 top NewCartoon
mol modcolor 0 top ColorID 1
mol addfile {aligned_trj_path} waitfor all

# set the inital view to frame min_rmsd_index
animate goto {min_rmsd_index}
"""

# Save the VMD script to a file
vmd_script_path = "visualization_script.vmd"
with open(vmd_script_path, "w") as script_file:
    script_file.write(vmd_script_content)

print(f"VMD script saved to {vmd_script_path}")

os.system(f"vmd -e {vmd_script_path}")