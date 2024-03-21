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
parser.add_argument("--image", action="store_true", help="For image saving")
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
        align.alignto(u, ref, select="protein and name CA")

        rmsd_val = rmsd(u.select_atoms("protein and name CA").positions, ref.select_atoms("protein and name CA").positions, superposition=False)
        rmsd_values.append(rmsd_val)
        W.write(u.atoms)

starting_frame_path = str(Path(f"{NAME}_cluster_center.pdb"))
if not args.image or not Path(starting_frame_path).exists():
    starting_frame_path = str(trj_path.parent / "min_rmsd_frame.pdb")
    min_rmsd_index = rmsd_values.index(min(rmsd_values))
    u.trajectory[min_rmsd_index]
    align.alignto(u, ref, select="protein and name CA")
    u.atoms.write(starting_frame_path)
else:
    starting_universe = mda.Universe(starting_frame_path)
    align.alignto(starting_universe, ref, select="protein and name CA")

    starting_universe.atoms.write(starting_frame_path)

# VMD script content
vmd_script_content = f"""
# VMD script for loading and visualizing the aligned trajectory and the structure with lowest RMSD

# Load the reference structure
mol new {ref_centered_path} type pdb
mol modstyle 0 top NewCartoon
mol modcolor 0 top ColorID 0
"""

if args.image:
    vmd_script_content += f"""

# Load the molecule and the aligned trajectory
mol new {starting_frame_path} type pdb
mol modstyle 0 top NewCartoon
mol modcolor 0 top ColorID 1
"""
else:
    vmd_script_content += f"""
    
# Load the molecule and the aligned trajectory
mol new {starting_frame_path} type pdb
mol modstyle 0 top NewCartoon
mol modcolor 0 top ColorID 1
mol addfile {aligned_trj_path} waitfor all

# set the inital view to frame min_rmsd_index
animate goto {min_rmsd_index}
"""

visualization_settings = f"""
# Additional commands for display settings and camera position

# Set display to GLSL rendering, enable culling, set projection to Orthographic, and disable depth cue
display rendermode GLSL
display culling on
display projection Orthographic
display depthcue off

# Set display background to white and turn off axes
color Display Background white
axes location Off
"""
if args.image:
    vmd_script_content += visualization_settings

def reproduce_camera_setting(rotation_matrix, global_matrix, scale_matrix, center_matrix):
    return f"""
molinfo 0 set {{rotate_matrix}} {{{rotation_matrix}}}
molinfo 1 set {{rotate_matrix}} {{{rotation_matrix}}}
molinfo 0 set {{global_matrix}} {{{global_matrix}}}
molinfo 1 set {{global_matrix}} {{{global_matrix}}}
molinfo 0 set {{scale_matrix}} {{{scale_matrix}}}
molinfo 1 set {{scale_matrix}} {{{scale_matrix}}}
molinfo 0 set {{center_matrix}} {{{center_matrix}}}
molinfo 1 set {{center_matrix}} {{{center_matrix}}}
"""

if "chignolin" in NAME:
    ROTATE_MATRIX = "{{-0.107528 0.135811 -0.149437 0} {-0.201915 -0.074426 0.0776492 0} {-0.0025197 0.168389 0.154847 0} {0 0 0 1}}"
    GLOBAL_MATRIX = "{{1 0 0 0.07} {0 1 0 -0.19} {0 0 1 0} {0 0 0 1}}"
    SCALE_MATRIX = "{{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}"
    CENTER_MATRIX = "{{1 0 0 0.670732} {0 1 0 0.371891} {0 0 1 0.702051} {0 0 0 1}}"
    vmd_script_content += reproduce_camera_setting(ROTATE_MATRIX, GLOBAL_MATRIX, SCALE_MATRIX, CENTER_MATRIX)



# Save the VMD script to a file
vmd_script_path = "visualization_script.vmd"
with open(vmd_script_path, "w") as script_file:
    script_file.write(vmd_script_content)

print(f"VMD script saved to {vmd_script_path}")

os.system(f"vmd -e {vmd_script_path}")