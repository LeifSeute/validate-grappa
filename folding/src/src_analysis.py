# %%
import MDAnalysis as mda
from MDAnalysis.analysis import diffusionmap, align, rms, dihedrals
from MDAnalysis.analysis.rms import rmsd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
from pathlib import Path
from typing import List
import numpy as np
import copy
import tempfile
import os

# disable warnings
import warnings
warnings.filterwarnings('ignore')

# %%
def get_traj_paths(trjdir:str, key:str="folding"):
    trjdir = Path(trjdir)
    trjs = trjdir.glob(f"*{key}*/*md_centered.trr")
    trjs = list(trjs)
    print(f"found {len(trjs)} trajectories")
    return trjs

# %%
def rmse_plot(traj_paths:List[str], ref_path:str, ref_rmsd:float=None, figpath:str=None, amber:bool=False, clustering:bool=False):
    """
    Creates an rmsd over time plot for every trajectory in traj_paths. assumes that the same folder contains a pep.gro file.
    """

    if len(traj_paths) == 0:
        print("No trajectories given")
        return
    
    FONTSIZE = 16
    FONT = 'Arial'
    #default:
    # FONT = 'DejaVu Sans'
    plt.rc('font', family=FONT)
    plt.rc('xtick', labelsize=FONTSIZE)
    plt.rc('ytick', labelsize=FONTSIZE)
    plt.rc('axes', labelsize=FONTSIZE)
    plt.rc('legend', fontsize=FONTSIZE)

    fig, ax = plt.subplots(figsize=(6, 3))  # Adjusted figure size for better aspect ratio

    us_per_frame = 2e-3 # microseconds per frame (we write frames every 2ns)
    c_palette = [
        "#1f77b4", # blue
        "#4daf4a", # green
        "#6acc64",
        "#d65f5f",
        "#956cb4",
        "#8c613c",
        "#dc7ec0",
        "#797979",
        "#d5bb67",
        "#82c6e2",
    ]

    # if amber:
    #     # invert:
    #     c_palette = ['#e41a1c'] + c_palette[::-1]
    u_ref = mda.Universe(str(ref_path))
    ref_c_alphas = u_ref.select_atoms('protein and name CA')

    max_ts = 0
    max_trjs = 5
    y_lim = 8
    for i, trj in enumerate(traj_paths):
        if i > max_trjs - 1:
            break
        pep_file = Path(trj).parent / "pep.gro"
        u = mda.Universe(str(pep_file), str(trj))
        traj_c_alphas = u.select_atoms('protein and name CA')

        if len(ref_c_alphas) != len(traj_c_alphas):
            raise ValueError(f"Number of C alpha atoms does not match: Reference={len(ref_c_alphas)}, Trajectory {i}={len(traj_c_alphas)}")

        prealigner = align.AlignTraj(u, u_ref, select="protein and name CA", in_memory=True).run()
        val = rms.RMSD(u, u_ref, select='name CA').run()

        df = pd.DataFrame(val.results.rmsd, columns=['Frame', 'Time (ns)', 'C-alphas'])
        times_micr_s = np.array(df['Frame']) * us_per_frame
        c_alpha_rmse = np.array(df['C-alphas'])
        c_alpha_rmse = pd.Series(c_alpha_rmse).rolling(100).mean().values  # Smoothed RMSD

        ax.plot(times_micr_s, c_alpha_rmse, alpha=1, color=c_palette[i % len(c_palette)], linewidth=1.5)

        t_max = times_micr_s[-1]
        max_ts = max(max_ts, t_max)

        # max_ts = 10 # NOTE

    # Set MAJOR ticks on the y-axis less frequently, e.g., every 2 angstroms
    # Major ticks will have labels
    ax.set_yticks(range(0, int(y_lim) + 1, 2))

    # Set MINOR ticks on the y-axis at every angstrom
    # Minor ticks will not have labels by default, but you can use them for grid lines
    ax.set_yticks(range(0, int(y_lim) + 1, 1), minor=True)

    if ref_rmsd:
        ax.plot([0, max_ts], [ref_rmsd, ref_rmsd], label='Lindorff Larsen et al.', c='black', linestyle='--', alpha=0.8)
    else:

        # Enable the grid only for y-axis, customizing its appearance
        ax.yaxis.grid(True, which='both', color='black', linestyle=':', linewidth=0.8, alpha=0.3)


    ax.set_ylim(0, y_lim)  # Adjust y-axis limits
    ax.set_xlim(0, max_ts)  # Adjust x-axis limits
    ax.set_xlabel("Time [μs]")
    ax.set_ylabel("C alpha RMSD [Å]")
    plt.tight_layout()  # Adjust layout to make room for the x-label
    plt.legend(loc='upper right', frameon=False)
    figpath = figpath or "rmsd_to_ref.png"
    plt.savefig(figpath, dpi=400, bbox_inches='tight')  # Save with tight bounding box to include labels

    if clustering:
        print('clustering...')
        u_concat = mda.Universe(str(pep_file), *traj_paths)

        print(f"num_frames: {len(u_concat.trajectory)}")

        n_frames_max = 5000
        print(f"pruning the trajectory to {n_frames_max} frames")

        idx_map = np.random.choice(len(u_concat.trajectory), n_frames_max, replace=False)

        # reduce the trajectory to those frames:
        with tempfile.NamedTemporaryFile(delete=False, suffix='.xtc') as temp_traj_file:
            temp_traj_path = temp_traj_file.name  # Store the temporary file name

            # Create a Writer for the temporary trajectory, writing only the selected frames
            with mda.Writer(temp_traj_path, n_atoms=u_concat.atoms.n_atoms) as W:
                for frame_idx in idx_map:
                    u_concat.trajectory[frame_idx]  # Set the current frame
                    W.write(u_concat.atoms)  # Write the current frame to the temporary file

            # Now that the selected frames are written, create a new Universe with the temporary file
            new_u = mda.Universe(pep_file, temp_traj_path)

        print(f"num_frames: {len(new_u.trajectory)}")

        aligner = align.AlignTraj(new_u, new_u, select='name CA', in_memory=True).run()

        matrix = diffusionmap.DistanceMatrix(new_u, select='name CA').run()
        dist_matrix = copy.deepcopy(matrix.dist_matrix)

        cutoff = 1.0 #A
        n_clusters = 3
        for i in range(n_clusters):
            neighbors = np.sum(dist_matrix < cutoff,axis=0)
            center = np.argmax(neighbors)
            print(f"Frame # {center} is center {i} with {max(neighbors)} neighbors")
            
            new_u.trajectory[center]
            center_c_alpha_positions = new_u.select_atoms('name CA').positions

            # align to reference:
            rmsd_center_to_ref = rmsd(center_c_alpha_positions, ref_c_alphas.positions, superposition=True)
            print(f"Distance to reference: {rmsd_center_to_ref:.2f} Å")
            # remove rows and cols from dist_matrix:
            neighbor_idxs = np.where(dist_matrix[center] < cutoff)[0]
            dist_matrix = np.delete(dist_matrix, neighbor_idxs, axis=0)
            dist_matrix = np.delete(dist_matrix, neighbor_idxs, axis=1)


            if i == 0:
                # write the whole protein from the center state as pdb file:
                # align the C alphas of the protein to the reference and save it as pdb:
                protein = new_u.select_atoms("protein")


                figpath = figpath.replace('_ref_rmsd.png', f'_cluster_center.pdb')
                protein.write(figpath)

                # now load it again, align and overwrite:
                u_cluster_center = mda.Universe(figpath)

                # align the C alphas of the protein to the reference and save it as pdb:
                align.alignto(u_cluster_center, u_ref, select="protein and name CA")

                u_cluster_center.atoms.write(figpath)

        os.remove(temp_traj_path)

        print("\n\n")