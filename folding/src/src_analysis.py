# %%
import MDAnalysis as mda
from MDAnalysis.analysis import diffusionmap, align, rms, dihedrals
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
from pathlib import Path
from typing import List
import numpy as np

# %%
def get_traj_paths(trjdir:str, key:str="folding"):
    trjdir = Path(trjdir)
    trjs = trjdir.glob(f"*{key}*/*md_centered.trr")
    trjs = list(trjs)
    print(f"found {len(trjs)} trajectories")
    return trjs

# %%
def rmse_plot(traj_paths:List[str], ref_path:str, ff_name:str=None, ref_rmsd:float=None, figpath:str=None):
    """
    Creates an rmsd over time plot for every trajectory in traj_paths. assumes that the same folder contains a pep.gro file.
    """

    if len(traj_paths) == 0:
        print("No trajectories given")
        return
    
    FONTSIZE = 16
    FONT = 'Arial'
    plt.rc('font', family=FONT)
    plt.rc('xtick', labelsize=FONTSIZE)
    plt.rc('ytick', labelsize=FONTSIZE)
    plt.rc('axes', labelsize=FONTSIZE)
    plt.rc('legend', fontsize=FONTSIZE)

    fig, ax = plt.subplots(figsize=(8, 4))  # Adjusted figure size for better aspect ratio

    us_per_frame = 2e-3 # microseconds per frame (we write frames every 2ns)
    c_palette = [
        "#4878d0",
        "#ee854a",
        "#6acc64",
        "#d65f5f",
        "#956cb4",
        "#8c613c",
        "#dc7ec0",
        "#797979",
        "#d5bb67",
        "#82c6e2",
    ]

    u_ref = mda.Universe(str(ref_path))
    ref_c_alphas = u_ref.select_atoms('protein and name CA')

    max_ts = 0
    for i, trj in enumerate(traj_paths):
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

        ax.plot(times_micr_s, c_alpha_rmse, alpha=0.8, color=c_palette[i % len(c_palette)])

        t_max = times_micr_s[-1]
        max_ts = max(max_ts, t_max)

    if ref_rmsd:
        ax.plot([0, max_ts], [ref_rmsd, ref_rmsd], label='Lindorff Larsen et al.', c='r', linestyle='--')

    ax.set_ylim(0, 6)  # Adjust y-axis limits
    ax.set_xlabel("Time [μs]")
    ax.set_ylabel("C alpha RMSD [Å]")
    plt.tight_layout()  # Adjust layout to make room for the x-label
    plt.legend(loc='upper right', frameon=False)
    figpath = figpath or "rmsd_to_ref.png"
    plt.savefig(figpath, dpi=400, bbox_inches='tight')  # Save with tight bounding box to include labels
