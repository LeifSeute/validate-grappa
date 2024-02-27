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



#%%
def get_traj_paths(trjdir:str, key:str):
    trjdir = Path(trjdir)
    trjs = trjdir.glob(f"*{key}*/*_protein.xtc")
    trjs = list(trjs)
    print(f"found {len(trjs)} trajectories")

#%%
def rmse_plot(traj_paths:List[str], ref_path:str, skip=1000000, ff_name:str='ff'):

    us_per_frame = 1/5000 # microseconds per frame
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
    max_ts = 0
    for i,trj in enumerate(traj_paths):
        print(i)
        u = mda.Universe(str(ref_path),str(trj))

        prealigner = align.AlignTraj(u, u_ref,select="protein and name CA", in_memory=True).run()
        val = rms.RMSD(u,u_ref,select='name CA').run()

        df = pd.DataFrame(val.rmsd, columns=['Frame',
                                    'Time (ns)',
                                    'C-alphas'])
        df['Time (μs)'] = df['Frame'].apply(lambda x: x *us_per_frame)
        df['C-alphas'] = df['C-alphas'].rolling(100).mean()

        plt.plot(df['Time (μs)'],df['C-alphas'],alpha=0.8,color=c_palette[i])
        t_max = len(df)*us_per_frame
        max_ts = t_max if t_max > max_ts else max_ts
    # plt.plot([0,max_ts],[1.0,1.0],label='ref RMSD',c='r')
    plt.plot([0,18],[1.0,1.0],label='ref RMSD',c='r')
    plt.scatter(18,1.6,c='k',s=40,label='Target')
    x_sin = np.arange(0,18,0.01)
    plt.plot(x_sin,9*np.cos(x_sin/13),c='k',linestyle='-.',label='Trend')
    plt.xlabel("Time [μs]")
    plt.yticks([0,1,2,3,4,5])
    plt.ylabel("RMSD [Å]")
    plt.title(f"RMSD to experimental structure with {ff_name}")
    plt.legend()
    plt.savefig(f"rmsd_to_ref_{ff_name.replace(' ','_').replace('/','_').replace('.','_')}.png",dpi=500)



### end of RMSD plot
#%%
# print(u)
# protein = u.select_atoms("protein")
# with mda.Writer(f"protein_{key}_{skip}.xtc", protein.n_atoms) as W:
#     for ts in u.trajectory[::skip]:
#         W.write(u)



# #%%
# # Ramachandran plot
# # Compute phi and psi angles
# r = u.select_atoms("protein")
# R = dihedrals.Ramachandran(r).run()

# #%%
# # Create a 2D histogram for the phi and psi angles
# plt.hist2d(R.angles[:,:,0].flatten(), R.angles[:,:,1].flatten(), bins=100, cmap='cividis',norm='log')
# # Label the axes with larger font size
# plt.xlabel('Phi [degrees]', fontsize=14)
# plt.ylabel('Psi [degrees]', fontsize=14)
# # Add title with larger font size
# plt.title(f"Ramachandran Plot with {key} ff", fontsize=16)
# plt.colorbar(label='# of snapshots')
# plt.savefig(f"Ramachandran_plot_{key}.png",dpi=300)
# # Display the plot
# plt.show()


# #%%
# u = mda.Universe(str(ref),f"protein_{key}_{skip}.xtc")
# print(u,len(u.trajectory))

# # %%
# aligner = align.AlignTraj(u, u, select='name CA',
#                           in_memory=True).run()
# # %%
# matrix = diffusionmap.DistanceMatrix(u, select='name CA').run()
# matrix.dist_matrix.shape

# # %%
# cutoff = 2 #A
# neighbors = np.sum(matrix.dist_matrix < cutoff,axis=0)
# center = np.argmax(neighbors)

# print(f"Frame # {center} is center 1 with {max(neighbors)} neighbors")


# # %%
# prealigner = align.AlignTraj(u, u_ref,select="protein and name CA", in_memory=True).run()
# u.trajectory[center]
# protein = u.select_atoms("protein")
# protein.write(f"center_{key}.pdb")

# #%%
# u.trajectory[center]
# protein = u.select_atoms("name CA")
# pos = protein.positions.copy()
# ref = u_ref.select_atoms("name CA")
# pos_ref = ref.positions.copy()
# print("RMSD between reference structure and chosen center:")
# print(rms.rmsd(pos,pos_ref))

# # %%
# plt.figure()
# plt.imshow(matrix.dist_matrix, cmap='cividis')
# plt.xlabel('Frame')
# plt.ylabel('Frame')
# plt.colorbar(label=r'RMSD ($\AA$)')
# plt.title(f"Pairwise RMSD of 10 simulations with {key} ff")
# plt.savefig(f"pairwise_rmsd_{key}_stride{skip}.png",dpi=300)

# %%
