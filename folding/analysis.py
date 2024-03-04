# %%
from src.src_analysis import get_traj_paths, rmse_plot
import argparse

parser = argparse.ArgumentParser(description='')

parser.add_argument('--ref', type=str, default='ref.pdb', help='reference structure', required=True)
parser.add_argument('--trjdir', type=str, default='trj', help='directory with trajectories', required=True)
parser.add_argument('--ref_rmsd', type=float, default=None, help='reference rmsd in angstroem')
parser.add_argument('--figpath', type=str, default="rmsd_to_ref.png", help='path to save figure')

args = parser.parse_args()

traj_paths = get_traj_paths(args.trjdir)
rmse_plot(traj_paths, ref_path=args.ref, ref_rmsd=args.ref_rmsd, figpath=args.figpath)