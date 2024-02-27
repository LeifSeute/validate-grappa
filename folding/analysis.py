# %%
from src.src_analysis import get_traj_paths, rmse_plot
import argparse

parser = argparse.ArgumentParser(description='')

parser.add_argument('--ref', type=str, default='ref.pdb', help='reference structure', required=True)
parser.add_argument('--trjdir', type=str, default='trj', help='directory with trajectories', required=True)
parser.add_argument('--key', type=str, default='grappa', help='key for the forcefield', required=True)
parser.add_argument('--skip', type=int, default=1000000, help='stride for the trajectory')
parser.add_argument('--ff', type=str, default='ff', help='forcefield name')

args = parser.parse_args()

traj_paths = get_traj_paths(args.trjdir, args.key)
rmse_plot(traj_paths, args.ref, args.skip, args.key, args.ff)