from src_pick_collapsed import pick_smallest
import argparse

parser = argparse.ArgumentParser(description='PICK THE FRAME WITH THE SMALLEST MAXIMUM DISTANCES BETWEEN TWO ATOMS OF THE PROTEIN. THEN VALIDATE THAT THE STATE IS IN THE BOX AND SAVE IT AS A GRO FILE')
parser.add_argument('--gro_path', '-gro', type=str, required=True, help='path to the gro file describing the protein')
parser.add_argument('--trr_path', '-trr', type=str, required=True, help='path to the trr file describing the trajectory')
parser.add_argument('--collapsed_gro_path', '-o', type=str, required=True, help='path to the output gro file')
parser.add_argument('--figdir', '-fig', type=str, default=None, help='path to the output figures')

args = parser.parse_args()

pick_smallest(args.gro_path, args.trr_path, args.collapsed_gro_path, args.figdir)