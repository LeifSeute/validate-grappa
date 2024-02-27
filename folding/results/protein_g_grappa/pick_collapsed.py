"""
PICK THE FRAME WITH THE SMALLEST MAXIMUM DISTANCES BETWEEN TWO ATOMS OF THE PROTEIN. THEN VALIDATE THAT THE STATE IS IN THE BOX AND SAVE IT AS A GRO FILE
"""
#%%
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import copy

# Load the trajectory and topology
u = mda.Universe('collapsed/md.gro', 'collapsed/md.trr')

protein = u.select_atoms("protein")

all_protein_positions = copy.deepcopy(np.array([
    protein.positions
    for _ in u.trajectory
]))

# calculate the maximum distance between any two atoms in the protein: (this is inefficient, upper triangular matrix etc would be better)
max_distances = np.max(np.max(np.linalg.norm(all_protein_positions[:, :, np.newaxis] - all_protein_positions[:, np.newaxis], axis=-1), axis=-1), axis=-1)

#%%
with mda.Writer("collapsed_structure.gro", protein.n_atoms) as W:
    W.write(protein)

print("Wrote collapsed_structure.gro")
#%%
print("VALIDATION: Check that the state is in the box")

# Plotting
plt.figure()
plt.hist(max_distances, bins=100)
plt.xlabel("Distance (Angstrom)")
plt.ylabel("Occurrences")
plt.title("Maxium distances for collapsed protein structure")
plt.savefig("max-distances.png")
plt.show()
# now plot the distances by frame:
plt.figure()
plt.scatter(np.arange(len(max_distances)), max_distances)
plt.xlabel("Frame")
plt.ylabel("Distance (Angstrom)")
plt.title("Max distances for collapsed protein structure")

plt.savefig("max-distances-by-frame.png")
plt.show()

min_frame = np.argmin(max_distances)

print(f"Smallest end-to-end distance: {max_distances[min_frame]:.2f} Angstrom in frame {min_frame}")

# save the frame with the smallest end-to-end distance as a gro file
ts = u.trajectory[min_frame]

#%%
    
u = mda.Universe('collapsed_structure.gro')

# claculate the end to end distance:
protein = u.select_atoms("protein")
positions = protein.positions
max_distance = np.max(np.max(np.linalg.norm(positions[:, np.newaxis] - positions, axis=-1), axis=-1))
assert np.isclose(max_distance, max_distances[min_frame], atol=0.1), f"Expected {max_distances[min_frame]} but got {max_distance}"

# scatter plot of 2 spatial dimensions:
plt.figure()
plt.scatter(protein.positions[:,0], protein.positions[:,1], color='red')
plt.xlabel("x (Angstrom)")
plt.ylabel("y (Angstrom)")
plt.title("Collapsed protein structure")
plt.show()

from mpl_toolkits.mplot3d import Axes3D
import copy

positions = copy.deepcopy(ts.positions)
protein_positions = copy.deepcopy(protein.positions)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(positions[:,0], positions[:,1], positions[:,2], color='blue', alpha=0.1)
ax.scatter(protein_positions[:,0], protein_positions[:,1], protein_positions[:,2], color='red')
# draw the bonds:
for i in range(len(protein_positions)-1):
    ax.plot(protein_positions[i:i+2,0], protein_positions[i:i+2,1], protein_positions[i:i+2,2], color='red')

ax.set_xlabel("x (nm)")
ax.set_ylabel("y (nm)")
ax.set_zlabel("z (nm)")
plt.title("Collapsed protein structure")
plt.show()

first_atom = protein_positions[0]
last_atom = protein_positions[-1]
print(f"End-to-end distance: {np.linalg.norm(last_atom - first_atom):.2f} Angstrom")


from pathlib import Path
import shutil

(Path(__file__).parent/'folding').mkdir(parents=True, exist_ok=True)

# copy the file in there as pep.gro
shutil.copy('collapsed_structure.gro', 'folding/pep.gro')

# copy the topology from collapsed/pep_grappa_initial
shutil.copy('collapsed/pep_grappa_initial.top', 'folding/pep.top')

