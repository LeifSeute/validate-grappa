#%%

import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import os

PDBNAME = '1ubq'

FONT = 'Arial'
FONTSIZE = 18

GRAPPA_COLOR = '#1f77b4'
AMBER_COLOR = '#e41a1c'

MAX_TIME = 100 # ns

FIGSIZE = (9, 3.5)

N = 1000  # Number of points to sample
dt_ns = 2  # Duration in nanoseconds to calculate RMSD

grappa_dir = Path(__file__).parent/'results'/f'{PDBNAME}_grappa/mdrun'

amber_dir = grappa_dir.parent.parent/f'{PDBNAME}_amber99sbildn/mdrun'

pep_file_grappa = str(grappa_dir/'pep.gro')
trj_file_grappa = str(grappa_dir/'md_centered.trr')

pep_file_amber = str(amber_dir/'pep.gro')
trj_file_amber = str(amber_dir/'md_centered.trr')

experimental_structure = f'{PDBNAME}.pdb'


plt.rc('font', family=FONT)
plt.rc('xtick', labelsize=FONTSIZE)
plt.rc('ytick', labelsize=FONTSIZE)
plt.rc('axes', labelsize=FONTSIZE)
plt.rc('legend', fontsize=FONTSIZE)

#%%

# Load the universes for both simulations and the reference structure
u_grappa = mda.Universe(pep_file_grappa, trj_file_grappa)
u_amber = mda.Universe(pep_file_amber, trj_file_amber)



fig, ax = plt.subplots(2, 1, figsize=(FIGSIZE[0], 2*FIGSIZE[1]))

# Select the atoms you're interested in, e.g., all backbone atoms of the protein
# Adjust the selection as needed


protein_grappa = u_grappa.select_atoms('protein and name CA')

protein_amber = u_amber.select_atoms('protein and name CA')
time_step_ns = u_grappa.trajectory.dt / 1e3  # Convert to ns

# Initialize the RMSD analysis
# The reference is the first frame of the trajectory by default
recalc = False
if recalc:
    rmsd_analysis = rms.RMSD(protein_grappa, protein_grappa, ref_frame=0)

    rmsd_analysis_amber = rms.RMSD(protein_amber, protein_amber, ref_frame=0)

    # Run the analysis
    rmsd_analysis.run()

    rmsd_analysis_amber.run()

    # moving average over dt_avg ns
    dt_avg = 0.05  # ns
    n_avg = int(dt_avg / time_step_ns)
    rmsd_analysis.results.rmsd[:, 2] = np.convolve(rmsd_analysis.results.rmsd[:, 2], np.ones(n_avg)/n_avg, mode='same')
    rmsd_analysis_amber.results.rmsd[:, 2] = np.convolve(rmsd_analysis_amber.results.rmsd[:, 2], np.ones(n_avg)/n_avg, mode='same')

    # rmsd_analysis.results.rmsd contains the results
    # Column 0 is the frame number, column 1 is the time, and column 2 is the RMSD

    times_grappa = rmsd_analysis.results.rmsd[:, 1]/1e3
    rmsd_grappa = rmsd_analysis.results.rmsd[:, 2]

    times_amber = rmsd_analysis_amber.results.rmsd[:, 1]/1e3
    rmsd_amber = rmsd_analysis_amber.results.rmsd[:, 2]

times_grappa = times_grappa[times_grappa <= MAX_TIME]
rmsd_grappa = rmsd_grappa[times_grappa <= MAX_TIME]

times_amber = times_amber[times_amber <= MAX_TIME]
rmsd_amber = rmsd_amber[times_amber <= MAX_TIME]

# delete the non-convoluted frames:
times_grappa = times_grappa[:-n_avg]
rmsd_grappa = rmsd_grappa[:-n_avg]
times_amber = times_amber[:-n_avg]
rmsd_amber = rmsd_amber[:-n_avg]


# Plot the RMSD over time
ax[0].plot(times_amber, rmsd_amber, label='FF99SBILDN', color=AMBER_COLOR)
ax[0].plot(times_grappa, rmsd_grappa, label='Grappa', color=GRAPPA_COLOR)
ax[0].set_ylim(0, 5.5)
ax[0].set_xlim(0, MAX_TIME)
ax[0].set_xlabel('Time [ns]')
ax[0].set_ylabel('RMSD [Å]')
# plt.title('RMSD to the First Frame')

n_ticks_x = 4
n_ticks_y = 2

ax[0].xaxis.set_major_locator(plt.MaxNLocator(n_ticks_x))
ax[0].yaxis.set_major_locator(plt.MaxNLocator(n_ticks_y))

ax[0].legend(frameon=False, loc='upper left')
# invert the order of the legend:
handles, labels = ax[0].get_legend_handles_labels()
ax[0].legend(handles[::-1], labels[::-1], frameon=False, loc='upper left')

# fig.save(f'rmsd_to_first_frame_{PDBNAME}.png', dpi=400, bbox_inches='tight')


if not Path('rmsd_data.npy').exists():

    u_grappa = mda.Universe(pep_file_grappa, trj_file_grappa)
    u_amber = mda.Universe(pep_file_amber, trj_file_amber)

    # Select the atoms you're interested in, e.g., all backbone atoms of the protein
    # Adjust the selection as needed
    protein_grappa = u_grappa.select_atoms('protein and name CA')
    protein_amber = u_amber.select_atoms('protein and name CA')

    # Sampling parameters
    dt_ps = 1000 * dt_ns  # Convert ns to ps

    # Time step in ps, assuming constant time step throughout the trajectory
    time_step = u_grappa.trajectory.dt
    frames_per_dt = int(dt_ps / time_step)

    # Store all RMSD curves
    all_rmsds = []

    all_rmsds_amber = []

    min_frames = min(len(u_grappa.trajectory), len(u_amber.trajectory))

    print(f'Total number of frames: {min_frames}\n')

    # pick N points with replacement:
    idxs = np.random.choice(range(min_frames - frames_per_dt), N, replace=False)
    for i, start_frame in enumerate(idxs):
        print(f'Processing Interval {i + 1}/{len(idxs)}', end='\r')
        u_grappa.trajectory[start_frame]  # Set the current frame as reference
        ref_positions = protein_grappa.positions.copy()  # Reference positions
        rmsds = []

        u_amber.trajectory[start_frame]  # Set the current frame as reference
        ref_positions_amber = protein_amber.positions.copy()  # Reference positions
        rmsds_amber = []
        
        # Calculate RMSD for the next nanosecond
        for i in range(start_frame, start_frame + frames_per_dt):
            u_grappa.trajectory[i]
            rmsd = rms.rmsd(ref_positions, protein_grappa.positions, superposition=True)
            rmsds.append(rmsd)

            u_amber.trajectory[i]
            rmsd_amber = rms.rmsd(ref_positions_amber, protein_amber.positions, superposition=True)
            rmsds_amber.append(rmsd_amber)
        
        all_rmsds.append(rmsds)
        all_rmsds_amber.append(rmsds_amber)

    # Convert to numpy array for easier manipulation
    all_rmsds = np.array(all_rmsds)
    all_rmsds_amber = np.array(all_rmsds_amber)

    # Calculate mean, 25th, and 75th percentiles
    rmsd_mean = np.mean(all_rmsds, axis=0)
    rmsd_25th = np.percentile(all_rmsds, 25, axis=0)
    rmsd_75th = np.percentile(all_rmsds, 75, axis=0)

    rmsd_mean_amber = np.mean(all_rmsds_amber, axis=0)
    rmsd_25th_amber = np.percentile(all_rmsds_amber, 25, axis=0)
    rmsd_75th_amber = np.percentile(all_rmsds_amber, 75, axis=0)

    # Time axis for plotting
    time_axis = np.arange(0, dt_ps, time_step)[:len(rmsd_mean)]


    print('\nSaving data...')

    # save everything in a npy file:
    np.save('rmsd_data.npy', {
        'time_axis': time_axis,
        'rmsd_mean': rmsd_mean,
        'rmsd_25th': rmsd_25th,
        'rmsd_75th': rmsd_75th,
        'rmsd_mean_amber': rmsd_mean_amber,
        'rmsd_25th_amber': rmsd_25th_amber,
        'rmsd_75th_amber': rmsd_75th_amber,
    })

else:

    data = np.load('rmsd_data.npy', allow_pickle=True)
    time_axis = data.item().get('time_axis')
    rmsd_mean = data.item().get('rmsd_mean')
    rmsd_25th = data.item().get('rmsd_25th')
    rmsd_75th = data.item().get('rmsd_75th')
    rmsd_mean_amber = data.item().get('rmsd_mean_amber')
    rmsd_25th_amber = data.item().get('rmsd_25th_amber')
    rmsd_75th_amber = data.item().get('rmsd_75th_amber')

time_axis = time_axis / 1000  # convert to ns

# Plotting
ax[1].plot(time_axis, rmsd_mean_amber, label='FF99SBILDN', color=AMBER_COLOR, linewidth=2.5)
ax[1].fill_between(time_axis, rmsd_25th_amber, rmsd_75th_amber, alpha=0.35, color=AMBER_COLOR, linewidth=0.0)
ax[1].plot(time_axis, rmsd_mean, label='Grappa', color=GRAPPA_COLOR, linewidth=2.5)
ax[1].fill_between(time_axis, rmsd_25th, rmsd_75th, alpha=0.35, color=GRAPPA_COLOR, linewidth=0.0)
ax[1].set_xlabel(r'$\Delta t$ [ns]')
ax[1].set_ylabel(r'RMSD $(t + \Delta t)$ [Å]')
# plt.title('RMSD Over Time Interval')

ax[1].legend(frameon=False, loc='lower right')
# invert the order of the legend:
handles, labels = ax[0].get_legend_handles_labels()
ax[1].legend(handles[::-1], labels[::-1], frameon=False, loc='lower right')

# only use n ticks:
n_ticks = 4
n_yticks = 2

ax[1].xaxis.set_major_locator(plt.MaxNLocator(n_ticks))
ax[1].yaxis.set_major_locator(plt.MaxNLocator(n_yticks))

fig.tight_layout(pad=2)

plt.savefig(f'rmsd_{PDBNAME}.png', dpi=400, bbox_inches='tight')
# plt.show()
# plt.close()

# %%
