# %%
# %%
import os 
import sys
from pathlib import Path
import scipy.io as sio
import numpy as np
import networkx as nx
import glob 
import argparse
import matplotlib.pyplot as plt

# %%
## import custom modules 
repoDir = Path().resolve().parent.parent ## dir of the repo root
print("repoDir:", repoDir)
configDir = repoDir / "config"
sys.path.append(str(configDir))

## define path variables 
from settings import DATA_PATH, CODE_PATH
sys.path.append(str(CODE_PATH["2-photon"]))

from utils.util_funcs import extract_value, strcmp, hdf5ToDict
from ephys.behavior import Behavior
from ephys.cell import Cell
from ephys.cell_ensemble import CellEnsemble
from ephys.events import Events
from ephys.events_process import EventsProcess
from ephys.session import Session

# %%
## define the data path for the session to analyze 
parser = argparse.ArgumentParser(description="Process neural and pupil data for a given session.")
parser.add_argument("--session_path", type=str, required=True, help="Path to the session data directory.")
args = parser.parse_args()
DATA_PATH_AROUSAL = Path(args.session_path)
# print("DATA_PATH:", DATA_PATH)

## make sure that data path exists
assert DATA_PATH_AROUSAL.exists(), f"Data path {DATA_PATH_AROUSAL} does not exist!"
# print(os.listdir(DATA_PATH_AROUSAL))

# %%
## load the .mat files 
behaviorFilePath = glob.glob(os.path.join(DATA_PATH_AROUSAL, "Beh*.mat"))[0]
behaviorFile = Behavior(behaviorFilePath)

headFile = behaviorFile.get_head_data() 
reyeFile = behaviorFile.get_reye_data()

## define the ephys dir 
sessionName = reyeFile["dirName"]
sessionName = extract_value(sessionName, 'dirName')
ephysDirPath = os.path.join(DATA_PATH_AROUSAL, sessionName)

## load the sortedUnits file
sortedUnitsFilePath = glob.glob(os.path.join(ephysDirPath, "Sort*.mat"))[0]
sortedUnitsFile = sio.loadmat(sortedUnitsFilePath, struct_as_record=True)
sortedUnitsFile = sortedUnitsFile['SortedUnits']
sortedUnitsFile = sortedUnitsFile.reshape(-1)

# %%
## create a list of Cell objects to pass to CellEnsemble class 
cells = []
for cell in sortedUnitsFile:
    cells.append(Cell(cell))

## create a CellEnsemble object
cellEnsemble = CellEnsemble(cells)

## create the Events object and make a Session object
eventsFilePath = glob.glob(os.path.join(ephysDirPath, "Events*.mat"))[0]
eventsFile = sio.loadmat(eventsFilePath, struct_as_record=True)["Events"]
eventsFile = eventsFile.reshape(-1)

events = [] 
for event in eventsFile:
    events.append(Events(event))

session = Session(events)
sessionDF = session.toDataFrame()

# %%
pupilDiameter = behaviorFile.get_pupilDiameter_data()
cellEnsembleDF = cellEnsemble.toDataFrame()

# %%
## Extract stimulus onset and offset times 
ms = 1e3 ## convert time in seconds to milliseconds
timeWindow = [-0.5, 0,5] ## time window in sec around stimulus onset to extract spike times

# %%
## define the frame numbers when the first and last stimuli were presented 
sessionFrameNums = reyeFile["TTs"][0, 0].reshape(-1) ## frame numbers when the stimuli were presented
firstEventFrame = sessionFrameNums[0] ## first frame number when the first stimulus was presented
lastEventFrame = sessionFrameNums[-1] ## last frame number when the last stimulus was presented

# %%
## extract the time for OpenEphys for the first and last stimulus presentation
OEFirstEventTime = session.get_event(0).soundcard_trigger_timestamp_sec ## OE time in seconds when the first stimulus was presented
OELastEventTime = session.get_event(-1).soundcard_trigger_timestamp_sec ## OE time in seconds when the last stimulus was presented 

## Calculate frame rate and create time alignment
print(f"First event frame: {firstEventFrame}, Last event frame: {lastEventFrame}")
print(f"First event time: {OEFirstEventTime:.3f}s, Last event time: {OELastEventTime:.3f}s")

# %%
## Calculate frames per second (frame rate)
framesPerSec = (lastEventFrame - firstEventFrame) / (OELastEventTime - OEFirstEventTime)
print(f"Calculated frame rate: {framesPerSec:.2f} frames/sec")
frameRate = extract_value(reyeFile['vid'][0][0]['framerate']) 
print(f"Actual frame rate: {frameRate}")

# %%
## Create frame number array for pupil data (assuming pupil data starts from frame 0)
pupil_frames = np.arange(len(pupilDiameter))  ## frame numbers for each pupil measurement

## Convert frame numbers to time (in sec) -> aligns frame numbers to the start of OEStart time
intercept = firstEventFrame - (framesPerSec * OEFirstEventTime)
timePupilEphys = (pupil_frames - intercept) / framesPerSec

## Convert to milliseconds if needed
timePupilEphys_ms = timePupilEphys * ms

# %%
## Create aligned pupil data structure
pupil_aligned = {
    'diameter': pupilDiameter,
    'frames': pupil_frames,
    'time_sec': timePupilEphys,
    'time_ms': timePupilEphys_ms,
    'frame_rate': framesPerSec
}

# print(f"Pupil data aligned: {len(pupilDiameter)} samples")
# print(f"Time range: {timePupilEphys[0]:.3f}s to {timePupilEphys[-1]:.3f} s")
# print(f"Duration: {timePupilEphys[-1] - timePupilEphys[0]:.3f} s")

# %%
## Extract spike data aligned to pupil timing
# Get spike times for all cells
all_spike_times = []
cell_ids = []
spike_counts_binned = []

for i, cell in enumerate(cells):
    spike_times = cell.spiketimes  # Assuming spike times are in seconds
    all_spike_times.extend([(t, i) for t in spike_times])  # (time, cell_id) pairs
    cell_ids.append(i)
    # print(f"Cell {i}: {len(spike_times)} spikes")

# %%
# Sort spikes by time
all_spike_times.sort(key=lambda x: x[0])
spike_times_only = [st[0] for st in all_spike_times]
spike_cell_ids = [st[1] for st in all_spike_times]

# print(f"Total spikes across all cells: {len(all_spike_times)}")
# print(f"Spike time range: {min(spike_times_only):.3f}s to {max(spike_times_only):.3f}s")

# %% [markdown]
# Binning the spikes as per the pupil sampling rate.

# %%
# Create binned spike data to match pupil sampling rate
pupil_sampling_interval = np.median(np.diff(timePupilEphys))  # Median interval between pupil samples
print(f"Pupil sampling interval: {pupil_sampling_interval*1000:.1f} ms")

# %%
# Bin spikes into pupil time bins
pupil_time_bins = timePupilEphys ## binning the spikes in these bins 
spike_rate_per_cell = np.zeros((len(cells), len(pupil_time_bins)))

for cell_idx in range(len(cells)):
    cell_spikes = [st[0] for st in all_spike_times if st[1] == cell_idx] ## incliding the spikes for the current cell[cell_idx], this is actually inefficient 
    if len(cell_spikes) > 0:
        # Count spikes in each pupil time bin
        spike_counts, _ = np.histogram(cell_spikes, bins=np.concatenate([pupil_time_bins, [pupil_time_bins[-1] + pupil_sampling_interval]])) ## adding the last bin limit
        spike_rate_per_cell[cell_idx, :] = spike_counts / pupil_sampling_interval  # Convert to firing rate (Hz) so divide by time unit of sec

# %%
# Calculate population firing rate per time bin
population_firing_rate = np.mean(spike_rate_per_cell, axis=0) ## per bin 
total_spike_count = np.sum(spike_rate_per_cell, axis=0) ## per bin 

# print(f"Population firing rate range: {np.min(population_firing_rate):.2f} to {np.max(population_firing_rate):.2f} Hz")

# %%
## Create comprehensive aligned data structure
neural_pupil_data = {
    'pupil': {
        'diameter': pupilDiameter,
        'time_sec': timePupilEphys,
        'time_ms': timePupilEphys_ms,
        'frames': pupil_frames
    },
    'spikes': {
        'all_spike_times': spike_times_only,
        'all_spike_cells': spike_cell_ids,
        'spike_rate_per_cell': spike_rate_per_cell,
        'population_firing_rate': population_firing_rate,
        'total_spike_count': total_spike_count,
        'time_bins': pupil_time_bins
    },
    'stimuli': {
        'event_times': [session.get_event(i).soundcard_trigger_timestamp_sec for i in range(len(sessionFrameNums))],
        'event_frames': sessionFrameNums
    },
    'frame_rate': framesPerSec
}

# %%
## Validation: Check alignment at event times of pupil frames
# print("\n=== VALIDATION: Pupil alignment at stimulus times ===")
event_times = neural_pupil_data["stimuli"]["event_times"]
event_frames = sessionFrameNums

# %% [markdown]
# Adding GC network and aligning it temporally. 

# %%
temp_network_path = glob.glob(str(Path(DATA_PATH_AROUSAL, "AGC_out*.mat")))[0]
assert(Path(temp_network_path)).exists()

# %%
AGC_data = hdf5ToDict(temp_network_path, include=["AGC"])

# %%
temporal_net = AGC_data["AGC"]

# %%
# print(nbins, temporal_net.shape)

# %%
## Examine the network data structure
# print("Network matrix shape:", temporal_net.shape)
# print("Number of bins:", nbins.item() if hasattr(nbins, 'item') else nbins)

# The temporal_net should be (n_cells, n_cells, n_time_bins)
n_cells_network = temporal_net.shape[0]
n_time_bins = temporal_net.shape[2] if len(temporal_net.shape) > 2 else 1 ## checking if there are multiple network mats

# print(f"Network cells: {n_cells_network}, Time bins: {n_time_bins}")
# print(f"Spike data cells: {len(cells)}")

# Check if network cells match spike data cells
if n_cells_network != len(cells):
    print(f"WARNING: Network has {n_cells_network} cells but spike data has {len(cells)} cells")
else:
    print("Network and spike data have matching cell counts")

# %%
# Calculate network time bin duration
temporal_network_time_start = 0 # in seconds -> this is derived from a separate algorithm ran separately
temporal_network_time_end = 1800 # in seconds
total_network_duration = temporal_network_time_end - temporal_network_time_start  # in seconds
n_time_bins = temporal_net.shape[2]  # number of time bins in the network data

# Each network time bin duration
network_bin_duration = total_network_duration / n_time_bins  # seconds per bin
network_bin_duration_ms = network_bin_duration * 1000  # milliseconds per bin

# print(f"Network temporal resolution:")
# print(f"  Total duration: {total_network_duration:.1f} seconds")
# print(f"  Network bins: {n_time_bins}")
# print(f"  Bin duration: {network_bin_duration_ms:.2f} ms")
# print(f"  Bin duration: {network_bin_duration:.4f} seconds")

# Compare with pupil sampling
pupil_sampling_interval_ms = pupil_sampling_interval * 1000
# print(f"\nComparison:")
# print(f"  Pupil sampling interval: {pupil_sampling_interval_ms:.2f} ms")
# print(f"  Network bin duration: {network_bin_duration_ms:.2f} ms")
# print(f"  Ratio (Network/Pupil): {network_bin_duration/pupil_sampling_interval:.2f}x")

# %%
## Create temporal alignment for network data
total_recording_duration = timePupilEphys[-1] - timePupilEphys[0]  # in seconds
network_time_bins = n_time_bins

network_time_array = np.linspace(temporal_network_time_start, temporal_network_time_end, network_time_bins) ## giving each network a time value for alignment

# print(f"Network temporal alignment:")
# print(f"  Total duration: {total_recording_duration:.1f} seconds")
# print(f"  Network bins: {network_time_bins}")
# print(f"  Bin duration: {network_bin_duration*1000:.2f} ms")
# print(f"  Network time range: {network_time_array[0]:.3f}s to {network_time_array[-1]:.3f}s")

# Compare with pupil sampling
# print(f"  Pupil samples: {len(timePupilEphys)}")
# print(f"  Pupil sampling interval: {pupil_sampling_interval*1000:.2f} ms")
# print(f"  Network/Pupil ratio: {network_time_bins/len(timePupilEphys):.1f}x")

# %%
## Calculate network measures over time
from scipy.stats import zscore

# Calculate various network measures for each time bin
network_measures = {}

# Check for any issues with network data
# print(f"Network data check:")
# print(f"  Contains NaN: {np.any(np.isnan(temporal_net))}")
# print(f"  Contains Inf: {np.any(np.isinf(temporal_net))}")
# print(f"  Data range: {np.nanmin(temporal_net):.6f} to {np.nanmax(temporal_net):.6f}")

# %%
## convert NaN to zeros 
temporal_net = np.nan_to_num(temporal_net, nan=0.0)

# recheck for any issues with network data
# print(f"Network data check:")
# print(f"  Contains NaN: {np.any(np.isnan(temporal_net))}")
# print(f"  Contains Inf: {np.any(np.isinf(temporal_net))}")
# print(f"  Data range: {np.nanmin(temporal_net):.6f} to {np.nanmax(temporal_net):.6f}")

# %% [markdown]
# Not creating the graph object just yet.

# %%
## basic network measures, for our graph 1 & 2 are the same (even 3)
# 1. Global network strength (sum of all connections)
network_strength = np.zeros(n_time_bins)
for t in range(n_time_bins):
    adj_matrix = temporal_net[:, :, t]
    network_strength[t] = np.nansum(adj_matrix)

# 2. Average connection strength (mean of all connections)
avg_connection_strength = np.zeros(n_time_bins)
for t in range(n_time_bins):
    adj_matrix = temporal_net[:, :, t]
    # Excluding diagonal (self-connections)
    off_diagonal = adj_matrix[~np.eye(adj_matrix.shape[0], dtype=bool)]
    avg_connection_strength[t] = np.nanmean(off_diagonal)

# 3. Network density (proportion of non-zero connections)
network_density = np.zeros(n_time_bins)
for t in range(n_time_bins):
    adj_matrix = temporal_net[:, :, t]
    off_diagonal = adj_matrix[~np.eye(adj_matrix.shape[0], dtype=bool)]
    # Count non-NaN, non-zero connections
    valid_connections = ~np.isnan(off_diagonal) & (off_diagonal != 0)
    total_possible = len(off_diagonal)
    network_density[t] = np.sum(valid_connections) / total_possible if total_possible > 0 else 0

# %%
# 4. Hub cell identification (cells with high out-degree)
outgoing_hub_strength = np.zeros((n_cells_network, n_time_bins))

for t in range(n_time_bins):
    # Out-degree: sum of outgoing connections for each cell
    outgoing_hub_strength[:, t] = np.nansum(temporal_net[:, :, t], axis=0)

# print(outgoing_hub_strength.shape)


# Find top hub cells (those with consistently high connectivity)
mean_hub_strength = np.nanmean(outgoing_hub_strength, axis=1)
top_hub_indices = np.argsort(mean_hub_strength)[-10:]  # Top 10 hub cells

# print(f"\nNetwork measures calculated:")
# print(f"  Global strength range: {np.nanmin(network_strength):.6f} to {np.nanmax(network_strength):.6f}")
# print(f"  Avg connection strength range: {np.nanmin(avg_connection_strength):.6f} to {np.nanmax(avg_connection_strength):.6f}")
# print(f"  Network density range: {np.min(network_density):.3f} to {np.max(network_density):.3f}")
# print(f"  Top outgoing hub cells (indices): {top_hub_indices}")
# print(f"  Top outgoing hub strengths: {mean_hub_strength[top_hub_indices]}")
# print(f"  Mean outgoing hub strength range: {np.nanmin(mean_hub_strength):.6f} to {np.nanmax(mean_hub_strength):.6f}")

# %%
## Interpolate network measures to match pupil sampling rate
from scipy import interpolate

# Interpolate network measures to pupil time points
network_strength_interp = interpolate.interp1d(network_time_array, network_strength, 
                                               kind='linear', bounds_error=False, fill_value=0)(timePupilEphys)

avg_connection_strength_interp = interpolate.interp1d(network_time_array, avg_connection_strength, 
                                                      kind='linear', bounds_error=False, fill_value=0)(timePupilEphys)

network_density_interp = interpolate.interp1d(network_time_array, network_density, 
                                              kind='linear', bounds_error=False, fill_value=0)(timePupilEphys)

# Interpolate hub cell activities to pupil time points
hub_activity_interp = np.zeros((len(top_hub_indices), len(timePupilEphys)))
for i, hub_idx in enumerate(top_hub_indices):
    hub_activity_interp[i, :] = interpolate.interp1d(network_time_array, outgoing_hub_strength[hub_idx, :], 
                                                      kind='linear', bounds_error=False, fill_value=0)(timePupilEphys)

# print(f"Network interpolation completed:")
# print(f"  Network measures interpolated to {len(timePupilEphys)} pupil time points")
# print(f"  Hub activities shape: {hub_activity_interp.shape}")
# print(f"  Interpolated network strength range: {np.nanmin(network_strength_interp):.6f} to {np.nanmax(network_strength_interp):.6f}")

# %%
## Update comprehensive aligned data structure with network data
neural_pupil_data['network'] = {
    # Raw network data
    'adjacency_matrix': temporal_net,  # (n_cells, n_cells, n_time_bins)
    'time_bins_raw': network_time_array,
    'bin_duration': network_bin_duration,
    
    # Network measures (interpolated to pupil sampling rate)
    'global_strength': network_strength_interp,
    'avg_connection_strength': avg_connection_strength_interp,
    'network_density': network_density_interp,
    
    # Hub cell data
    'hub_indices': top_hub_indices,
    'hub_activity': hub_activity_interp,  # (n_hubs, n_timepoints)
    'hub_mean_strength': mean_hub_strength[top_hub_indices],
    
    # Raw measures (not interpolated)
    'global_strength_raw': network_strength,
    'avg_connection_strength_raw': avg_connection_strength,
    'network_density_raw': network_density,
    'hub_strength_raw': outgoing_hub_strength
}

# %%
# Add cross-modal correlations
neural_pupil_data['correlations'] = {}

# Calculate correlations between pupil and network measures
from scipy.stats import pearsonr

# Remove NaN values for correlation calculation
valid_indices = ~(np.isnan(network_strength_interp) | np.isnan(pupilDiameter))
if np.sum(valid_indices) > 10:
    pupil_clean = pupilDiameter[valid_indices]
    network_strength_clean = network_strength_interp[valid_indices]
    
    corr_pupil_network, p_pupil_network = pearsonr(pupil_clean, network_strength_clean)
    neural_pupil_data['correlations']['pupil_network_strength'] = {
        'correlation': corr_pupil_network,
        'p_value': p_pupil_network
    }
    
    # Correlation between pupil and network density
    network_density_clean = network_density_interp[valid_indices]
    corr_pupil_density, p_pupil_density = pearsonr(pupil_clean, network_density_clean)
    neural_pupil_data['correlations']['pupil_network_density'] = {
        'correlation': corr_pupil_density,
        'p_value': p_pupil_density
    }
    
#     print(f"Cross-modal correlations calculated:")
#     print(f"  Pupil vs Network Strength: r={corr_pupil_network:.3f}, p={p_pupil_network:.3f}")
#     print(f"  Pupil vs Network Density: r={corr_pupil_density:.3f}, p={p_pupil_density:.3f}")

# print(f"\nComprehensive neural-pupil-network data structure updated:")
# print(f"  Pupil data points: {len(neural_pupil_data['pupil']['diameter'])}")
# print(f"  Spike data points: {len(neural_pupil_data['spikes']['time_bins'])}")
# print(f"  Network measures: {len(neural_pupil_data['network']['global_strength'])}")
# print(f"  All aligned to the same temporal resolution")

# %%
## Create comprehensive neural-pupil-network visualization
# plt.figure(figsize=(20, 16))

# # Extract data for plotting
# time_sec = neural_pupil_data['pupil']['time_sec']
# pupil_diameter = neural_pupil_data['pupil']['diameter'] 
# population_fr = neural_pupil_data['spikes']['population_firing_rate']
# network_strength = neural_pupil_data['network']['global_strength']
# network_density = neural_pupil_data['network']['network_density']
# hub_activity = neural_pupil_data['network']['hub_activity']
# event_times = neural_pupil_data['stimuli']['event_times']

# # Plot 1: Pupil diameter
# plt.subplot(6, 1, 1)
# plt.plot(time_sec, pupil_diameter, 'b-', alpha=0.8, linewidth=1, label='Pupil diameter')
# for i, stim_time in enumerate(event_times[:20]):
#     plt.axvline(stim_time, color='red', alpha=0.3, linestyle='--', linewidth=0.8)
# plt.ylabel('Pupil Diameter\n(pixels)')
# plt.title('Multi-Modal Neural Activity: Pupil, Spikes, and Granger Causal Networks', fontsize=16)
# plt.legend(loc='upper right')
# plt.grid(True, alpha=0.3)
# plt.xlim([time_sec[0], time_sec[-1]])

# # Plot 2: Population firing rate
# plt.subplot(6, 1, 2)
# plt.plot(time_sec, population_fr, 'g-', alpha=0.8, linewidth=1, label='Population firing rate')
# plt.fill_between(time_sec, population_fr, alpha=0.3, color='green')
# for i, stim_time in enumerate(event_times[:20]):
#     plt.axvline(stim_time, color='red', alpha=0.3, linestyle='--', linewidth=0.8)
# plt.ylabel('Population\nFiring Rate (Hz)')
# plt.legend(loc='upper right')
# plt.grid(True, alpha=0.3)
# plt.xlim([time_sec[0], time_sec[-1]])

# # Plot 3: Network global strength
# plt.subplot(6, 1, 3)
# plt.plot(time_sec, network_strength, 'purple', alpha=0.8, linewidth=1, label='Network global strength')
# plt.fill_between(time_sec, network_strength, alpha=0.3, color='purple')
# for i, stim_time in enumerate(event_times[:20]):
#     plt.axvline(stim_time, color='red', alpha=0.3, linestyle='--', linewidth=0.8)
# plt.ylabel('Network\nGlobal Strength')
# plt.legend(loc='upper right')
# plt.grid(True, alpha=0.3)
# plt.xlim([time_sec[0], time_sec[-1]])

# # Plot 4: Network density
# plt.subplot(6, 1, 4)
# plt.plot(time_sec, network_density, 'orange', alpha=0.8, linewidth=1, label='Network density')
# plt.fill_between(time_sec, network_density, alpha=0.3, color='orange')
# for i, stim_time in enumerate(event_times[:20]):
#     plt.axvline(stim_time, color='red', alpha=0.3, linestyle='--', linewidth=0.8)
# plt.ylabel('Network\nDensity')
# plt.legend(loc='upper right')
# plt.grid(True, alpha=0.3)
# plt.xlim([time_sec[0], time_sec[-1]])

# # Plot 5: Hub cell activities (top 5)
# plt.subplot(6, 1, 5)
# colors = plt.cm.Set3(np.linspace(0, 1, 5))
# hub_indices = neural_pupil_data['network']['hub_indices']
# for i in range(min(5, len(hub_indices))):
#     plt.plot(time_sec, hub_activity[i, :], color=colors[i], alpha=0.8, linewidth=1, 
#              label=f'Hub cell {hub_indices[i]}')
# for i, stim_time in enumerate(event_times[:20]):
#     plt.axvline(stim_time, color='red', alpha=0.3, linestyle='--', linewidth=0.8)
# plt.ylabel('Hub Cell\nActivity')
# plt.legend(loc='upper right', ncol=5, fontsize=8)
# plt.grid(True, alpha=0.3)
# plt.xlim([time_sec[0], time_sec[-1]])

# # Plot 6: Cross-correlations
# plt.subplot(6, 1, 6)
# # Calculate sliding window correlations between all measures
# window_size = min(5000, len(time_sec) // 20)  # 5% of data or 5000 points
# stride = window_size // 4

# correlations_pupil_network = []
# correlations_spikes_network = []
# correlations_pupil_spikes = []
# correlation_times = []

# for i in range(0, len(time_sec) - window_size, stride):
#     end_idx = i + window_size
#     if end_idx <= len(time_sec):
#         # Get data windows
#         pupil_window = pupil_diameter[i:end_idx]
#         spikes_window = population_fr[i:end_idx]
#         network_window = network_strength[i:end_idx]
        
#         # Remove NaN values
#         valid_idx = ~(np.isnan(pupil_window) | np.isnan(spikes_window) | np.isnan(network_window))
        
#         if np.sum(valid_idx) > 10:  # Need enough valid points
#             pupil_clean = pupil_window[valid_idx]
#             spikes_clean = spikes_window[valid_idx]
#             network_clean = network_window[valid_idx]
            
#             # Calculate correlations
#             corr_pn, _ = pearsonr(pupil_clean, network_clean)
#             corr_sn, _ = pearsonr(spikes_clean, network_clean)
#             corr_ps, _ = pearsonr(pupil_clean, spikes_clean)
            
#             correlations_pupil_network.append(corr_pn)
#             correlations_spikes_network.append(corr_sn)
#             correlations_pupil_spikes.append(corr_ps)
#             correlation_times.append(time_sec[i + window_size//2])

# plt.plot(correlation_times, correlations_pupil_network, 'blue', linewidth=2, alpha=0.8, 
#          label='Pupil ↔ Network')
# plt.plot(correlation_times, correlations_spikes_network, 'green', linewidth=2, alpha=0.8, 
#          label='Spikes ↔ Network')
# plt.plot(correlation_times, correlations_pupil_spikes, 'red', linewidth=2, alpha=0.8, 
#          label='Pupil ↔ Spikes')

# plt.axhline(y=0, color='black', linestyle='-', alpha=0.3)
# for i, stim_time in enumerate(event_times[:20]):
#     plt.axvline(stim_time, color='red', alpha=0.3, linestyle='--', linewidth=0.8)

# plt.xlabel('Time (seconds)')
# plt.ylabel('Correlation\nCoefficient')
# plt.legend(loc='upper right')
# plt.grid(True, alpha=0.3)
# plt.xlim([time_sec[0], time_sec[-1]])
# plt.ylim([-1, 1])

# plt.tight_layout()
# plt.show()

# # Print comprehensive summary
# print(f"\n=== COMPREHENSIVE MULTI-MODAL ALIGNMENT SUMMARY ===")
# print(f"Data duration: {time_sec[-1] - time_sec[0]:.1f} seconds")
# print(f"Sampling rate: {1/pupil_sampling_interval:.1f} Hz")
# print(f"Total data points: {len(time_sec)}")
# print(f"Number of stimuli: {len(event_times)}")
# print(f"Number of cells: {len(cells)}")
# print(f"Hub cells identified: {len(hub_indices)}")

# print(f"\nCross-modal correlations:")
# for key, corr_data in neural_pupil_data['correlations'].items():
#     print(f"  {key}: r={corr_data['correlation']:.3f}, p={corr_data['p_value']:.3f}")

# print(f"\nNetwork measures:")
# print(f"  Global strength: {np.nanmin(network_strength):.3f} to {np.nanmax(network_strength):.3f}")
# print(f"  Network density: {np.min(network_density):.3f} to {np.max(network_density):.3f}")
# print(f"  Top hub cell activities: {np.min(hub_activity):.3f} to {np.max(hub_activity):.3f}")

# # %%
# print(neural_pupil_data.keys())

# %%
import pickle

output_path = Path(DATA_PATH_AROUSAL, "neural_pupil_network_aligned_temporally.pkl")
# Save the dictionary
with open(output_path, 'wb') as file:
    pickle.dump(neural_pupil_data, file)


