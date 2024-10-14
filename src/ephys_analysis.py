## import libraries 
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def roi_indices(iscell, prob=None):
    """Return the indices of the ROIs that are cells."""
    if prob is not None:
        return np.where(iscell[:, 1] >= prob)[0]
    return np.where(iscell[:, 0] == 1)[0]

def extract_cells(data, indices):
    """Extract the data for the specified cell indices."""
    return data[indices]

def correct_neuropil(cell_data, npil_data, npil_coeff=0.7):
    """Correct the neuropil contamination in the cell data using the neuropil data."""
    return cell_data - npil_coeff * npil_data

def compute_dff(data, f0=None):
    """Compute the dF/F for the data using the specified baseline."""
    if f0 is None:
        f0 = np.nanmean(data, axis=1)
    return (data - f0[:, np.newaxis]) / f0[:, np.newaxis]

## define function to plot spikes
def plot_spikes_subset(data, time_points=50, cell_indices=None, y_lim=None):
    """Plot subset of spikes for a subset of cells. This is to see how the spikes look after dF/F correction."""
    if cell_indices is not None:
        cell_indices = np.array(list(cell_indices))
        for i in cell_indices:
            if i >= data.shape[0]:
                raise ValueError("Cell index out of range")
        someCells = np.array(cell_indices)
    else:
        someCells = np.array(list(np.arange(10)))
    
    num_rows = int(len(someCells) // np.sqrt(len(someCells)))
    num_cols = int(np.ceil(len(someCells) / num_rows))
    # print(num_rows, num_cols)

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(20, 20))  
    axes = axes.flatten()

    for i, cell_idx in enumerate(someCells):
        axes[i].plot(data[int(cell_idx), :time_points])
        axes[i].set_title("Cell: " + str(someCells[i] + 1))
        axes[i].set_xlabel("Time (s)")
        axes[i].set_ylabel("dF/F")
        if y_lim is not None:
            axes[i].set_ylim(y_lim)

    plt.tight_layout()
    plt.show()

def plot_trace(data, cell_indices=None, y_lim=None):
    """Plot the traces of the specified cells. These are the entire dF/F traces."""
    if cell_indices is not None:
        cell_indices = np.array(list(cell_indices))
        for i in cell_indices:
            if i >= data.shape[0]:
                raise ValueError("Cell index out of range")
        someCells = np.array(cell_indices)
    else:
        someCells = np.array(list(np.arange(10)))
    
    num_rows = int(len(someCells) // np.sqrt(len(someCells)))
    num_cols = int(np.ceil(len(someCells) / num_rows))
    # print(num_rows, num_cols)

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(20, 20))  
    axes = axes.flatten()

    for i, cell_idx in enumerate(someCells):
        axes[i].plot(data[int(cell_idx)])
        axes[i].set_title("Cell: " + str(someCells[i] + 1))
        axes[i].set_xlabel("Time (s)")
        axes[i].set_ylabel("dF/F")
        if y_lim is not None:
            axes[i].set_ylim(y_lim)

    plt.tight_layout()
    plt.show()

def subset_data(data, num_points):
    """Return a subset of the data with the specified number of time points (frames)."""
    num_cells = data.shape[0]
    data_subset = np.zeros((num_cells, num_points))

    for i in range(num_cells):
        data_subset[i] = data[i][:num_points]

    return data_subset