## import standard libraries 
import os 
import sys 
import glob 
import numpy as np 
import scipy.io as sio
import matplotlib.pyplot as plt

## add path to custom libraries
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

## import the custom libraries 
from utils.extract_vals import extract_value 
from utils.strcmp import strcmp
from behavior import Behavior
from events import Events
from session import Session
from cell import Cell 
from cell_ensemble import CellEnsemble
from events_process import EventsProcess
import ephys_analysis


## define variables 
EPHYS_DIR_PATH = "/Volumes/projects/5XFAD/Rig3Phys/2025-02-03_8-37-45_mouse-3359/2025-02-03_08-37-56mouse-3359/"

## load the sorted units file
sortedUnitsFilePath = glob.glob(os.path.join(EPHYS_DIR_PATH, "Sort*.mat"))[0]
sortedUnitsFile = sio.loadmat(sortedUnitsFilePath, struct_as_record=True) 
sortedUnitsFile = sortedUnitsFile['SortedUnits']
sortedUnitsFile = sortedUnitsFile.reshape(-1) 

## create a list of Cell objects to pass to CellEnsemble class 
cells = []
for cell in sortedUnitsFile:
    cells.append(Cell(cell))

## create a CellEnsemble object
cellEnsemble = CellEnsemble(cells)

cellEnsembleDF = cellEnsemble.toDataFrame()

## for sanity check let's plot the spike raster for all the cells 
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
for cellnum, cellSpike in zip(cellEnsembleDF["cellnum"], cellEnsembleDF["spiketimes"]):
    y = np.ones_like(cellSpike) * cellnum
    ax.plot(cellSpike, y, 'k.', markersize=1)

ax.set_xlabel('Time (ms)')
ax.set_ylabel('Cell Number')
plt.show()
plt.close()

## create a matrix of the spike times for each cell
ms = 1e3 
spikeTimes = cellEnsembleDF["spiketimes"].values
spikeTimes = spikeTimes * ms ## convert to ms

minTime = 0
maxtime = np.max([np.max(cellSpike) for cellSpike in spikeTimes]) + 1000 ## added 1000 ms so that 1 sec padding added 
numCells = len(spikeTimes)
spikeMatrix = np.zeros((numCells, int(maxtime - minTime)))

for cellnum, cellSpike in enumerate(spikeTimes):
    spikeMatrix[cellnum, cellSpike.astype(int)] = 1

## run convolution of the entire data 
spikeMatrixConv = ephys_analysis.convolve_spikes(spikeMatrix, sigma=1000) ## why sigma = 100?

## plot the convolved spike matrix
plt.plot(spikeMatrixConv[0, :])

## subset the data 




