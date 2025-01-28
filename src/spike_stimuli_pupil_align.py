## import standard libraries 
import os 
import sys 
import glob 
import numpy as np 
import scipy.io as sio

## add path to custom libraries -> can directly import the libraries, no need to add path 
from utils.extract_vals import extract_value 
from behavior import Behavior
from events import Events
from session import Session
from cell import Cell 
from cell_ensemble import CellEnsemble

## define variables 
BONSAI_DIR_PATH = "/Volumes/projects/ToneDecodingArousal/2022-05-09_13-52-25_mouse-0951"

## load the .mat files 
behaviorFilePath = glob.glob(os.path.join(BONSAI_DIR_PATH, "Beh*.mat"))[0]
behaviorFile = Behavior(behaviorFilePath)

headFile = behaviorFile.get_head_data() 
reyeFile = behaviorFile.get_reye_data()

## define the ephys dir 
sessionName = reyeFile["dirName"]
sessionName = extract_value(sessionName, 'dirName')
ephysDirPath = os.path.join(BONSAI_DIR_PATH, sessionName)

## load the sortedUnits file
sortedUnitsFilePath = glob.glob(os.path.join(ephysDirPath, "Sort*.mat"))[0]
sortedUnitsFile = sio.loadmat(sortedUnitsFilePath, struct_as_record=True)
sortedUnitsFile = sortedUnitsFile['SortedUnits']
sortedUnitsFile = sortedUnitsFile.reshape(-1)

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

pupilDiameter = behaviorFile.get_pupilDiameter_data()
cellEnsembleDF = cellEnsemble.toDataFrame()

## for sanity check let's plot the spike raster for all the cells 
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 1, figsize=(10, 10))
for cellnum, cellSpike in zip(cellEnsembleDF["cellnum"], cellEnsembleDF["spiketimes"]):
    y = np.ones_like(cellSpike) * cellnum
    ax.plot(cellSpike, y, 'k.', markersize=1)

ax.set_xlabel('Time (s)')
ax.set_ylabel('Cell Number')
plt.show()




