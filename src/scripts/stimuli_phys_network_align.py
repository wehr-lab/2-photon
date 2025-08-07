import os 
import sys
from pathlib import Path
import scipy.io as sio
import numpy as np
import networkx as nx
import h5py 
import glob 
import matplotlib.pyplot as plt

## import custom modules 
repoDir = Path(__file__).resolve().parent.parent.parent ## dir of the repo root
# print("repoDir:", repoDir)
configDir = repoDir / "config"
sys.path.append(str(configDir))

## define path variables 
from settings import DATA_PATH, CODE_PATH
sys.path.append(str(CODE_PATH))

from utils.util_funcs import extract_value 
from ephys.behavior import Behavior
from ephys.cell import Cell
from ephys.cell_ensemble import CellEnsemble
from ephys.events import Events
from ephys.events_process import EventsProcess
from ephys.session import Session

## load .mat files and typecast to python classes 
DATA_PATH = Path(DATA_PATH["toneDecode"], "2022-05-10_14-02-30_mouse-0956")
# print("DATA_PATH:", DATA_PATH)

## load the .mat files 
behaviorFilePath = glob.glob(os.path.join(DATA_PATH, "Beh*.mat"))[0]
behaviorFile = Behavior(behaviorFilePath)

headFile = behaviorFile.get_head_data() 
reyeFile = behaviorFile.get_reye_data()

## define the ephys dir 
sessionName = reyeFile["dirName"]
sessionName = extract_value(sessionName, 'dirName')
ephysDirPath = os.path.join(DATA_PATH, sessionName)

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
sessionDF = session.toDataFrame()

pupilDiameter = behaviorFile.get_pupilDiameter_data()
cellEnsembleDF = cellEnsemble.toDataFrame()

## Extract stimulus onset and offset times 
ms = 1e3 ## convert time in seconds to milliseconds
timeWindow = [-500, 500] ## time window in milliseconds to extract spikes around the stimulus onset

## define the frame numbers when the first and last stimuli were presented 
sessionFrameNums = reyeFile["TTs"][0][0][0] ## frame numbers when the stimuli were presented
firstEventFrame = sessionFrameNums[0] ## first frame number when the first stimulus was presented
lastEventFrame = sessionFrameNums[-1] ## last frame number when the last stimulus was presented

## extract the time for OpenEphys for the first and last stimulus presentation
OEFirstEventTime = session.get_event(0).soundcard_trigger_timestamp_sec ## OE time in seconds when the first stimulus was presented
OELastEventTime = session.get_event(-1).soundcard_trigger_timestamp_sec ## OE time in seconds when the last stimulus was presented 

## to track the frame numbers with OE time, find the slope and intercept of FrameNumber vs OE time
slope = (lastEventFrame - firstEventFrame) / (OELastEventTime - OEFirstEventTime)
intercept = firstEventFrame - (slope * OEFirstEventTime)








print("run complete!")