## import standard libraries 
import os 
import sys 
import glob 
import numpy as np 
import scipy.io as sio

## add path to custom libraries -> can directly import the libraries, no need to add path 
from utils.extract_vals import extract_value 
from utils.strcmp import strcmp
from behavior import Behavior
from events import Events
from session import Session
from cell import Cell 
from cell_ensemble import CellEnsemble
from events_process import EventsProcess

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

## load the notebook file 
notebookFilePath = glob.glob(os.path.join(ephysDirPath, "note*.mat"))[0]
notebookFile = sio.loadmat(notebookFilePath, struct_as_record=True)
stimLog = notebookFile["stimlog"].reshape(-1)

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

## for sanity check let's plot the spike raster for all the cells 
# import matplotlib.pyplot as plt

# fig, ax = plt.subplots(1, 1, figsize=(10, 10))
# for cellnum, cellSpike in zip(cellEnsembleDF["cellnum"], cellEnsembleDF["spiketimes"]):
#     y = np.ones_like(cellSpike) * cellnum
#     ax.plot(cellSpike, y, 'k.', markersize=1)

# ax.set_xlabel('Time (s)')
# ax.set_ylabel('Cell Number')
# plt.show()

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

## Store event specific information in a EventsProcess object. This is different from the Events object as it will store the frame numbers, spike times, pupil diameter, etc. for each event. 

eventsProcess = []
for trial in range(len(session)):
    if strcmp(session.get_event(trial).type, ('tone', 'whitenoise', 'silentsound')):
        if hasattr(session.get_event(trial), 'soundcard_trigger_timestamp_sec'):
            eventStartTime = session.get_event(trial).soundcard_trigger_timestamp_sec * ms ## convert time in seconds to milliseconds
        else:
            raise ValueError('Event does not have soundcard_trigger_timestamp_sec attribute.')

        ## now store all the event specific information in the EventsProcess object
        currEvent = EventsProcess(session.get_event(trial)) ## initialize with the Events object
        if not strcmp(session.get_event(trial).type, 'tone'):
            currEvent.event.frequency = np.nan 
        
        ## event (trial) specific start and stop times 
        currEvent.OEStartTime = np.ceil(eventStartTime + timeWindow[0]) ## start time in milliseconds
        currEvent.OEEndTime = np.ceil(eventStartTime + timeWindow[1])
        currEvent.freqStart = np.ceil(eventStartTime) ## start time in of stimulus in milliseconds
        currEvent.freqEnd = np.ceil(currEvent.freqStart + currEvent.event.duration) ## end time of stimulus in milliseconds; note that the duration is already in milliseconds

        ## find the frame number for the start and end time of the stimulus
        currEvent.frameRate = extract_value(reyeFile['vid'][0][0]['framerate']) ## lots of nesting, but checked and its correct
        frameRateInMs = (1/ms) * currEvent.frameRate ## frame rate in milliseconds
        currEvent.firstFrame = np.int32(np.ceil(sessionFrameNums[trial] + timeWindow[0] * frameRateInMs)) ## frame number for the start time of the stimulus
        currEvent.lastFrame = np.int32(np.ceil(sessionFrameNums[trial] + timeWindow[1] * frameRateInMs))

        ## extracting the spiketimes now 
        for cellNum in range(len(cellEnsemble)):
            currCell = cellEnsemble.get_cell(cellNum)
            currEvent.spikeTimes.append(currCell.spiketimes[(currCell.spiketimes >= currEvent.OEStartTime) & (currCell.spiketimes <= currEvent.OEEndTime)]) ## extract the spike times within the time window for each cell and append

        ## extract the pupil diameter for the event
        currEvent.pupilDiameter = [pupilDiameter[currEvent.firstFrame:currEvent.lastFrame]] ## extract the pupil diameter for the event

        ## time of frames 
        currEvent.frameNums = np.arange(currEvent.firstFrame, currEvent.lastFrame + 1) ## frame numbers for the event
        currEvent.frameTimes = np.int32(np.ceil((currEvent.frameNums - intercept) / slope)) ## convert frame numbers to OE time in milliseconds 

        ## store timeWindow and frameWindow
        currEvent.timeWindow = timeWindow
        currEvent.frameWindow = [currEvent.firstFrame, currEvent.lastFrame]

        ## append the currEvent to the eventsProcess list
        eventsProcess.append(currEvent)





