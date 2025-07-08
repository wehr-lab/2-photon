import os 
import sys
from pathlib import Path
import glob
import re 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

## define the path to the data
parentPath = Path("/Volumes/Projects/2P5XFAD/JarascopeData/wehr3378")
suite2pPathStr = "suite2p/plane0"
trackPathStr = "track2p/track2p-03-27-002-004--Loc14"
trackPath = os.path.join(parentPath, trackPathStr) 

matchMatPath = glob.glob(os.path.join((trackPath), "*match*.npy"))
matchMat = np.load(matchMatPath[0], allow_pickle=True) ## this is the shape of the match matrix from the track2p algo 

nCellsTracked, nSessionsTracked = matchMat.shape

sessionsTracked = {}
for num, session in enumerate(sorted(os.listdir(os.path.join(trackPath, "matched_suite2p")))):
    if not re.match(r"^\d", session):
        continue
    sessionsTracked[num] = session


## load iscell from suite2p output 
iscellPaths = {}
for i, p in sessionsTracked.items():
    iscellPaths[i] = os.path.join(parentPath, p, suite2pPathStr, "iscell.npy")

## count of iscells for each sessions
iscellCounts = {}
for i, p in sessionsTracked.items():
    iscell = np.load(iscellPaths[i], allow_pickle=True)
    # iscellCounts[i] = np.nonzero(iscell[:, 0])[0]
    iscellCounts[p] = np.sum(iscell[:, 0])
    print(f"Session {p}: {iscellCounts[p]} cells")


## making sure the day 1 recording is the one referenced in the match matrix
firstVec = (matchMat[:, 0] == None) ## where are the none values i.e. all the none cells 
secVec = (matchMat[:, 1:] != None) ## are any cells not none here? 
weirdCells = np.where(firstVec & np.any(secVec == True, axis=1)) ## this is the weird cells
print(f"Weird cells: {weirdCells}") ## if empty then good 

## cerate a dataframe for cell maps 
cellMap = pd.DataFrame(matchMat, columns=sessionsTracked.values(), index=range(nCellsTracked))

print(cellMap)

## track the cells lost across sessions 
cellsLost = {}

for i, cell in enumerate(cellMap[cellMap.columns[0]].values): ## cause this is our master cell map
    if cell != None:
        rowVals = cellMap.iloc[i].values
        try:
            firstNone = np.where(rowVals == None)[0][0]
            if (cellMap.columns[firstNone - 1], cellMap.columns[firstNone]) not in cellsLost.keys():
                cellsLost[(cellMap.columns[firstNone - 1], cellMap.columns[firstNone])] = []
                cellsLost[(cellMap.columns[firstNone - 1], cellMap.columns[firstNone])].append(cell)
            else:
                cellsLost[(cellMap.columns[firstNone - 1], cellMap.columns[firstNone])].append(cell)
        except:
            continue 
print(cellsLost.keys()) 




