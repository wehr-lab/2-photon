## import libraries 
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import ephys_analysis 
import behavior_analysis

## define path variables -> note to make these paths genaralizable by doing sys.argv[]
DATA_PATH = "/Volumes/Projects/2P5XFAD/JarascopeData/wehr2867/10-08-24-001/suite2p/plane0/" ## path to where the data is mounted; can also sys.argv[1]
# DATA_PATH_SAVE = "/Users/praveslamichhane/Desktop/alzheimers/data/wehr2867"


sessions = [] ## list of paths to the data for different sessions
dFFs = {}
for session in sessions:
    ## load data 
    F, Fneu, spks, stat, ops, iscell = ephys_analysis.load2P(session) 

    fCells = ephys_analysis.extractCells(F, iscell, prob=0.8)
    fneuCells = ephys_analysis.extractCells(Fneu, iscell, prob=0.8)

    ## neuropil correction
    fCorrected = ephys_analysis.correct_neuropil(fCells, fneuCells, npil_coeff=0.7)

    ## remove negative values from the corrected fluorescence values
    fCorrectedShifted = ephys_analysis.absFluoresce(fCorrected)

    ## another plot to visualize the time series of the corrected fluorescence values and to select a threshold
    # plt.plot(fCorrectedShifted[0, :])

    ## draw a horizontal line at the 10th percentile
    # plt.axhline(np.percentile(fCorrectedShifted[0, :], 20), color="red")
    # plt.show() 

    ## compute dF/F using a sliding window
    dFF = ephys_analysis.computeDFFSlide(fCorrectedShifted, window=100, percentile=20)


    # ## save data for further processing 
    # np.save(os.path.join(DATA_PATH_SAVE, "dF_F_100824_sliding_f_nought.npy"), dFF)

    dFFs[session] = dFF

fcMatrix = {}
for session, dFF in dFFs.items():
    ## compute the correlation matrix
    fcMatrix[session] = np.corrcoef(dFF)

