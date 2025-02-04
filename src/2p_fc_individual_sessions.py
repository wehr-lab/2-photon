## import libraries 
import os
import numpy as np
import matplotlib.pyplot as plt
import ephys_analysis 
import behavior_analysis
import seaborn as sns
import igraph as ig

sessionPaths = ["/Volumes/projects/2P5XFAD/JarascopeData/wehr3204/01-17-25-002/suite2p/plane0", "/Volumes/projects/2P5XFAD/JarascopeData/wehr3204/01-22-25-002/suite2p/plane0", "/Volumes/projects/2P5XFAD/JarascopeData/wehr3204/01-30-25-002/suite2p/plane0"]

dFFs = {} ## this will store the dFF values for each session

for session in sessionPaths:
    ## load data 
    F, Fneu, spks, stat, ops, iscell = ephys_analysis.load2P(session) 

    print("Data loaded for session: ", session)

    ## NOTE: data loading is taking a long time. This is because the data is being loaded from a network drive. This will be the same when we work on talapas. Maybe there is a better way? 

    indicesCells = ephys_analysis.getIndices(iscell, prob=0.5)
    fCells = ephys_analysis.extractCells(F, indicesCells)
    fneuCells = ephys_analysis.extractCells(Fneu, indicesCells)

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
    dFF = ephys_analysis.computeDFFSlide(fCorrectedShifted, window=200, percentile=20) ## note that the window size and the percentile are parameters are totally arbitrary for now 
    print("Sliding window dF/F computed for session: ", session)

    # ## save data for further processing 
    # np.save(os.path.join(DATA_PATH_SAVE, "dF_F_100824_sliding_f_nought.npy"), dFF)
    dFFConvolve = ephys_analysis.convolve_spikes(dFF, sigma=10)
    dFFs[session] = dFFConvolve

fcMatrix = {}
fcMatrixBin = {}
for session, dFF in dFFs.items():
    ## compute the correlation matrix
    fcMatrix[session] = np.corrcoef(dFF)
    print("Functional connectivity matrix computed for session: ", session)

    ## binarize the correlation matrix
    fcMatrixBin[session] = ephys_analysis.binarize_matrix(fcMatrix[session], threshold=0.1)
    ## remove autocorrelations
    np.fill_diagonal(fcMatrixBin[session], 0)

## plot the functional connectivity matrices
# fig, axs = plt.subplots(1, 3, figsize=(15, 5))
for i, (session, fc) in enumerate(fcMatrixBin.items()):
    sns.clustermap(fc, method="ward", cmap="coolwarm", square=True,  xticklabels=False, yticklabels=False)
    plt.title(session.split("/")[-3])



## now stuff about the graph theory analysis -> I am using the igraph library for this


