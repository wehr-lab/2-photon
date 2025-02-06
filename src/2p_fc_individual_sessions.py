## import libraries 
import os
import numpy as np
import matplotlib.pyplot as plt
import ephys_analysis 
import behavior_analysis
import seaborn as sns
import igraph as ig

# sessionPaths = ["/Volumes/projects/2P5XFAD/JarascopeData/wehr3204/01-17-25-002/suite2p/plane0", "/Volumes/projects/2P5XFAD/JarascopeData/wehr3204/01-22-25-002/suite2p/plane0", "/Volumes/projects/2P5XFAD/JarascopeData/wehr3204/01-30-25-002/suite2p/plane0"]

sessionPaths = ["/Volumes/projects/2P5XFAD/JarascopeData/wehr3204/01-17-25-002/suite2p/plane0"]

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

    ## compute the 20th percentile of the correlation matrix
    matToProcess = fcMatrix[session]
    upperTriangular = matToProcess[np.triu_indices_from(matToProcess, k=1)] 
    thresholdPercentile = np.percentile(upperTriangular.flatten(), 60) 
    ## binarize the correlation matrix at the threshold percentile 
    fcMatrixBin[session] = ephys_analysis.binarize_matrix(matToProcess, threshold=thresholdPercentile)
    ## remove autocorrelations
    np.fill_diagonal(fcMatrixBin[session], 0)

## plot the functional connectivity matrices
# fig, axs = plt.subplots(1, 3, figsize=(15, 5))
for i, (session, fc) in enumerate(fcMatrixBin.items()):
    sns.clustermap(fc, method="ward", cmap="coolwarm", square=True,  xticklabels=False, yticklabels=False)
    plt.title(session.split("/")[-3])



## now stuff about the graph theory analysis -> I am using the igraph library for this

## create a graph object from the functional connectivity matrix
graphs = [ig.Graph.Adjacency(fcMatrixBin[session].tolist()) for session in sessionPaths]

graphMatrics = {}

for i, graph in enumerate(graphs):
    ## compute the degree of each node
    storeDict = {}
    storeDict["degree"] = graph.degree()
    ## compute the clustering coefficient of each node
    storeDict["clustering"] = graph.transitivity_local_undirected()
    ## compute the betweenness centrality of each node
    storeDict["betweenness"] = graph.betweenness()
    ## compute the eigenvector centrality of each node
    storeDict["eigenvector"] = graph.eigenvector_centrality()
    storeDict["path_length"] = graph.average_path_length()
    graphMatrics[i] = storeDict


## plot all the graph metrics
degrees = [graphMatrics[i]["degree"] for i in range(len(graphMatrics))]
clusterings = [graphMatrics[i]["clustering"] for i in range(len(graphMatrics))]
pathLengths = [graphMatrics[i]["path_length"] for i in range(len(graphMatrics))]

fig, axs = plt.subplots(1, 3, figsize=(15, 5))
axs[0].boxplot(degrees)
axs[0].set_title("Degree")
axs[1].boxplot(clusterings)
axs[1].set_title("Clustering")
axs[2].boxplot(pathLengths)
axs[2].set_title("Path Length")
plt.show()


