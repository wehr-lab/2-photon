import numpy as np 
import networkx as nx

def makeSymmetric(G):
    """
    Make a square matrix symmetric by averaging it with its transpose.
    
    Parameters:
    G (numpy.ndarray): The input square matrix.
    
    Returns:
    numpy.ndarray: The symmetric version of the input matrix.
    """
    return (G + G.T) / 2

def makeBinary(G, threshold=0.5):
    """
    Convert a matrix to binary based on a threshold.
    
    Parameters:
    G (numpy.ndarray): The input matrix.
    threshold (float): The threshold value for binarization.
    
    Returns:
    numpy.ndarray: The binary version of the input matrix.
    """
    return (G >= threshold).astype(int)


def safeRichClub(G, k=None):
    """
    Wrapper for the rich-club coefficient calculation in NetworkX to return only the RC value for a specific k.
    """ 
    

def richClubDegreeN(G, k):
    pass 

def matDensity(mat):
    """Calculate the connection density of a 3D connectivity matrix."""
    total_connections = mat.size
    non_zeros = np.count_nonzero(mat)
    density = non_zeros / total_connections
    return density