import numpy as np 
import h5py

def extract_value(arr, fieldname=None):
    '''Only run this function after importing numpy as np. 
    
    input: arr, a numpy array or list
    output: the very first value from the array or list i.e. arr[0][0]...
    '''
    if fieldname == 'spiketimes':
        if isinstance(arr, np.ndarray):
            return np.array(arr[0])
        elif isinstance(arr, list):
            return arr  # Return the list as is
        else:
            return None  # Handle unsupported type

    if fieldname == 'PupilDiameter':
        return arr[0][0][:, 0] ## I was lazy to write a proper function. This is a temporary solution 

    while True:
        ## check if the current item is a NumPy array or a list
        if isinstance(arr, (np.ndarray, list)):
            if len(arr) > 0:  # Ensure the array or list is not empty
                arr = arr[0]
            else:
                return None  # Return None for empty arrays or lists
        # Check if the value is a NumPy integer or float
        elif isinstance(arr, (np.integer, np.floating, np.str_)):
            return arr.item()  # Convert to plain Python type
        # Check if the value is a plain Python integer or float
        elif isinstance(arr, (int, float)):
            return arr
        else:
            return None  # Return None for unsupported types


def strcmp(str1, str2):
    if str2 is not str or type(str2) in [list, tuple]:
        return str1 in str2 
    else:
        return str1 == str2

def hdf5ToDict(hdf5FilePath):
    '''Convert a HDF5 file (including groups and datasets) into a nested dictionary.

    Args:
        hdf5FilePath: full path to the HDF5 file

    Returns:
        A nested dictionary mirroring the HDF5 hierarchy, with datasets loaded into arrays.
    '''
    def _recurse(h5obj):
        data = {}
        for key, item in h5obj.items():
            if isinstance(item, h5py.Dataset):
                data[key] = item[()]
            elif isinstance(item, h5py.Group):
                data[key] = _recurse(item)
        return data

    with h5py.File(hdf5FilePath, 'r') as f:
        return _recurse(f)


def get_unique_dtypes(arr):
    '''Get unique data types from a numpy array or list.

    Args:
        arr: a numpy array or list

    Returns:
        A set of unique data types found in the array.
    '''
    if isinstance(arr, np.ndarray):
        return set(arr.dtype.type for x in np.ravel(arr))
    elif isinstance(arr, list):
        return set(type(x) for x in np.ravel(arr))
    else:
        raise TypeError("Input must be a numpy array or a list.")