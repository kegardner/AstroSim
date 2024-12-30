


from collections.abc import Callable
import csv
import os
import sys
from typing import Any

import numpy as np
from result import Ok, Err, Result


def memoize(closure: Callable[[], list[list[Any]]], filenames : list[str] | str, force_recalc_cla = None) -> Result[np.ndarray, str]:
    """
    
    A function for memoizing 2d collections into csv files.
    Force a recalc/rewrite of all uses of memoize within a program using CLA "--no-cache"

    Parameters
    ----------
    closure : func -> list of list of Any
            a lambda like "lambda someFunc(data)" where "data" is being presented from the enclosing context
            output of this function _must_ be directly compatible with CSVs.
    filenames : list of str or str
            the filename(s) of the output data. When data is a collection of 2d arrays/lists/etc, this MUST be a collection of strings. If the closure only returns one 2d array, it may be a string, and the closure may return the data as the raw 2d array (as opposed to wrapping it in a list)
    force_recalc_cla : str | None, default = None
            an optional parameter that, when populated, directs the memoizer to check the command line arguments for the presence of that string. If present, the closure is called and the results are written regardless of whether the corresponding files are present or not.
    Examples:
        >>> def someFunc(data : Any) -> list[list[int]]:
                return [[1,2,3],[4,5,6]]

        >>> memoize(lambda : closure("hello"), "output_hello")
            [[1,2,3],[4,5,6]]

        memoize will look to see if "output_hello.csv" exists, running the closure if it does not (writing the results of that call to "output_hello.csv"), and reading the results in if it does. In both cases, it will return the same results.
    """
    

    print(f"Checking Cache State: {filenames}")
    if type(filenames) is not list:
        filenames = [filenames]
    if ((not all(map(os.path.isfile, filenames))) 
        or sys.argv.count("--no-cache") > 0 
        or (force_recalc_cla is not None and sys.argv.count(force_recalc_cla) > 0)):
        print("Missing data, running closure...")
        #convert to np array here, since thats what we always need anyway
        try:
            output_arrays = np.array(closure())
            result = _writeOutput(output_arrays, filenames)
        except Exception as e:
            result = Err(f"Output of closure was not compatible with numpy array conversion: {e}")

    else:
        result = _readCachedResult(filenames)
    return result

def _writeOutput(output, filenames) -> Result[np.ndarray, str]:
    """
    
    """
    
    def _write(data, filename):
        print(f"writing {filename}")
        with open(filename, "w+") as f :
            w = csv.writer(f, delimiter="|",lineterminator="\n")
            for row in data:
                w.writerow(row)
    try:
        print(output.shape)
        match len(output.shape):
            case 2: 
                _write(output, filenames)
            case 3: 
                list(map(lambda data, name : _write(data, name), output, filenames)) 
            case a:
                return Err(f"Output must be of dimension 2 or 3, but is of dimension {a}")
    except Exception as e:
        return Err(f"Writing of output(s) failed with: {e}")
    
    return Ok(output)

def _readCachedResult(filenames) -> Result[np.ndarray, str]:
    arrays = []
    try:
        for file in filenames:
            with open(file, "r") as f :
                r = csv.reader(f, delimiter="|",lineterminator="\n")
                temp = []
                for row in r:
                    temp.append(row)
                arrays.append(temp)
    except Exception as e:
        return Err(f"Reading of cached results failed with: {e}")
    
    return Ok(np.array(arrays))