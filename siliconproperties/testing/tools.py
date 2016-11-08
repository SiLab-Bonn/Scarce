import inspect
import os
import tables as tb
from collections import OrderedDict
import numpy as np
import itertools

FIXTURE_FOLDER = 'fixtures'

def _call_function_with_args(function, **kwargs):
    ''' Calls the function with the given kwargs
    and returns the result in a numpy array  '''

    # Create all combinations of arguments from list parameters
    # This is ugly but avoids recursion
    call_values = []
    fixed_arguments = []
    fixed_arguments_pos = []
    for index, values in enumerate(kwargs.values()):
        if isinstance(values, list):
            call_values.extend([values])
        else:
            fixed_arguments.append(values)
            fixed_arguments_pos.append(index)
    call_values = list(itertools.product(*call_values))
    
    data = []
    
    # Call functions with all parameter combinations
    for call_value in call_values:
        actual_call_value = list(call_value)
        for index, fixed_arg_pos in enumerate(fixed_arguments_pos):
            actual_call_value.insert(fixed_arg_pos, fixed_arguments[index])
        call_args = {key: value for key, value in zip(kwargs.keys(), actual_call_value)}
        data.append(function(**call_args))
        
    return data


def create_fixture(function, **kwargs):
    ''' Calls the function with the given kwargs values and stores the result.
    
    Numpy arrays are given as one parameter, lists parameters are looped with repeated
    function calls.
    '''

    # Check if all parameters are defined
    func_args = inspect.getargspec(function)[0]
    if not all([a in kwargs for a in func_args]):
        raise RuntimeError('Not all function arguments values defined')
    
    data = _call_function_with_args(function, **kwargs)
        
    # Store function return values in compressed pytable array
    data = np.array(data)
    with tb.open_file(os.path.join(FIXTURE_FOLDER, '%s.h5' % str(function.__name__)), 'w') as out_file:
        data_array = out_file.create_carray(out_file.root, name='Data', 
                                           title='%s return values' % function.__name__, atom=tb.Atom.from_dtype(data.dtype), 
                                           shape=data.shape, filters=tb.Filters(complib='blosc', complevel=5, fletcher32=False))
        data_array[:] = data


def check_with_fixture(function, **kwargs):
    ''' Calls the function with the given kwargs values and compares the result with the fixture.
    
    Numpy arrays are given as one parameter, lists parameters are looped with repeated
    function calls.
    '''

    with tb.open_file(os.path.join(FIXTURE_FOLDER, '%s.h5' % str(function.__name__)), 'r') as in_file:
        data_fixture = in_file.root.Data[:]
    
    data = np.array(_call_function_with_args(function, **kwargs))
    
    return np.allclose(data_fixture, data)
