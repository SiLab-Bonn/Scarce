import inspect
import os
import tables as tb
import numpy as np
import itertools

from scarce import constant


def compare_arrays(array_1, array_2, threshold=0.01):
    ''' Compares two arrays and returns the number of
        data points in percent where the value difference
        in any dimension between the two arrays exceeds
        a threshold.
        A threshold value of 0.01 corresponds to a relative
        error of 1% given by the data of array_1.
    '''

    if array_1.shape != array_2.shape:
        raise ValueError('The shape of the arrays does not match.')

    bad_selection = np.abs(array_1[0] - array_2[0]) > np.abs(array_1[0] * threshold)

    for index in range(1, len(bad_selection.shape)):
        bad_selection = np.logical_or(bad_selection, np.abs(array_1[index] - array_2[index]) > np.abs(array_1[index] * threshold))

    bad_points = np.where(bad_selection)
    n_points = array_1.ravel().shape[0] / 2.

    return float(bad_points[0].shape[0]) / n_points


def _call_function_with_args(function, **kwargs):
    ''' Calls the function with the given kwargs
    and returns the result in a numpy array. All combinations
    of functions arguments in a list are used for multiple
    function calls.'''

    # Create all combinations of arguments from list parameters
    # This is ugly but avoids recursion and does effectively
    # a nested loop of n parameters:
    # for par_1 in pars_1:
    #  for par_2 in pars_2:
    #    ...
    #    for par_n in pars_n:
    #      function(par_1, par_2, ..., par_n)

    call_values = []  # Arguments with permutations
    fixed_arguments = []  # Constant arguments
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
        call_args = {
            key: value for key, value in zip(kwargs.keys(), actual_call_value)}
        data.append(function(**call_args))

    return np.array(data)


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
    with tb.open_file(os.path.join(constant.FIXTURE_FOLDER, '%s.h5' % str(function.__name__)), 'w') as out_file:
        data_array = out_file.create_carray(out_file.root, name='Data',
                                            title='%s return values' % function.__name__, atom=tb.Atom.from_dtype(data.dtype),
                                            shape=data.shape, filters=tb.Filters(complib='blosc', complevel=5, fletcher32=False))
        data_array[:] = data


def check_with_fixture(function, **kwargs):
    ''' Calls the function with the given kwargs values and compares the result with the fixture.

    Numpy arrays are given as one parameter, lists parameters are looped with repeated
    function calls.
    '''

    with tb.open_file(os.path.join(constant.FIXTURE_FOLDER, '%s.h5' % str(function.__name__)), 'r') as in_file:
        data_fixture = in_file.root.Data[:]

    data = _call_function_with_args(function, **kwargs)

    return np.allclose(data_fixture, data)
