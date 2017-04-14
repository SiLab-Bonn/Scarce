r"""The tools module offers helper functions.

"""

import dill
import gzip
import logging

import numpy as np
from scipy import interpolate
_LOGGER = logging.getLogger(__name__)


def save(obj, filename, compress=True):
    ''' Save object to file'''
    _LOGGER.info('Saving to %s', filename)
    if compress:
        with gzip.open(filename, mode='wb') as out_file:
            dill.dump(obj, out_file, protocol=0)
    else:
        with open(filename, mode='wb') as out_file:
            dill.dump(obj, out_file, protocol=0)


def load(filename, compress=True):
    ''' Load object from file'''
    _LOGGER.info('Loading %s', filename)
    if compress:
        with gzip.open(filename, mode='rb') as in_file:
            return dill.load(in_file)
    else:
        with open(filename, mode='rb') as in_file:
            return dill.load(in_file)


def apply_async(pool, fun, args=None, **kwargs):
    ''' Run fun(*args, **kwargs) in different process.

    fun can be a complex function since pickling is not done with the
    cpickle module as multiprocessing.apply_async would do, but with
    the more powerfull dill serialization.
    Additionally kwargs can be given and args can be given'''
    payload = dill.dumps((fun, args, kwargs))
    return pool.apply_async(_run_with_dill, (payload,))


def _run_with_dill(payload):
    ''' Unpickle payload with dill.

    The payload is the function plus arguments and keyword arguments.
    '''
    fun, args, kwargs = dill.loads(payload)
    if args:
        return fun(*args, **kwargs)
    else:
        return fun(**kwargs)


def time_data_interpolate(T, data, t, axis=-1, fill_value=0., kind='linear'):
    ''' Interpolates data recorded at time points T to time point t.
    '''
    results = np.full(shape=(t.shape[0], data.shape[1]), fill_value=np.nan)

    for i in range(data.shape[1]):  # Loop over e-h pairs
        f = interpolate.interp1d(T[:, i], data[:, i], axis=axis,
                                 bounds_error=False, fill_value=fill_value,
                                 kind=kind)
        results[:, i] = f(t)
    return results
    
