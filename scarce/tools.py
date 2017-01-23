r"""The tools module offers helper functions.

"""

import gzip
import dill
import logging
_LOGGER = logging.getLogger(__name__)


def save(obj, filename):
    ''' Save object to file'''
    _LOGGER.info('Saving to %s', filename)
    with gzip.open(filename, mode='wb') as out_file:
        dill.dump(obj, out_file, protocol=0)


def load(filename):
    ''' Load object from file'''
    _LOGGER.info('Loading %s', filename)
    with gzip.open(filename, mode='rb') as in_file:
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
