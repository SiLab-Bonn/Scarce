r"""The tools module offers helper functions.

"""

import gzip
import dill as pickle
import logging
_LOGGER = logging.getLogger(__name__)


def save(obj, filename):
    _LOGGER.info('Saving to %s', filename)
    with gzip.open(filename, mode='wb') as out_file:
        pickle.dump(obj, out_file, protocol=0)


def load(filename):
    _LOGGER.info('Loading %s', filename)
    with gzip.open(filename, mode='rb') as in_file:
        return pickle.load(in_file)
