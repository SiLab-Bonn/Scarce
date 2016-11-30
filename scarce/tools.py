r"""The tools module offers helper functions.

"""

import gzip
import dill as pickle


def save(obj, filename):
    with gzip.open(filename, mode='wb') as out_file:
        pickle.dump(obj, out_file, protocol=0)


def load(filename):
    with gzip.open(filename, mode='rb') as in_file:
        return pickle.load(in_file)
