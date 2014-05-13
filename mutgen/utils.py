"Utility functions"

import itertools as it
from operator import itemgetter

def iunzip(iterable):
    """Iunzip is the same as zip(*iter) but returns iterators, instead of
    expand the iterator. Mostly used for large sequence"""
    _tmp, iterable = it.tee(iterable, 2)
    iters = it.tee(iterable, len(_tmp.next()))
    return (it.imap(itemgetter(i), itbl) for i, itbl in enumerate(iters))


