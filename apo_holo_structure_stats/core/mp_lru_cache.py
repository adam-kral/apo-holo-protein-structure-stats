
class MultiProcessingLRUCache:
    """

    Dictionary keys are from Python 3.7 guaranteed to be ordered.
    """
    def __init__(self, fn, manager, max_size):
        self.fn = fn
        self.d = manager.dict()
        self.max_size = max_size

    def __call__(self, *args, **kwargs):
        hashed_args = hash(args + tuple(kwargs.items()))

        try:
            res = self.d[hashed_args]
            # move current key to top (was LRU)
            # by deleting it and again adding the key (adding is done below)
            del self.d[hashed_args]
        except KeyError:
            res = self.fn(*args, **kwargs)

        self.d[hashed_args] = res

        # keep cache size close to max size
        # close, because we're not using any locks
        # hopefully it will work
        while len(self.d) > self.max_size:
            try:
                oldest_key = next(self.d.keys())
                del self.d[oldest_key]
            except StopIteration:
                pass

        return res
