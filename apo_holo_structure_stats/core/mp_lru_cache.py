# todo nefunguje - nehandluje hash kolize!!
import logging
import multiprocessing.dummy

logger = logging.getLogger(__name__)


class MultiProcessingLRUCache:
    """

    Dictionary keys are from Python 3.7 guaranteed to be ordered.
    """
    def __init__(self, fn, manager, max_size, args_to_key_fn=lambda x: x):
        self.fn = fn

        self.data = manager.dict()
        self.recently_used = manager.dict()  # small dict (to keep recently used keys), because emulating it with data dict would overload IPC

        self.max_size = max_size
        self.key_fn = args_to_key_fn

    def __call__(self, *args, **kwargs):
        key = tuple(map(self.key_fn, args)) + tuple(((self.key_fn(k), self.key_fn(v)) for k, v in kwargs.items()))

        try:
            res = self.data[key]
            logger.debug(f'CACHE HIT {key}, {list(self.data.keys())}')
        except KeyError:
            res = self.fn(*args, **kwargs)

            self.recently_used[key] = True
            self.data[key] = res

            logger.debug(f'CACHE MISS {key}, {list(self.data.keys())}')
        else:
            # move current key to top (was most recently used)
            # by deleting it and again adding the key (adding is done below)
            try:
                del self.recently_used[key]
            except KeyError:
                pass
            self.recently_used[key] = True

        # keep cache size close to max size
        # close, because we're not using any locks
        # hopefully it will work
        while len(self.data) > self.max_size:
            # the requirement for recently_used keys is that for every key in every moment, there must be the same key
            # in self.data. Delay between adding to data and then to recently_used is allowed. However may cause cache
            # misses not sure how many.
            try:
                oldest_key = next(iter(self.recently_used.keys()))
            except StopIteration:
                continue

            try:
                del self.data[oldest_key]
            except KeyError:
                pass
            try:
                del self.recently_used[oldest_key]  # tohle se přece smaže vždycky, když se smaže data?!! I kdyby někdo jinej, tak to smaže on. Ledaže by se mazaly ještě jinde?
                # nevim proc, ale v jenom tryi to nefungovalo, vzniknul deadlock, kterej chapu, ale nechapu, jak mohly nastat ty podminky pro něj (nesmazal se data..
                # mozna vim, v tom elsu nahore se to smazalo a pak vratilo.. Ale bylo to teda rychly
            except KeyError:
                pass

        return res

    # neni moc hezky, ale nechci menit a ted ani nevim jak tu serializaci v SerializableAnalyzer
    def __getattr__(self, item):
        if 'fn' not in vars(self):
            # for pickle, somehow `fn` is not present when (un)?pickling
            # so stop potential infinite recursion
            raise AttributeError

        return getattr(self.fn, item)
