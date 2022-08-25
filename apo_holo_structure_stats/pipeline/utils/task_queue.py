import itertools
from collections import deque
from concurrent.futures import Executor, ProcessPoolExecutor


# podle me by chunksize mel efekt i u threadpoolu, protoze je min futures i _workItems

# nebo použít chunksize?? - to by zmensilo mozna i memory usage, hlavne ale overhead komunikace ne?
# sem pridat chunksize, pro processpool?
from functools import partial
from traceback import format_exc


def submit_tasks(executor: Executor, window_size: int, fn, *iterables):
    """ Assuming all task take approximately the same time. Effectively executor.map, but done in batches/windows to
    reduce the number of queued Future objects (seems to consume a lot of memory). """

    args_iterator = zip(*iterables)

    # initially submit `window_size` tasks, and put them into the queue
    fs_queue = deque(
        (executor.submit(fn, *args) for args in itertools.islice(args_iterator, window_size)),
        maxlen=window_size,
    )

    # pop tasks and for each completed add a new one
    try:
        while fs_queue:
            # following comment copied from Executor.map
            # Careful not to keep a reference to the popped future
            yield fs_queue.popleft()
            fs_queue.append(executor.submit(fn, *next(args_iterator)))
    except StopIteration:
        # no more tasks to submit
        pass

    # collect remaining tasks (after all tasks have been submitted)
    while fs_queue:
        # Careful not to keep a reference to the popped future
        yield fs_queue.popleft()


""" Following three functions extracted from stdlib concurrent.futures.process. """
def _get_chunks(*iterables, chunksize):
    """ Iterates over zip()ed iterables in chunks. """
    it = zip(*iterables)
    while True:
        chunk = tuple(itertools.islice(it, chunksize))
        if not chunk:
            return
        yield chunk


def _process_chunk(fn, chunk):
    """ Processes a chunk of an iterable passed to map.

    Runs the function passed to map() on a chunk of the
    iterable passed to map.

    This function is run in a separate process.

    """
    return [fn(*args) for args in chunk]


def _chain_from_iterable_of_lists(iterable):
    """
    Specialized implementation of itertools.chain.from_iterable.
    Each item in *iterable* should be a list.  This function is
    careful not to keep references to yielded objects.
    """
    for element in iterable:
        element.reverse()
        while element:
            yield element.pop()

""" Exception and handling are like in Pebble library. """

class RemoteTraceback(Exception):
    """Traceback wrapper for exceptions in remote process.
    Exception.__cause__ requires a BaseException subclass.
    """

    def __init__(self, traceback):
        self.traceback = traceback

    def __str__(self):
        return self.traceback


class RemoteExceptionVehicle(object):
    """Pickling wrapper for exceptions in remote process."""

    def __init__(self, exception, traceback):
        self.exception = exception
        self.traceback = traceback

    def __reduce__(self):
        return rebuild_exception, (self.exception, self.traceback)


def rebuild_exception(exception, traceback):
    exception.__cause__ = RemoteTraceback(traceback)

    return exception


def process_execute(function, *args, **kwargs):
    """Runs the given function returning its results or exception."""
    try:
        return function(*args, **kwargs)
    except Exception as error:
        error.traceback = format_exc()
        return RemoteExceptionVehicle(error, error.traceback)


class FutureLike:
    def __init__(self, result_or_exception):
        self.result_or_exception = result_or_exception

    def result(self):
        if isinstance(self.result_or_exception, Exception):
            raise self.result_or_exception
        # for compatibility with ThreadPoolExecutor (for testing)
        elif isinstance(self.result_or_exception, RemoteExceptionVehicle):
            raise rebuild_exception(self.result_or_exception.exception, self.result_or_exception.traceback)

        return self.result_or_exception


def submit_short_tasks(executor: ProcessPoolExecutor, window_size: int, p_chunk_size: int, fn, *iterables):
    # wrap fn in exception-handling fn
    fn = partial(process_execute, fn)

    chunk_futures = submit_tasks(executor,
                                 window_size,
                                 partial(_process_chunk, fn),  # 'vectorify' the function to return a list of results
                                 _get_chunks(*iterables, chunksize=p_chunk_size))

    results_or_exceptions = (chunk_f.result() for chunk_f in chunk_futures)
    return (FutureLike(res_or_exc) for res_or_exc in _chain_from_iterable_of_lists(results_or_exceptions))

