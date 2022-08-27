import concurrent
import itertools
from collections import deque
from concurrent.futures import Executor, ProcessPoolExecutor, Future

from functools import partial
from traceback import format_exc


class SemiBlockingQueueExecutor(Executor):
    """ Saves memory when thousands of tasks are submitted to ProcessPoolExecutor or ThreadPoolExecutor.

    In ThreadPoolExecutor and ProcessPoolExecutor internally created WorkItem and Future objects which take up memory.
    Block `submit` when too many futures are pending."""
    def __init__(self, executor: Executor, max_pending_futures: int):
        self.executor = executor
        self.max_pending_futures = max_pending_futures
        self.pending_futures = set()

    def submit(self, fn, *args, **kwargs) -> Future:
        # Enforce upper bound on number of internally created WorkItem adn Future objects which take up memory
        # block until a free slot
        if len(self.pending_futures) >= self.max_pending_futures:
            done, self.pending_futures = concurrent.futures.wait(self.pending_futures,
                                                                 return_when=concurrent.futures.FIRST_COMPLETED)
        future = self.executor.submit(fn, *args, **kwargs)
        self.pending_futures.add(future)
        return future


def submit_tasks(executor: Executor, window_size: int, fn, *iterables):
    """ Assuming all task take approximately the same time. Effectively executor.map, but done in batches/windows to
    reduce the number of queued Future objects (seems to consume a lot of memory). """

    queue_executor = SemiBlockingQueueExecutor(executor, window_size)
    yield from queue_executor.map(fn, *iterables)


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

