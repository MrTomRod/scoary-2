import multiprocessing as mp


def init():
    mgr = mp.Manager()
    ns = mgr.Namespace()
    counter = mgr.Value('i', 0)
    lock = mgr.Lock()
    return mgr, ns, counter, lock
