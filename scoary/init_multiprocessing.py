import multiprocessing as mp

mp.set_start_method('spawn')


def init():
    mgr = mp.Manager()
    ns = mgr.Namespace()
    counter = mgr.Value('i', 0)
    lock = mgr.Lock()
    return mgr, ns, counter, lock
