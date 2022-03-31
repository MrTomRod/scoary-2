import multiprocessing as mp

mgr = mp.Manager()
ns = mgr.Namespace()
counter = mgr.Value('i', 0)
lock = mgr.Lock()
