from unittest import TestCase
import time
from scoary.progressbar import *


class Test(TestCase):
    def test_print_progress(self):
        start_time = datetime.now()
        n_tot = 20
        for i in range(n_tot + 1):
            time.sleep(0.05)
            msg = f'{i}: {" a" * i}'
            print_progress(i, n_tot, message=msg, start_time=start_time, message_width=30, default_width=120)

    def test_print_edge(self):
        start_time = datetime.now()
        msg = f'{0}: {" a" * 0}'
        print_progress(0, 0, message=msg, start_time=start_time, message_width=30, default_width=120)
