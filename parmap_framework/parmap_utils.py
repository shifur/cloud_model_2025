import time
import sys


def warn(s):
    print(s, file=sys.stderr)


class Timer:
    def __enter__(self):
        self.start = time.perf_counter()
        return self

    def __exit__(self, *args):
        self.end = time.perf_counter()
        self.interval = self.end - self.start
        warn("timer (sec): " + str(self.interval))