"""
parmap.py -- Utility function that does a parallel Map (map/collect) of a Function onto a list of Items
Degree of parallelism is controlled by numPartitions parameter.
Parallel framework is chosen by mode parameter.
"""

import sys
import pysparkling, multiprocessing, cloudpickle, pickle
from parmap_framework.parmap_utils import Timer, warn

# Legal parallel execution modes
MODES = [
    "seq",  # sequential
    "par",  # pool.map using multiprocessing library
    "sparkling",  # spark map/collect using pysparkling library (multicore only)
    "spark",  # spark map/collect using spark cluster
    "dask"    # map/collect using dask distributed scheduler
]

def warn(s):
    print(s, file=sys.stderr)

# def warn(s): print >>sys.stderr, s

class Parmap:
    """Do a parallel Map of a function onto a list of items and then collect and return the list of results.
    The function receives a list of values or S3 URL's and can perform I/O as a side effect (e.g. write files/objects and return their URLs).
    """

    def __init__(
        self,
        mode="par",  # one of the six modes above
        numWorkers=None,  # number of workers, degree of parallelism (None means auto)
        master=None,  # URL for master scheduler
        config={},  # configuration dict for additional info:  credentials, etc.
        context=None,  # reuse existing scheduler context:  for spark or sparkling, a Context object; for dask, a Client object
    ):
        if mode not in MODES:
            warn(
                'parmap: Bad mode arg, using "par" (local multicore) instead: %s' % mode
            )
            mode = "par"
        self.mode = mode
        self.numWorkers = numWorkers
        self.master = master
        self.config = config
        self.context = context

        if mode == "seq":
            pass

        elif mode == "par":
            if self.context is None:
                worker = config.get("worker", None)
                if worker is not None and worker == "thread":
                    self.context = multiprocessing.pool.ThreadPool(numWorkers)
                else:
                    self.context = multiprocessing.Pool(numWorkers)

        elif mode == "sparkling":
            if self.context is None:
                try:
                    worker = config.get("worker", None)
                    if worker is not None and worker == "thread":
                        self.context = pysparkling.Context(
                            multiprocessing.pool.ThreadPool(numWorkers),
                            serializer=cloudpickle.dumps,
                            deserializer=pickle.loads,
                        )
                    else:
                        self.context = pysparkling.Context(
                            multiprocessing.Pool(numWorkers),
                            serializer=cloudpickle.dumps,
                            deserializer=pickle.loads,
                        )
                except:
                    raise Exception("parmap: mode sparkling, cannot get Context object")

        elif mode == "spark":
            import pyspark

            if type(context) == type(pyspark.SparkContext):
                return
            if self.context is None:
                try:
                    appName = config["appName"]
                except:
                    appName = None
                try:
                    self.context = pyspark.SparkContext(master=master, appName=appName)
                except:
                    raise Exception("parmap; mode spark, cannot get Context object")

        elif mode == "dask":
            from dask.distributed import Client, LocalCluster

            if type(context) == type(Client):
                return
            if master is None:
                master = config.get("master", None)
            if master is not None:
                try:
                    self.context = Client(
                        master,
                        serializers=["dask", "cloudpickle"],
                        deserializers=["dask", "pickle"],
                    )
                    self.context.upload_file("/home/bdwilson/code/AIST/framework/parmap.py")
                    warn("parmap: mode dask, using cluster scheduler at: %s" % master)
                except:
                    raise Exception(
                        "parmap: mode dask, could not open cluster client at: %s"
                        % master
                    )
            else:
                #                self.context = Client( LocalCluster(numWorkers) )    # makes default LocalCluster using numWorkers
                self.context = (
                    Client()
                )  # makes default LocalCluster using number of cores
                warn("parmap: mode dask, no client URL so using LocalCluster()")

    def __call__(
        self,
        fn,
        items,
        numPartitions=None,
    ):
        if numPartitions is None:
            numPartitions = self.numWorkers
        mode = self.mode
        ctx = self.context
        warn(
            "\nparmap %s: Running in mode %s with numPartitions %s"
            % (str(fn), mode, str(numPartitions))
        )

        if mode == "seq":
            return list(map(fn, items))

        elif mode == "par":
            return ctx.map(fn, items)  # here a Pool of processes

        elif mode == "sparkling" or mode == "spark":
            itemRDD = ctx.parallelize(items, numPartitions)
            return itemRDD.map(fn).collect()

        elif mode == "dask":
            futures = ctx.map(fn, items, key="parmap", retries=1)
            return ctx.gather(futures)

def test(fn, n):
    from math import pow

    items = list(range(1, n))
    for mode in ["seq", "par", "sparkling", "dask"]:
        parmap = Parmap(mode, numWorkers=4, config={"appName": "parmap_factorial_test"})
        with Timer():
            results = parmap(fn, items, numPartitions=4)
        print(results[-1])


def cubed(i):
    from math import pow

    return pow(i, 3)


if __name__ == "__main__":
    from math import factorial

    fn = sys.argv[1]
    n = int(sys.argv[2])
    if fn == "cubed":
        fn = cubed
    else:
        fn = factorial
    test(fn, n)
