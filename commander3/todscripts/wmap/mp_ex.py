from multiprocessing import Pool

class Adder:
    """I'm using this class in place of a monte carlo simulator"""

    def add(self, a, b):
        return a + b

def setup(x, y, z):
    """Sets up the worker processes of the pool. 
    Here, x, y, and z would be your global settings. They are only included
    as an example of how to pass args to setup. In this program they would
    be "some arg", "another" and 2
    """
    global adder
    adder = Adder()

def job(a, b):
    """wrapper function to start the job in the child process"""
    return adder.add(a, b)

if __name__ == "__main__":
    num = 10000000
    args = list(zip(range(num), range(num, 2*num)))
    # args == [(0, 10), (1, 11), ..., (8, 18), (9, 19)]

    with Pool(initializer=setup, initargs=["some arg", "another", 2]) as pool:
        # runs jobs in parallel and returns when all are complete
        results = pool.starmap(job, args)

    #print(results) # prints [10, 12, ..., 26, 28] 
