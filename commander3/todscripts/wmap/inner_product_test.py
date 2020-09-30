from multiprocessing import Process, Value, Array
import numpy as np


def worker1(arr):
    arr[0] += 1


def worker2(arr):
    arr[1] += 1

array = Array("i", [0,1])

p1 = Process(target=worker1, args=[array])
p2 = Process(target=worker2, args=[array])

p1.start()
p2.start()
p1.join()
p2.join()

print(array[:])
