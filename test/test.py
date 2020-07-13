import pandas as pd
from Levenshtein import distance
import time

d = pd.read_csv("levenshteinTests.csv").values


def f():
    for rw in d:
        distance(rw[0], rw[1]) == rw[2]

t = time.time()
f()
print(time.time()-t)
