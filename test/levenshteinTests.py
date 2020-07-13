import pandas as pd
import numpy as np
import Levenshtein
import random

random.seed(12345)

d = pd.read_csv("../data/asjp19wide.csv", index_col=0)


words = d.values[~d.isnull()]
words = np.concatenate([w.split('-') for w in words])

tests = pd.DataFrame(columns=['word1', 'word2', 'LD'])


for i in range(1000):
    if i % 100 == 0:
        print(i)
    w1, w2 = random.sample(list(words), 2)
    tests.loc[i] = [w1, w2, Levenshtein.distance(w1, w2)]


tests.to_csv('levenshteinTests.csv', index=False)
