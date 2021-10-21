import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import itertools

name = ["Hash", "RangeHash", "Graph", "HashGraph", "HNSW", "CosineGraph"]
dataset = sys.argv[1]
topK = sys.argv[2]

#recall-time graph
for name_ in name :
  data = np.loadtxt(dataset + "_" + name_+topK+".txt", dtype = 'float32')
  time_ = data[:,1]
  recall_ = data[:,0]

  duplist = list(zip(time_, recall_))
  recall, time = np.array([]), np.array([])
  for j in sorted(set(recall_)) :
    mintime = min([t for t,r in duplist if r==j])
    time = np.append(time, mintime)
    recall = np.append(recall, j)
  plt.plot(recall, time)
plt.xlabel('Recall(%)')
plt.ylabel('Time(ms)')
plt.legend(name, loc='upper left')
plt.show()

#nSearch-recall graph
for name_ in name :
  data = np.loadtxt(dataset + "_" + name_+topK+".txt", dtype = 'float32')
  recall_ = data[:,0]
  nsearch_ = data[:,2]
  #recall-time graph
  
  duplist = list(zip(recall_, nsearch_))
  recall, nsearch = np.array([]), np.array([])
  for j in sorted(set(nsearch_)) :
    maxrecall = max([r for r,n in duplist if n==j])
    recall = np.append(recall, maxrecall)
    nsearch = np.append(nsearch, j)
  plt.plot(nsearch, recall)
plt.xlabel('#nSearch')
plt.ylabel('Recall(%)')
plt.legend(name, loc='upper left')
plt.show()
