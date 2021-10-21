import numpy as np
import math
import sys
from random import randint

def vnball(n) :
  pipow = math.pi**(n/2.0)
  gammapow = math.gamma(n/2.0+1)
  return pipow/gammapow

def vnsphere(n) :
  return (n+1) * vnball(n+1)

def gammafunction(n) :
  return (math.gamma(n/2.0) / math.gamma((n-1)/2.0)) * (math.gamma(n/2.0) / math.gamma((n+1)/2.0))


def essReference(ver, d, maxdim, mindim=1) :
  if (d==1) :
    n = range(max(mindim, 2), maxdim+1)
    dimval = np.array([gammafunction(i) for i in n])
  else :
    print "d is not 1 : wrong"
  
  if (mindim <= d) :
    dimval = np.append(np.zeros(d-mindim+1), dimval)
  return dimval

def gencoordinates(m, n):
  seen= set()
  x,y = randint(m,n), randint(m,n)
  while len(seen) < ((n+1-m)*(n-m))/2:
    while ((x,y) in seen) or ((y,x) in seen) :
      x,y = randint(m,n), randint(m,n)
      while (x==y) :
        x,y = randint(m,n), randint(m,n)
    seen.add((x,y))
    yield (x,y)
  return


def computeEss(data, ver, d=1) :
  p = d+1 #2
  n = data.shape[1] #D
  if(p>n): return 0

  vectors = data - data.mean(axis=0, keepdims=True)
  N = vectors.shape[0]
  sampler = gencoordinates(0, N-1)

  weight = 0
  vol = 0
  for it in range(3000) :
    sample = next(sampler)
    i,j = sample[0], sample[1]
    weight += np.linalg.norm(vectors[i]) * np.linalg.norm(vectors[j])
    mat = np.vstack((vectors[i], vectors[j]))
    mat_t = np.transpose(mat)
    vol += np.sqrt(np.linalg.det(np.dot(mat, mat_t)))

  return vol/weight

def essLocalDimEst(data, ver='a', d=1) :
  essval = computeEss(data, ver, d)
  mindim=1
  maxdim=20
  dimvals = essReference(ver, d, maxdim, mindim)
  while (essval > dimvals[maxdim-1]) :
    mindim = maxdim+1
    maxdim = 2*(maxdim-1)
    dimvals = np.append(dimvals, essReference(ver, d, maxdim, mindim))
  i = np.where(dimvals[mindim:maxdim+1]>=essval)[0][0]
  de_integer = mindim+i
  if not (dimvals[de_integer-1]<essval and dimvals[de_integer]>=essval):
    print "min/maxdim : " + str(mindim) + " / " + str(maxdim)
    print "i : " + str(i)
    print "essval : " + str(essval)
    print "dimvals[de_integer-1] : " + str(dimvals[de_integer-1])
    print "dimvals[de_integer] : " + str(dimvals[de_integer])
    assert false

  de_fractional = (essval-dimvals[de_integer])/(dimvals[de_integer+1]-dimvals[de_integer])
  de = de_integer + de_fractional
  return de
 
#data = np.loadtxt('20000_100_normal.txt')
#data = np.loadtxt('movielens.txt')
#data = np.loadtxt('netflix_movie.txt')
#data = np.loadtxt('yelp.txt')
#data = np.loadtxt('music.txt')
try :
  data = np.loadtxt(sys.argv[1])
  print "Intrinsic Dimension : " + str(essLocalDimEst(data))
except ValueError:
  print("Please remove the first line of the data file - N D") 

 

