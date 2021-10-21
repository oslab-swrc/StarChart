// SPDX-FileCopyrightText: Copyright Copyright 2021 Seoul National University
//
// SPDX-License-Identifier: MIT License

//#define ORIGINALIMP
#ifdef ORIGINALIMP
class StreamTracker {
public:
  float* tracker;
  int* ids;
  int maxSize;
  int size;
  float smallest;
  int smallestID;
  flat_hash_set<int> trackerSet;

  StreamTracker(int n) {
    tracker = new float[n];
    ids = new int[n];
    for(int i=0; i<n; i++) {
      tracker[i] = numeric_limits<float>::lowest();
      ids[i] = -1;
    }
    maxSize = n;
    smallestID = 0;
    smallest = numeric_limits<float>::lowest();
    size = 0;
  }
  ~StreamTracker() {
    delete[] tracker;
    delete[] ids;
  }
  void update(int id, float sim) {
    if(size == maxSize && sim <= smallest) {
      return;
    }
    else if(size == maxSize && sim > smallest) {
      if(trackerSet.find(id) != trackerSet.end())
        return;
      trackerSet.erase(ids[smallestID]);
      ids[smallestID] = id;
      tracker[smallestID] = sim;
      trackerSet.insert(id);
    }
    else if(size < maxSize) {
      if(trackerSet.find(id) != trackerSet.end())
        return;
      ids[smallestID] = id;
      tracker[smallestID] = sim;
      size++;
      trackerSet.insert(id);
    }
    float t = tracker[0];
    int ti = 0;
    for(int i = 0; i < maxSize; i++) {
      if(tracker[i] < t) {
        t = tracker[i];
        ti = i;
      }
    }
    smallest = t;
    smallestID = ti;
  }
  operand get() {
    float l = tracker[0];
    int li = 0;
    assert(size > 0);
    for(int i = 0; i < maxSize; i++) {
      if(tracker[i] > l) {
        l = tracker[i];
        li = i;
      }
    }
    assert(ids[li] != -1);
    operand ret = make_pair(ids[li], l);
    trackerSet.erase(ids[li]);
    tracker[li] = numeric_limits<float>::lowest();
    ids[li] = -1;
    smallestID = li;
    smallest = numeric_limits<float>::lowest();
    size--;
    return ret;
  }
  

  // bool isLargest(float sim) {
  //   for(int i = 0; i < maxSize; i++) {
  //     if(tracker[i] > sim)
  //       return false;
  //   }
  //   return true;
  // }
  // operand check() {
  //   float l = tracker[0];
  //   int li = 0;
  //   assert(size > 0);
  //   for(int i = 0; i < maxSize; i++) {
  //     if(tracker[i] > l) {
  //       l = tracker[i];
  //       li = i;
  //     }
  //   }
  //   operand ret = make_pair(ids[li], l);
  //   return ret;
  // }
  void deleteItem(int id) {
    for(int i = 0; i < maxSize; i++) {
      if(ids[i] == id) {
        smallestID = i;
        smallest = numeric_limits<float>::lowest();
        trackerSet.erase(ids[i]);
        ids[i] = -1;
        tracker[i] = numeric_limits<float>::lowest();
        size--;
        break;
      }
    }
    assert(size >=0);
  }
};
#else
#include <set>

struct compare_descending_seperate {
  bool operator() (const operand& lhs, const operand& rhs) const {
    if (lhs.second > rhs.second && lhs.first != rhs.first) 
      return true;
    else
      return false;
  }
};

class StreamTracker {
public:
  set<operand, compare_descending_seperate> s;
  int size;
  StreamTracker(int n) {
    size = n;
  }
  ~StreamTracker() {
  }
  void update(int id, float sim) {
    if(s.size() == size && s.rbegin()->second >= sim) {
      return;
    }
    else if(s.size() == size && s.rbegin()->second < sim) {
      if(s.find(make_pair(id, sim)) != s.end())
        return;
      s.erase(prev(s.end()));
      s.insert(make_pair(id, sim));
    }
    else if(s.size() < size) {
      if(s.find(make_pair(id, sim)) != s.end())
        return;
      s.insert(make_pair(id, sim));
    }
  }
  operand get() {
    operand ret = *(s.begin());
    s.erase(s.begin());
    return ret;
  }
};
#endif
