using phmap::flat_hash_map;
using phmap::flat_hash_set;
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
  bool isLargest(float sim) {
    for(int i=0; i<maxSize; i++) {
      if(tracker[i] > sim)
        return false;
    }
    return true;
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
    for(int i=0; i<maxSize; i++) {
      if(tracker[i] < t) {
        t = tracker[i];
        ti = i;
      }
    }
    smallest = t;
    smallestID = ti;
  }
  operand check() {
    float l = tracker[0];
    int li = 0;
    assert(size >0);
    for(int i=0; i<maxSize; i++) {
      if(tracker[i] > l) {
        l = tracker[i];
        li = i;
      }
    }
    operand ret = make_pair(ids[li], l);
    return ret;
  }
  operand get() {
    float l = tracker[0];
    int li = 0;
    assert(size >0);
    for(int i=0; i<maxSize; i++) {
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

    //cout << ids[li] << " li " << li << " " << tracker[li] << "\n";
    //cout << ret.first << "rf " << ret.second << "rs \n";
    return ret;
  }
  void deleteItem(int id) {
    for(int i=0; i<maxSize; i++) {
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
class GraphSearchTJ : public BaseSearch {
public:
  Graph* fg;
  float **dirs;
  int M;

  GraphSearchTJ(string name, int k, int nSearch) : BaseSearch(name, k, nSearch) {
    fg = new Graph();
  }
  ~GraphSearchTJ() {
  }
  void param_register(vector<int> params) {
    cout << "[Parameters]\n--" << name << "\n";
    assert(params.size() == 1);
    for (size_t i = 0; i < params.size(); i++) {
      if (i == 0) {
        M = params[i];
        cout << "M : " << M << "\n";
      }
    }
  }
  void preprocess(){
    auto t1 = Clock::now();
    auto t2 = Clock::now();
    timers[0] += compute_clock(t2, t1);
    ready = true;
  }
  float inline getSimilarity(int id, float *current_query) {
    float sim = inner_product(data[id], current_query, D);
    return sim;
  }
  vector<operand> search(float* current_query){
    vector<operand> answer;
    vector<int> rank;
    rank.reserve(N);
    answer.reserve(N);
    for(int i=0; i<N; i++) {
      float sim = getSimilarity(i, current_query);
      answer.push_back(make_pair(i, sim));  
      rank.push_back(-1);
    }
    sort(answer.begin(), answer.end(), compare_descending);
    for(int i=0; i<N; i++) {
      rank[answer[i].first] = i;
    }

    assert(ready == true);
    flat_hash_set<int> visited;
    flat_hash_set<int> inserted;
    inserted.reserve(4096);
    visited.reserve(512);
    StreamTracker st(4);
    auto t1 = Clock::now();
    //int init = rand() % N;
    int init = 0;
    st.update(init, getSimilarity(init, current_query));
    int walkcount = 0;
    for(int i = 0; i < 128; i++) {
      walk(init, visited, inserted, current_query, rank, st);
      //walk(init, visited, inserted, current_query);
      if(searched > nSearch) {
        walkcount = i;
        break;
      }
    }
    auto t2 = Clock::now();
    partial_sort(top_k_candidates.begin(), top_k_candidates.begin() + K, top_k_candidates.end(), compare_descending);
    auto t3 = Clock::now();
    cout << "Walk : " << walkcount << "\n";
    timers[1] += compute_clock(t2, t1);
    timers[2] += compute_clock(t3, t2);
    return top_k_candidates;
  }

  void inline insertCandidate(flat_hash_set<int> &visited, int id, float sim) { ///, float* current_query) {
    if(visited.find(id) == visited.end()) {
      visited.insert(id);
      top_k_candidates.emplace_back(id, sim);
    }
  }

  void walk(int &seed, flat_hash_set<int> &visited, flat_hash_set<int> &inserted, float *current_query, vector<int> &rank, StreamTracker &st) {
    operand cur = st.get();
    float curSim = cur.second;//getSimilarity(seed, current_query);
    int curSeed = cur.first;//seed;
    //int curSeed = seed;
    //float curSim  = getSimilarity(seed, current_query);
    assert(curSeed != -1);
    bool found = true;
    int searchEdge = 1024;
    int obufSize = 0;

    while(found) {
      visited.insert(curSeed);
      st.deleteItem(curSeed);
      insertCandidate(inserted, curSeed, curSim);
      
      //if(wc == 0)
      //   cout << "\n";
      //cout << "CurSim: " << curSim << " / " << wc << " /  Rank : " << rank[curSeed] << "\n";
      
      // vector<int> &printneighbor = fg->graph[curSeed];
      // for (int i=0; i<searchEdge; i++) {
      //   if(rank[printneighbor[i]] <= 50)
      //   cout << i << "th edge : " << rank[printneighbor[i]] << " / " << (bool)(visited.find(printneighbor[i]) != visited.end()) << "\n";
      // }
      // wc++;
      found = false;
      vector<operand> obuf;
      obuf.reserve(obufSize);
      vector<int> &neighborList = (*fg->graph)[curSeed];
      float neighborMaxSim = numeric_limits<float>::lowest();
      int maxNeighbor = -1;

      for (int i=0; i<searchEdge; i++) {
        int candidate = neighborList[i];
        float neighborSim = getSimilarity(candidate, current_query);
        searched++;
        if(i < obufSize)
          obuf.emplace_back(candidate, neighborSim);
        if(visited.find(candidate) == visited.end()) {
          st.update(candidate, neighborSim);          
        }
        if(neighborSim > curSim && visited.find(candidate) == visited.end()) {
          curSim = neighborSim;
          curSeed = candidate;
          found = true;
        }
        if(neighborSim > neighborMaxSim && visited.find(candidate) == visited.end()) {
          neighborMaxSim = neighborSim;
          maxNeighbor = candidate;
        }
        if(found && i >= 8) {
          //cout << "Move at " << i << "\n";
          break;
        }
      }
      if(!found) {
        seed = maxNeighbor;
        for(int i=0; i<obufSize; i++) {
          insertCandidate(inserted, obuf[i].first, obuf[i].second);
        }
      }
    }
  }
  void print_time() {
    cout << "[GraphSearchTJ] " << (timers[1] + timers[2]+timers[3] + timers[4]) /queryCount << "(ms)\n";
    cout << "Preprocessing : " << timers[0] << "(ms)\n";
    cout << "Search " << (timers[1] + timers[2]) / queryCount << "(ms)\n";
    cout << "---[1] Walk Front : " << timers[1] / queryCount << "(ms)\n";
    cout << "---[2] Walk : " << timers[2] / queryCount << "(ms)\n";
    cout << "---[3] Expand : " << timers[3] / queryCount << "(ms)\n";
    cout << "---[4] Expand: " << timers[4]/queryCount << "(ms)\n";
   }
  void print_information() {
    cout << "[GraphSearchTJ] \n";
  }
}; 