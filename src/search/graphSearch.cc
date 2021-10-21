using phmap::flat_hash_set;

class GraphSearch : public BaseSearch {
public:
  Graph* fg;
  int M;
  GraphSearch(string name, int k, int nSearch) : BaseSearch(name, k, nSearch) {
    fg = new Graph();
  }
  ~GraphSearch() {
    delete fg;
  }
  void param_register(vector<int> params) {
    cout << "[Parameters]\n--" << name << "\n";
    assert(params.size() == 1);   // there's only M
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
    return inner_product(data[id], current_query, D);
  }

  vector<operand> search(float* current_query){
    assert(ready == true);
    auto t1 = Clock::now();
    int seed = walk(current_query);
    auto t2 = Clock::now();
    expand(seed, current_query);
    auto t3 = Clock::now();
    timers[1] += compute_clock(t2, t1);
    timers[2] += compute_clock(t3, t2);
    return top_k_candidates;
  }
  int walk(float* current_query) {
    int maxCheck = 64;
    int curSeed = 0;
    float curSim = getSimilarity(curSeed, current_query);
    bool cont = true;
    while(cont) {
      cont = false;
      vector<int> &neighborList = (*fg->graph)[curSeed];
      for (size_t i = 0; i < min(maxCheck, (int)neighborList.size()); i++) {
        float neighborSim = getSimilarity(neighborList[i], current_query);
        if(neighborSim > curSim) {
          curSim = neighborSim;
          curSeed = neighborList[i];
          cont = true;
        }
      }
    }
    return curSeed;
  }

  void expand(int seed,  float *current_query) {
    auto t0 = Clock::now();
    int maxEdge = 64;
    int ef = 64;
    dbl_pq_d candidateSet;
    pq_d top_k;
    flat_hash_set<int> visited;
    visited.insert(seed);
    searched++;
    float sim = getSimilarity(seed, current_query);
    candidateSet.insert(make_pair(seed, sim));
    top_k_candidates.emplace_back(seed, sim);
    while(candidateSet.size() > 0 && searched < nSearch) {
      operand c = candidateSet.getMax();
      candidateSet.deleteMax();
      vector<int> &neighborList = (*fg->graph)[c.first];
      for (size_t i = 0; i < min(maxEdge, (int)neighborList.size()); i++) {
        int candidate = neighborList[i];
        if(visited.find(candidate) == visited.end()) {
          visited.insert(candidate);
          float furthestSimilarity = candidateSet.getMin().second;
          float neighborSimilarity = getSimilarity(candidate, current_query);
          searched++;
          if ((neighborSimilarity > furthestSimilarity) || (candidateSet.size() < ef)) {
            if (candidateSet.size() == ef) 
              candidateSet.deleteMin();
            candidateSet.insert(make_pair(candidate, neighborSimilarity));
            top_k_candidates.emplace_back(candidate, neighborSimilarity);
          }
        }
      }
    }
    auto t1 = Clock::now();
    partial_sort(top_k_candidates.begin(), top_k_candidates.begin() + K, top_k_candidates.end(), compare_descending);
    auto t2 = Clock::now();
    
    timers[3] += compute_clock(t1, t0);
    timers[4] += compute_clock(t2, t1);
  }
 
  void print_time() {
    cout << "[GraphSearch] " << (timers[1] + timers[2]) /queryCount << "(ms)\n";
    cout << "Preprocessing : " << timers[0] << "(ms)\n";
    cout << "Search " << (timers[1] + timers[2]) / queryCount << "(ms)\n";
    cout << "---[1] Walk : " << timers[1] / queryCount << "(ms)\n";
    cout << "---[2] Expand : " << timers[2] / queryCount << "(ms)\n";
    cout << "---[2-1] Expand-Main : " << timers[3] / queryCount << "(ms)\n";
    cout << "---[2-2] Expand-Sort : " << timers[4] / queryCount << "(ms)\n";
   }
  void print_information() {
    cout << "[GraphSearch] \n";
  }
}; 
