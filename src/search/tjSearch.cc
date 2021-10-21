using phmap::flat_hash_map;
using phmap::flat_hash_set;
class GraphSearchTJBase : public BaseSearch {
public:
  Graph* fg;
  int M;
  int lastSeed;
  GraphSearchTJBase(string name, int k, int nSearch) : BaseSearch(name, k, nSearch) {
    fg = new Graph();
    lastSeed = -1;
  }
  ~GraphSearchTJBase() {
    delete fg;
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
    lastSeed = -1;
    vector<operand> answer;

    assert(ready == true);
    vector<bool> visited(N, false);
    vector<bool> inserted(N, false);
    auto t1 = Clock::now();
    //int init = rand() % N;
    int init = 0;
    for(int i = 0; i < 1024; i++) {
      walk(init, visited, inserted, current_query);
      if(searched > nSearch || init == -1) {
        break;
      }
    }
    auto t2 = Clock::now();
    partial_sort(top_k_candidates.begin(), top_k_candidates.begin() + K, top_k_candidates.end(), compare_descending);
    auto t3 = Clock::now();
    timers[1] += compute_clock(t2, t1);
    timers[2] += compute_clock(t3, t2);
    return top_k_candidates;
  }

  void inline insertCandidate(vector<bool> &visited, int id, float sim) { ///, float* current_query) {
    if(!visited[id]) {
      visited[id] = true;
      top_k_candidates.emplace_back(id, sim);
    }
  }
  void walk(int &seed, vector<bool> &visited, vector<bool> &inserted, float *current_query) {
    int curSeed = seed;
    float curSim  = getSimilarity(seed, current_query);
    bool found = true;
    int searchEdge = 256;
    while(found) {
      visited[curSeed] = true;
      insertCandidate(inserted, curSeed, curSim);
      found = false;
      vector<int> &neighborList = (*fg->graph)[curSeed];
      float neighborMaxSim = numeric_limits<float>::lowest();
      int maxNeighbor = -1;
      int deg = min((int)neighborList.size(), searchEdge);
      for (int i=0; i<deg; i++) {
        int candidate = neighborList[i];
        //cout << "cand : "<< candidate << "\n";
        float neighborSim = getSimilarity(candidate, current_query);
        searched++;
        if(neighborSim > curSim) {
          if(visited[candidate])
            continue;
          curSim = neighborSim;
          curSeed = candidate;
          found = true;
        }
        if(neighborSim > neighborMaxSim) {
          if(visited[candidate])
            continue;
          neighborMaxSim = neighborSim;
          maxNeighbor = candidate;
        }
        //if(found && i>= 16)
        //  break;
      }
      if(!found) {
        seed = maxNeighbor;
      }
    }
  }
  void print_time() {
    cout << "[" << name << "] " << (timers[1] + timers[2]+timers[3] + timers[4]) /queryCount << "(ms)\n";
    cout << "Preprocessing : " << timers[0] << "(ms)\n";
    cout << "Search " << (timers[1] + timers[2]) / queryCount << "(ms)\n";
    cout << "---[1] Walk Front : " << timers[1] / queryCount << "(ms)\n";
    cout << "---[2] Walk : " << timers[2] / queryCount << "(ms)\n";
    cout << "---[3] Expand : " << timers[3] / queryCount << "(ms)\n";
    cout << "---[4] Expand: " << timers[4]/queryCount << "(ms)\n";
    //cout << "---[2-2] Expand-Sort : " << timers[4] / queryCount << "(ms)\n";
   }
  void print_information() {
    cout << "[GraphSearchTJ] \n";
  }
}; 
