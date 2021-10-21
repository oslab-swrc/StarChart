using phmap::flat_hash_map;
using phmap::flat_hash_set;
//#define LOG
class MemorySearch : public BaseSearch {
public:
  unordered_map<int, unordered_map<int, vector<int>>>* in_graph;
  int in_walk;
  int walk;
  Graph* fg;
  int M;
  int mode;
  int ct1;
  int ct2;
  int itc;
  #ifdef LOG
  vector<operand> answer;
  vector<int> rank;
  #endif
  MemorySearch(string name, int k, int nSearch, int mode_) : BaseSearch(name, k, nSearch) {
    fg = new Graph();
    mode = mode_;
    ct1 = 0;
    ct2 = 0;
    in_walk = 0;
    walk = 0;
    itc = 0;
  }
  ~MemorySearch() {
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
    #ifdef LOG
    /* Statistics */
    rank.clear();
    answer.clear();
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
    #endif
    /* End Statistics */
    flat_hash_set<int> visited;
    flat_hash_set<int> inserted;
    inserted.reserve(4096);
    visited.reserve(1024);
    auto t1 = Clock::now();
    int seed = 0;
    walk_in_graph(seed, current_query, inserted, visited);
    auto t2 = Clock::now();
    partial_sort(top_k_candidates.begin(), top_k_candidates.begin() + K, top_k_candidates.end(), compare_descending);
    auto t3 = Clock::now();
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
  vector<int> intersect(vector<int> &v1, vector<int> &v2) {
    vector<int> res;
    res.reserve(128);
    flat_hash_set<int> orig; 
    for(int i=0; i<v1.size(); i++) {
      orig.insert(v1[i]);
    }
    for(int i=0; i<v2.size(); i++) {
      if(orig.find(v2[i]) != orig.end()) {
        res.push_back(v2[i]);
      }
      if(res.size() == 128)
        break;
    }
    itc++;
    return res;
  }
  void walk_in_graph(int seed, float* current_query, flat_hash_set<int> &inserted, flat_hash_set<int>& visited) {
    int curNode = seed;
    float curSim = getSimilarity(curNode, current_query);
    visited.insert(curNode);
    while(curNode != -1) {
      walk++;
      insertCandidate(inserted, curNode, curSim);
      visited.insert(curNode);
      vector<int>& neighbor = (*fg->graph)[curNode];
      int curNeighborNode = -1;
      float curNeighborSim = numeric_limits<float>::lowest();
      #ifdef LOG
      cout << "Current Out-visit : " << rank[curNode] << "\n";
      for(int i=0; i<neighbor.size(); i++) {
        if(rank[neighbor[i]] < 50) {
          cout << "High Rank Neighbors : " << rank[neighbor[i]] << " / " << i << "\n";
        }
      }
      #endif
      /* Decide where to Start */
      bool found = false;
      bool cont = true;
      for (int i = 0; i < 32; i++) {
        searched++;
        ct1++;
        float sim = getSimilarity(neighbor[i], current_query);
        if (sim > curNeighborSim) {
          if(visited.find(neighbor[i]) != visited.end())
            continue;
          curNeighborNode = neighbor[i];
          curNeighborSim = sim;
          found = true;
        }
        if(found && i >= 15)
          break;
      }
      #ifdef LOG
      cout << "Local Max of Minigraph for " << rank[curNode] << "\n";
      for (auto it = (*in_graph)[curNode].begin(); it != (*in_graph)[curNode].end(); ++it) {
        vector<int> &neighborList = it->second;
        float sim = getSimilarity(it->first, current_query);
        bool isLocalMax = true;
        for (size_t i = 0; i <(int)neighborList.size(); i++) {
          int neighborNode = neighborList[i];
          float neighborSim = getSimilarity(neighborNode, current_query);
          if(neighborSim > sim) { 
            isLocalMax = false;
            break;
          }
        }
        if(isLocalMax) {
          cout << "LocalMax: " << rank[it->first] << "\n";
        }
      }
      #endif
      /* Start walk from NeighborSeed */
      while(cont) {
        cont = false;
        in_walk++;
        #ifndef EMULATE
        vector<int> &neighborList = (*in_graph)[curNode][curNeighborNode];
        #endif
        #ifdef EMULATE
        vector<int> neighborList = intersect((*fg->graph)[curNode],(*fg->graph)[curNeighborNode]); //(*in_graph)[curNode][curNeighborNode];
        #endif
        #ifdef LOG
        cout << "Current In-visit : " << rank[curNeighborNode] << "\n";
        #endif
        for (size_t i = 0; i <(int)neighborList.size(); i++) {
          int neighborNode = neighborList[i];
          ct2++;
          float neighborSim = getSimilarity(neighborNode, current_query);
          searched++;
          if(neighborSim > curNeighborSim) {
            if (visited.find(neighborNode) != visited.end())
              continue;
            curNeighborNode = neighborNode;
            curNeighborSim = neighborSim;
            cont = true;
          }
          if (cont && i >= 15) 
            break;
        }
      }
      curNode = curNeighborNode;
      curSim = curNeighborSim;
      if (searched >= nSearch)
        break;
    }
  }
  void print_time() {
    cout << "[MemorySearch] " << (timers[1] + timers[2]) /queryCount << "(ms)\n";
    cout << "Preprocessing : " << timers[0] << "(ms)\n";
    cout << "Search " << (timers[1] + timers[2]) / queryCount << "(ms)\n";
    cout << "---[1] Walk Front : " << timers[1] / queryCount << "(ms)\n";
    cout << "---[2] Walk : " << timers[2] / queryCount << "(ms)\n";
    cout << "Average in_graph walk : " << (float)in_walk / walk << "\n";
    cout << "Average out walk : " << (float)walk / queryCount << "\n";
    cout << "Average out : " << (float)ct1 / queryCount << "\n";
    cout << "Average in : "  << (float)ct2 / queryCount << "\n";

   }
  void print_information() {
    cout << "[MemorySearch] \n";
  }
}; 
