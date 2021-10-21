using phmap::flat_hash_map;
using phmap::flat_hash_set;
class GraphSearchYJEmulated : public BaseSearch {
public:
  vector<flat_hash_set<int>> sorted_graph;
  int in_walk;
  int walk;
  Graph* fg;
  int M;
  int mode;
  int ct1;
  int ct2;
  int ct3;
  int load_front_M;
  int check_width;
  int intersection_width;
  vector<int> rank;

  GraphSearchYJEmulated(string name, int k, int nSearch, int mode_, int load_front_M_, int check_width_, int intersection_width_) : BaseSearch(name, k, nSearch) {
    fg = new Graph();
    mode = mode_;
    ct1 = 0;
    ct2 = 0;
    ct3 = 0;
    in_walk = 0;
    walk = 0;
    load_front_M = load_front_M_;
    check_width = check_width_;
    intersection_width = intersection_width_;
  }
  ~GraphSearchYJEmulated() {
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
    sorted_graph.reserve(N);
    for(int i = 0; i < N; i++) {
      vector<int> &neighbor = (*fg->graph)[i];
      flat_hash_set<int> s;
      for(int j = 0; j < check_width; j++) {
        s.insert(neighbor[j]);
      }
      sorted_graph.push_back(s);
    }
    auto t2 = Clock::now();
    timers[0] += compute_clock(t2, t1);
    ready = true;
  }
  float inline getSimilarity(int id, float *current_query) {
    float sim = inner_product(data[id], current_query, D);
    return sim;
  }
  vector<operand> search(float* current_query){
    #ifdef TEST
    vector<operand> answer;
    rank.clear();
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
    vector<bool> visited(N, false);
    vector<bool> inserted(N, false);
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
  void inline insertCandidate(vector<bool> &visited, int id, float sim) { ///, float* current_query) {
    if(!visited[id]) {
      visited[id] = true;
      top_k_candidates.emplace_back(id, sim);
    }
  }
  void walk_in_graph(int seed, float* current_query, vector<bool> &inserted, vector<bool>& visited) {
    int curNode = seed;
    float curSim = getSimilarity(curNode, current_query);
    while(curNode != -1) {
      walk++;
      insertCandidate(inserted, curNode, curSim);
      visited[curNode] = true;
      // cout << "CurNode : "<< curNode << "\n";
      vector<int>& neighbor = (*fg->graph)[curNode];
      int curNeighborNode = -1;
      float curNeighborSim = numeric_limits<float>::lowest();
      /* Decide where to Start */
      bool found = false;
      bool cont = true;
      // flat_hash_map<int, int> curNodeRank;
      // for(int i = 0; i < neighbor.size(); i++) {
      //   curNodeRank.insert(make_pair(neighbor[i], i));
      // }
      load_front_M = 32;
      for (int i = 0; i < load_front_M; i++) {
        searched++;
        ct1++;
        float sim = getSimilarity(neighbor[i], current_query);
        if (sim > curNeighborSim) {
          if(!visited[neighbor[i]]) {
            curNeighborNode = neighbor[i];
            curNeighborSim = sim;
            found = true;
          }
        }
        if(found && i >= 15)
          break;
      }
      /* Start walk from NeighborSeed */
      while(cont) {
        cont = false;
        in_walk++;  
        vector<int> &neighborList = (*fg->graph)[curNeighborNode];
        int reali = 0;
        //cout << "InNode : "<< curNeighborNode << "\n";
        //assert(neighborList.size() >= intersection_width);
        for (size_t i = 0; i < intersection_width; i++) {
          int neighborNode = neighborList[i];
          if(sorted_graph[curNode].count(neighborNode) == 0)
            continue;
          ct2++;
          float neighborSim = getSimilarity(neighborNode, current_query);
          // cout << "neighbor id : " << neighborNode << " , neighborsim : "  << neighborSim << " (curneighborSim : " << curNeighborSim << ") sorted_graph found? : " << sorted_graph[curNode].count(neighborList[i]) << "\n";
          //insertCandidate(inserted, neighborNode, neighborSim);
          searched++;
          if(neighborSim > curNeighborSim) {
            // cout << "visited : " << visited[neighborNode] << "\n";
            if(!visited[neighborNode]) {
              curNeighborNode = neighborNode;
              curNeighborSim = neighborSim;
              cont = true;
            }
          }
          // if(reali == 63)
          //   break;
          reali++;
        }
        //if(cont && curNodeRank[curNeighborNode] >= 512) {
          //cout << "curNeighborNode(" << curNeighborNode << ") is not in minigraph\n";
          //ct3++;
          //break;
        //}
      }
      curNode = curNeighborNode;
      curSim = curNeighborSim;
      if (searched >= nSearch)
        break;
    }
  }

  void print_time() {
    cout << "[" << name << "] " << (timers[1] + timers[2]) /queryCount << "(ms)\n";
    cout << "Preprocessing : " << timers[0] << "(ms)\n";
    cout << "Search " << (timers[1] + timers[2]) / queryCount << "(ms)\n";
    cout << "---[1] Walk Front : " << timers[1] / queryCount << "(ms)\n";
    cout << "---[2] Walk : " << timers[2] / queryCount << "(ms)\n";
    cout << "Average in_graph walk : " << (float)in_walk / walk << "\n";
    cout << "Average out walk : " << (float)walk / queryCount << "\n";
    cout << "Per Node In-walk : "<< (float)ct2 / in_walk << "\n";
    cout << "Average out : " << (float)ct1 / queryCount << "\n";
    cout << "Average in : "  << (float)ct2 / queryCount << "\n";
    cout << "Averge kick out : " << (float)ct3 / queryCount << "\n";
   }
  void print_information() {
    cout << "[GraphSearchYJEmulated] \n";
  }
}; 
