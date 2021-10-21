using phmap::flat_hash_map;
using phmap::flat_hash_set;
//#define SPARSIFYTYPE1
#define EXPAND
class GraphSearchFinal2 : public BaseSearch {
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
  int ct3t;
  int load_front_M;
  int check_width;
  int intersection_width;
  vector<int> rank;
  int seed;
  GraphSearchFinal2(string name, int k, int nSearch, int mode_, int load_front_M_, int check_width_, int intersection_width_) : BaseSearch(name, k, nSearch) {
    fg = new Graph();
    mode = mode_;
    ct1 = 0;
    ct2 = 0;
    ct3 = 0;
    ct3t = 0;
    in_walk = 0;
    walk = 0;
    load_front_M = load_front_M_;
    check_width = check_width_;
    intersection_width = intersection_width_;
  }
  ~GraphSearchFinal2() {
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
    seed = 0;
    cout << "Preprocess done\n";
    timers[0] += compute_clock(t2, t1);
    ready = true;
  }
  float inline getSimilarity(int id, float *current_query) {
    float sim = inner_product(data[id], current_query, D);
    return sim;
  }
  vector<operand> search(float* current_query){
    vector<bool> visited(N, false);
    vector<bool> inserted(N, false);
    auto t1 = Clock::now();
    walk_in_graph(seed, current_query, inserted, visited);
    auto t2 = Clock::now();
    partial_sort(top_k_candidates.begin(), top_k_candidates.begin() + K, top_k_candidates.end(), compare_descending);
    auto t3 = Clock::now();
    timers[1] += compute_clock(t2, t1);
    timers[2] += compute_clock(t3, t2);
    if(queryCount % 50 == 0) {
      cout << "Current Query : " << queryCount << "\n";
    }
    return top_k_candidates;
  }
  void inline insertCandidate(vector<bool> &visited, int id, float sim) {
    if(!visited[id]) {
      visited[id] = true;
      top_k_candidates.emplace_back(id, sim);
    }
  }
  int skipcnt = 0;
  void walk_in_graph(int seed, float* current_query, vector<bool> &inserted, vector<bool>& visited) {
    StreamTracker st(8);
    int curNode = seed;
    float curSim = getSimilarity(curNode, current_query);
    while(curNode != -1) {
      walk++;
      insertCandidate(inserted, curNode, curSim);
      assert(visited[curNode] == false);
      visited[curNode] = true;
      vector<int>& neighbor = (*fg->graph)[curNode];
      #ifdef SPARSIFYTYPE1
      neighbor = GraphConstruct::sparsify(neighbor, curNode, true);
      ct3 += neighbor.size();
      ct3t++;
      #endif
      int curNeighborNode = -1;
      float curNeighborSim = numeric_limits<float>::lowest();
      int secondLargestNeighbor = -1;
      float secondLargestNeighborSim = numeric_limits<float>::lowest();
      /* Decide where to Start */
      bool found = false;
      bool cont = true;
      load_front_M = 32;     
      bool skip = neighbor.size() <= 128;  //SKIP
      if (skip) skipcnt++;
      for (int i = 0; i < min(load_front_M, (int)neighbor.size()); i++) {
        ct1++;
        float sim = getSimilarity(neighbor[i], current_query);
        if (sim > curNeighborSim) {
          if(!visited[neighbor[i]]) {
            searched++;
            curNeighborNode = neighbor[i];
            curNeighborSim = sim;
            found = true;
          }
        }
        if(curNeighborSim > sim && sim > secondLargestNeighborSim) {
          if(!visited[neighbor[i]]) {
            secondLargestNeighbor = neighbor[i];
            secondLargestNeighborSim = sim;
          }
        }
        if(found && i >= 15 && !skip)   //SKIP
          break;
      }

      assert(curNeighborNode != -1);
      st.update(curNeighborNode, curNeighborSim);
      if(secondLargestNeighbor != -1)
        st.update(secondLargestNeighbor, secondLargestNeighborSim);
      /* Start walk from NeighborSeed */
      int maxNeighbor = -1;
      float maxSim = numeric_limits<float>::lowest();
      while(cont && !skip) {  //SKIP
        cont = false;
        in_walk++;  
        int maxM = 32;
        maxSim = numeric_limits<float>::lowest();
        vector<operand> mainneighbor;
        for(int i = 0; i < neighbor.size(); i++) {
          if (neighbor[i] != curNeighborNode)
            mainneighbor.push_back(make_pair(neighbor[i], 0));
        }
        #pragma omp parallel for
        for(int i = 0; i < mainneighbor.size(); i++) {
          float sim = inner_product(data[curNeighborNode], data[mainneighbor[i].first], D);
          mainneighbor[i].second = sim;
        }
        int deg = (int)min(maxM, (int)mainneighbor.size());
        partial_sort(mainneighbor.begin(), mainneighbor.begin() + deg, mainneighbor.end(), compare_descending);        
        deg = (int)min(maxM, (int)mainneighbor.size());
        for(int i = 0; i < deg; i++) {
          int neighborNode = mainneighbor[i].first;
          ct2++;
          float neighborSim = getSimilarity(neighborNode, current_query);
          searched++;
          if(neighborSim > curNeighborSim) {
            if(!visited[neighborNode]) {
              curNeighborNode = neighborNode;
              curNeighborSim = neighborSim;
              cont = true;
            }
          }
          #ifdef EXPAND
          if (neighborSim > maxSim) {
            if(!visited[neighborNode]) {
              maxNeighbor = neighborNode;
              maxSim = neighborSim;
            }
          }
          #endif
        }
      }
      if (!skip) {  //SKIP
        st.update(curNeighborNode, curNeighborSim);
        st.update(maxNeighbor, maxSim);
        insertCandidate(inserted, curNeighborNode, curNeighborSim);
      }
      operand nextNode = st.get();
      curNode = nextNode.first;
      curSim = nextNode.second;
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
    cout << "Averge kick out : " << (float)ct3 / ct3t << "\n";
   }
  void print_information() {
    cout << "[GraphSearchYJEmulated] \n";
    cout << skipcnt/(float)queryCount << "\n";
  }
}; 
