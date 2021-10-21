using phmap::flat_hash_map;
using phmap::flat_hash_set;

#ifdef FIGUREOUR
extern int nSearch_ = 1000;
#endif

class GraphSearchYJBase : public BaseSearch {
public:
  vector<flat_hash_map<int, vector<int>>>* in_graph;
  int in_walk;
  int walk;
  Graph* fg;
  int M;
  int mode;
  int ct1;
  int ct2;
  int ct3;
  int ct3t;
  int realwalk;
  int skipwalk;

  int load_front_M;
  int seed;
  int skipdegree;
  int skiplimit;
  int itc1;
  int itc2;

  GraphSearchYJBase(string name, int k, int nSearch, int mode_, int load_front_M_, int skipdegree_, int skiplimit_) : BaseSearch(name, k, nSearch) {
    fg = new Graph();
    mode = mode_;
    ct1 = 0;
    ct2 = 0;
    ct3 = 0;
    ct3t = 0;
    in_walk = 0;
    walk = 0;
    realwalk = 0;
    skipwalk = 0;
    itc1 = 0;
    itc2 = 0;
    load_front_M = load_front_M_;
    seed = -1;
    skipdegree = skipdegree_;
    skiplimit = skiplimit_;
  }
  ~GraphSearchYJBase() {
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
    #ifdef FIGUREOUR
    nSearch = nSearch_;
    for (int i = 0; i < 10; i++)
      timers[i] = 0;
    searched = 0;
    totalSearched = 0;
    queryCount = 0;
    #endif
    // int maxNeighborNum = numeric_limits<int>::lowest();
    // int maxID = -1;
    // cout << "Finding maxOutDegreeNode starts, fgsize : " << (*fg->graph).size() << "\n";
    // for (int i = 0; i < (*fg->graph).size(); i++){
    //   if ((int)(*fg->graph)[i].size() > maxNeighborNum) {
    //     maxNeighborNum = (*fg->graph)[i].size();
    //     maxID = i;
    //    }
    // }
    // seed = maxID;
    // cout << "MaxOutDegreeNode : " << maxID << "\n";
    seed = 0;
    auto t2 = Clock::now();
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
    StreamTracker st(8);
    float curSim = getSimilarity(curNode, current_query);
    while(curNode != -1) {
      auto t0 = Clock::now();
      walk++;
      insertCandidate(inserted, curNode, curSim);
      visited[curNode] = true;
      vector<int>& neighbor = (*fg->graph)[curNode];
      int neighborSize = neighbor.size();
      int curNeighborNode = -1;
      float curNeighborSim = numeric_limits<float>::lowest();
      int secondLargestNeighbor = -1;
      float secondLargestNeighborSim = numeric_limits<float>::lowest();
      bool found = false;
      bool cont = true;
      int deg = min(load_front_M, neighborSize);
      int maxM = 32;
      auto t1 = Clock::now();
      timers[3] += compute_clock(t1, t0);
      for (int i = 0; i < deg; i++) {
        searched++;
        float sim = getSimilarity(neighbor[i], current_query);
        //ivisit.insert(neighbor[i]);
        bool cond1 = sim > curNeighborSim;
        bool cond2 = sim > secondLargestNeighborSim;
        if((cond1 || cond2) && !visited[neighbor[i]]) {
            if(cond1) {
              secondLargestNeighborSim = curNeighborSim; 
              secondLargestNeighbor = curNeighborNode;
              curNeighborSim = sim;
              curNeighborNode = neighbor[i];
              found = true;
            }
            else if(cond2) {
              secondLargestNeighbor = neighbor[i];
              secondLargestNeighborSim = sim;
            }
        }
      }
      auto t2 = Clock::now();
      timers[4] += compute_clock(t2,t1);
      if (searched >= nSearch)
        break;

      if (curNeighborNode == -1) {
        if (st.size > 0) {
          operand get = st.get();
          curNode = get.first;
          curSim = get.second;
          auto t11 = Clock::now();
          timers[5] += compute_clock(t11, t2);
          continue;
        }
        else
          break;
      }
      auto t12 = Clock::now();
      st.update(curNeighborNode, curNeighborSim);
      if(secondLargestNeighbor != -1)
        st.update(secondLargestNeighbor, secondLargestNeighborSim);
      auto t3 = Clock::now();
      timers[6] += compute_clock(t3, t12);

      auto t4 = Clock::now();

      /* Start walk from NeighborSeed */
      realwalk++;
      // int budget = neighborSize/4;
      // bool end = false;
      // int counter = 0;
      while(cont) {
        //ivisit.insert(curNeighborNode);
        auto t5 = Clock::now();

        cont = false;
        in_walk++;  
        secondLargestNeighbor = -1;
        secondLargestNeighborSim = numeric_limits<float>::lowest();
        vector<int> &neighborList = (*in_graph)[curNode][curNeighborNode];
        // deg = (int)min(maxM, (int)neighborList.size());
        deg = neighborList.size();
        auto t6 = Clock::now();
        timers[8] += compute_clock(t6, t5);
        for(int i = 0; i < deg; i++) {
          itc2++;
          // counter++;
          int neighborNode = neighborList[i];
          //if((ivisit.find(neighborNode) == ivisit.end())) {
            ct2++;
            searched++;
            // budget--;

            float neighborSim = getSimilarity(neighborNode, current_query);
            //ivisit.insert(neighborNode);
            bool cond1 = neighborSim > curNeighborSim;
            bool cond2 = neighborSim > secondLargestNeighborSim;
            if (cond1 || cond2) {
              if(cond1 && !visited[neighborNode]) {
                curNeighborSim = neighborSim;
                curNeighborNode = neighborNode;
                cont = true;
              }
              else if(cond2 && !visited[neighborNode]) {
                secondLargestNeighbor = neighborNode;
                secondLargestNeighborSim = neighborSim;
              }
            }
          //}
          // if (counter >= budget * 2 && !cont) {
          //   end = true;
          //   break;
          // }
        }
        auto t7 = Clock::now();
        timers[9] += compute_clock(t7, t6);

        // if(!cont /*&& budget > 0 && !end*/ && (secondLargestNeighbor != -1)) {
        //   cont = true;
        //   st.update(curNeighborNode, curNeighborSim);
        //   curNeighborNode = secondLargestNeighbor;
        //   curNeighborSim = secondLargestNeighborSim;
        //   // secondLargestNeighbor = -1;
        //   // secondLargestNeighborSim = numeric_limits<float>::lowest();
        // }
        // auto t8 = Clock::now();
      }
      auto t9 = Clock::now();
      timers[10] += compute_clock(t9, t4);

      st.update(curNeighborNode, curNeighborSim);
      if(secondLargestNeighbor != -1)
        st.update(secondLargestNeighbor, secondLargestNeighborSim);
      insertCandidate(inserted, curNeighborNode, curNeighborSim);
      operand nextNode = st.get();
      curNode = nextNode.first;
      curSim = nextNode.second;
      auto t10 = Clock::now();
      timers[11] += compute_clock(t10, t9);

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
    cout << "Average Large walk : " << (float)walk / queryCount << "\n";
    cout << "Average Real walk : " << (float)realwalk / queryCount << "\n";
    cout << "Average steps in IG walk : " << (float)in_walk / realwalk << "\n";
    cout << "Skip Walk : "<< (float)skipwalk / walk << "\n";
    cout << "Per Node In-walk : "<< (float)ct2 / in_walk << "\n";
    cout << "Average Walk Spent on Out : " << (float)ct1 / queryCount << "\n";
    cout << "Average Walk Spent on In : "  << (float)ct2 / queryCount << "\n";
    cout << "\n\n";
    cout << " ~ before load front M search : " << (timers[3] / queryCount) << "(ms)\n";
    cout << "load front M search : " << (timers[4] / queryCount) << "(ms)\n";
    cout << "st.get : " << (timers[5] / queryCount) << "(ms)\n";
    cout << "st.update : " << (timers[6] / queryCount) << "(ms)\n";
    cout << "skipwalk (st.get) : " << (timers[7] / queryCount) << "(ms)\n";
    cout << "(*in_graph)[curNode][curNeighborNode] : " << (timers[8] / queryCount) << "(ms)\n";
    cout << "ingraph iteration : " << (timers[9] / queryCount) << "(ms)\n";
    cout << "Whole while statement: " << (timers[10] / queryCount) << "(ms)\n";
    cout << "st update ant get: " << (timers[11] / queryCount) << "(ms)\n";
   }
  void print_information() {
    cout << "[GraphSearchYJBase] \n";
  }
}; 
