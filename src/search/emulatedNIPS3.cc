using phmap::flat_hash_map;
using phmap::flat_hash_set;
//#define SPARSIFYTYPE1
//#define SPARSIFYTYPE2
//#define LOG
class GraphSearchFinal : public BaseSearch {
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
  int realwalk;
  int skipwalk;
  int load_front_M;
  vector<int> rank;
  vector<int> seedset;
  int seed;
  int itc1;
  int itc2;
  int skipdegree;
  int skiplimit;
  GraphSearchFinal(string name, int k, int nSearch, int mode_, int load_front_M_, int skipdegree_, int skiplimit_) : BaseSearch(name, k, nSearch) {
    fg = new Graph();
    mode = mode_;
    ct1 = 0;
    ct2 = 0;
    ct3 = 0;
    ct3t = 0;
    in_walk = 0;
    walk = 0;
    skipwalk = 0;
    realwalk = 0;
    itc1 = 0;
    itc2 = 0;
    load_front_M = load_front_M_;
    skipdegree = skipdegree_;
    skiplimit = skiplimit_;
  }
  ~GraphSearchFinal() {
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
    timers[0] += compute_clock(t2, t1);
    // vector<operand> deg;
    // cout << "Finding maxOutDegreeNode starts, fgsize : " << (*fg->graph).size() << "\n";
    // for (int i = 0; i < (*fg->graph).size(); i++){
    //   deg.push_back(make_pair(i, (float)(*fg->graph)[i].size()));
    // }
    // sort(deg.begin(), deg.end(), compare_descending);
    //seed = deg[0].first;
    ready = true;
  }
  float inline getSimilarity(int id, float *current_query) {
    float sim = inner_product(data[id], current_query, D);
    return sim;
  }

  vector<operand> search(float* current_query){
    #ifdef LOG
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
    for(int i = 0; i < N; i++)
      rank[answer[i].first] = i;
    for(int i = 0; i < 50; i++) {
      vector<int>& neighbor = (*fg->graph)[answer[i].first];
      cout << "Myneighbor : " << i << " / " << neighbor.size() << "\n";
      for (int j = 0; j < min(512, (int)neighbor.size()); j++) {
        if(rank[neighbor[j]] < 100)
          cout << " i : " << i << " friend : " << j <<"  = " <<  rank[neighbor[j]] << "\n";
      }
    }
    #endif
    StreamTracker st(8);
    
    vector<bool> visited(N, false);
    vector<bool> inserted(N, false);
    auto t1 = Clock::now();
    walk_in_graph(seed, current_query, inserted, visited, st);
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
  void walk_in_graph(int seed, float* current_query, vector<bool> &inserted, vector<bool>& visited, StreamTracker &st) {
    int curNode = seed;
    float curSim = getSimilarity(curNode, current_query);
    while(curNode != -1) {
      walk++;
      insertCandidate(inserted, curNode, curSim);
      visited[curNode] = true;
      vector<int> neighbor = (*fg->graph)[curNode];
      #ifdef LOG
      int mid = -1;
      int maxrank = 9999999;
      for (int i=0; i<neighbor.size(); i++) {
        if(rank[neighbor[i]] <= 50 && !visited[neighbor[i]])
         cout << "unvisited good neighbor : " << rank[neighbor[i]] << "\n";
        if(rank[neighbor[i]] < maxrank) {
         mid = neighbor[i];
         maxrank = rank[neighbor[i]];
        }
      }
      cout << "MaxRank : " << maxrank << " / " << mid << " \n";
      cout << "Curnode : " << curNode << " / rank : " << rank[curNode] << " / " << neighbor.size() << "\n";
      #endif
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
      int deg = min(load_front_M, (int)neighbor.size());
      int maxM = 32;
      int ipass = 15;
      bool skip = false;
      flat_hash_set<int> ivisit;

      skip = false;
      if (neighbor.size() <= skiplimit) {
        deg =  min((int)neighbor.size(), skipdegree);
        skip = true;
      }
      for (int i = 0; i < deg ; i++) {
        itc1++;
        searched++;
        ct1++;
        float sim = getSimilarity(neighbor[i], current_query);
        ivisit.insert(neighbor[i]);
        bool cond1 = sim > curNeighborSim;
        bool cond2 = sim > secondLargestNeighborSim;
        if((cond1 || cond2) && !visited[neighbor[i]]) {
          bool cond1 = sim > curNeighborSim;
          bool cond2 = sim > secondLargestNeighborSim;
          if (cond1 || cond2) {
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
      }
      #ifdef LOG
      if (curNeighborNode != -1)
       cout << "Selected Curnode : " << curNeighborNode << " / rank : " << rank[curNeighborNode] << "\n";
      #endif
      if (curNeighborNode == -1) {
        if (st.size > 0) {
          operand get = st.get();
          curNode = get.first;
          curSim = get.second;
          continue;
        }
        else if (st.size == 0)
          break;
      }
      st.update(curNeighborNode, curNeighborSim);
      if(secondLargestNeighbor != -1)
        st.update(secondLargestNeighbor, secondLargestNeighborSim);

      if(skip) {
        operand nextNode = st.get();
        curNode = nextNode.first;
        curSim = nextNode.second;
        skipwalk++;
        if (searched >= nSearch)
          break;
        continue;
      }
      
      /* Start walk from NeighborSeed */
      realwalk++;
      int budget = neighbor.size()/4; //neighbor.size()/2;
      // int budget = 0;
      while(cont) {
        ivisit.insert(curNeighborNode);
        #ifdef LOG
        cout << "Curneighbor : " << curNeighborNode << " / rank : " << rank[curNeighborNode] << "\n";
        #endif
        cont = false;
        in_walk++;  
        secondLargestNeighbor = -1;
        secondLargestNeighborSim = numeric_limits<float>::lowest();

        ///////////
        vector<operand> mainneighbor;
        for(int i = 0; i < neighbor.size(); i++) {
          if (neighbor[i] != curNeighborNode)
            mainneighbor.push_back(make_pair(neighbor[i], 0));
        }
        #pragma omp parallel for num_threads(16)
        for(int i = 0; i < mainneighbor.size(); i++) {
          float sim = inner_product(data[curNeighborNode], data[mainneighbor[i].first], D);
          mainneighbor[i].second = sim;
        }
        int deg = (int)min(maxM, (int)mainneighbor.size());
        partial_sort(mainneighbor.begin(), mainneighbor.begin() + deg, mainneighbor.end(), compare_descending);        
        #ifdef SPARSIFYTYPE2
        mainneighbor = GraphConstruct::sparsify(mainneighbor, curNeighborNode, deg);
        #endif
        ////////////

        deg = (int)min(maxM, (int)mainneighbor.size());
        for(int i = 0; i < deg; i++) {
          itc2++;
          int neighborNode = mainneighbor[i].first;
          if(!visited[neighborNode] && (ivisit.find(neighborNode) == ivisit.end())) {
            ct2++;
            searched++;
            budget--;

            float neighborSim = getSimilarity(neighborNode, current_query);
            ivisit.insert(neighborNode);
            bool cond1 = neighborSim > curNeighborSim;
            bool cond2 = neighborSim > secondLargestNeighborSim;
            if (cond1 || cond2) {
              if(cond1) {
                curNeighborSim = neighborSim;
                curNeighborNode = neighborNode;
                cont = true;
              }
              else if(cond2) {
                secondLargestNeighbor = neighborNode;
                secondLargestNeighborSim = neighborSim;
              }
            }
          }
          //if(cont && i >= 15)
          //  break;
        }
        if(!cont && budget > 0 && (secondLargestNeighbor != -1)) {
          cont = true;
          st.update(curNeighborNode, curNeighborSim);
          curNeighborNode = secondLargestNeighbor;
          curNeighborSim = secondLargestNeighborSim;
          secondLargestNeighbor = -1;
          secondLargestNeighborSim = numeric_limits<float>::lowest();
        }
      }
      #ifdef LOG
      cout << "Final Local Max : " << curNeighborNode << " / rank : " << rank[curNeighborNode] << "\n";
      if (secondLargestNeighbor != -1)
       cout << "Second larget neighbor : " << secondLargestNeighbor << " / rank : " << rank[secondLargestNeighbor] << "\n";
      #endif

      st.update(curNeighborNode, curNeighborSim);
      if(secondLargestNeighbor != -1)
        st.update(secondLargestNeighbor, secondLargestNeighborSim);
      insertCandidate(inserted, curNeighborNode, curNeighborSim);
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
    cout << "Average Large walk : " << (float)walk / queryCount << "\n";
    cout << "Average Real walk : " << (float)realwalk / queryCount << "\n";
    cout << "Average steps in IG walk : " << (float)in_walk / realwalk << "\n";
    cout << "Skip Walk : "<< (float)skipwalk / walk << "\n";
    cout << "Per Node In-walk : "<< (float)ct2 / in_walk << "\n";
    cout << "Average Walk Spent on Out : " << (float)ct1 / queryCount << "\n";
    cout << "Average Walk Spent on In : "  << (float)ct2 / queryCount << "\n";
    cout << "Averge Degree : " << (float)ct3 / ct3t << "\n";
    cout << "Out walk Iter Portion: " << (float) ct1 / itc1 << "\n";
    cout << "In walk Iter Portion: " << (float) ct2  / itc2 << "\n";
    cout << "Skipdegree : " << skipdegree << "\n";
    cout << "SkipLimit : " << skiplimit << "\n";

  }
  void print_information() {
    cout << "[GraphSearchYJEmulated] \n";
  }
}; 
