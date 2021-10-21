using phmap::flat_hash_map;
using phmap::flat_hash_set;
//#define SPARSIFYTYPE1
//#define ADVANCEDSEED
class GraphSearchFinal : public BaseSearch {
public:
  vector<flat_hash_set<int>> sorted_graph;
  int in_walk;
  int walk;
  Graph* fg;
  int M;
  int ct1;
  int ct2;
  int ct3;
  int ct3t;
  int realwalk;
  int skipwalk;
  int load_front_M;
  vector<int> rank;
  int seed;
  GraphSearchFinal(string name, int k, int nSearch, int load_front_M_) : BaseSearch(name, k, nSearch) {
    fg = new Graph();
    ct1 = 0;
    ct2 = 0;
    ct3 = 0;
    ct3t = 0;
    in_walk = 0;
    walk = 0;
    skipwalk = 0;
    realwalk = 0;
    load_front_M = load_front_M_;
  }
  ~GraphSearchFinal() {
  	delete fg;
  }
  void param_register(vector<int> params) {
  }
  void preprocess(){
    auto t1 = Clock::now();
    auto t2 = Clock::now();
    seed = 0;
    timers[0] += compute_clock(t2, t1);
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
    StreamTracker st(16);
    seed = 0;
    vector<bool> visited(N, false);
    vector<bool> inserted(N, false);
    auto t1 = Clock::now();
    walk_in_graph(seed, current_query, inserted, visited, st);
    auto t2 = Clock::now();
    int magic = 4;
    //partial_sort(buffer.begin(), buffer.begin() + K * magic, buffer.end(), compare_descending);
    sort(buffer.begin(), buffer.end(), compare_descending_strict);
    int last = -1;
    int iter = ((int)buffer.size());
    for(int i=0; i< iter; i++) {
      if(buffer[i].first != last) {
        top_k_candidates.push_back(buffer[i]);
        last = buffer[i].first;
      }
      if(top_k_candidates.size() >= K)
        break;
    }
    //unique(top_k_candidates.begin(), top_k_candidates.begin() + K * magic);
    auto t3 = Clock::now();  
    flat_hash_set<int> valid;
    for(int i=0; i<K; i++) {
      if(valid.count(top_k_candidates[i].first) != 0) {
        cout << "[CRITICAL] : Increase Threshold in PARTIAL SORT \n";
        cout << "Current nSearch : "<< nSearch << "\n";
        for(int j=0; j<K; j++) {
          cout << "j: " << top_k_candidates[j].first << " / " << top_k_candidates[j].second <<"\n";
        }
        assert(false);
      }
      valid.insert(top_k_candidates[i].first);
    }

    timers[1] += compute_clock(t2, t1);
    timers[2] += compute_clock(t3, t2);
    if(queryCount % 50 == 0) {
      cout << "Current Query : " << queryCount << "\n";
    }
    return top_k_candidates;
  }
  void inline insertCandidate(vector<bool> &visited, int id, float sim) {
    //if(!visited[id]) {
    //  visited[id] = true;
    buffer.emplace_back(id, sim);
    //}
  }
  void walk_in_graph(int seed, float* current_query, vector<bool> &inserted, vector<bool>& visited, StreamTracker &st) {
    int curNode = seed;
    float curSim = getSimilarity(curNode, current_query);
    st.update(curNode, curSim);
    while(curNode != -1 && searched < nSearch) {
      walk++;
      operand c = st.get();
      curNode = c.first;
      curSim = c.second;
      #ifdef LOG
      cout << "Curnode : " << curNode << " / rank : " << rank[curNode] << "\n";
      #endif
      insertCandidate(inserted, curNode, curSim);
      assert(visited[curNode] == false);
      visited[curNode] = true;
      vector<int> neighbor = (*fg->graph)[curNode];
      #ifdef LOG
      for (int i=0; i<neighbor.size(); i++) {
       if(rank[neighbor[i]] <= 50 && !visited[neighbor[i]])
         cout << "unvisited good neighbor : " << rank[neighbor[i]] << "\n";
      }
      #endif
      #ifdef SPARSIFYTYPE1
      neighbor = GraphConstruct::sparsify(neighbor, curNode, true);
      ct3 += neighbor.size();
      ct3t++;
      #endif
      #ifdef ADVANCEDSEED
      vector<int> myNeighbor = GraphConstruct::sparsify_seed(neighbor, curNode, neighbor.size(), 256);
      #endif
      int curNeighborNode = -1;
      float curNeighborSim = numeric_limits<float>::lowest();
      int secondLargestNeighbor = -1;
      float secondLargestNeighborSim = numeric_limits<float>::lowest();
      /* Decide where to Start */
      bool found = false;
      bool cont = true;
      bool skip = false;
      int maxM = 64;
      flat_hash_set<int> frontcheck;
      int deg = min(load_front_M, (int)neighbor.size());
      if (neighbor.size() <= 256) {
        //deg =  neighbor.size(); 
        deg = min((int)neighbor.size(), 192);
        skip = true;
      }
      #ifdef ADVANCEDSEED
      else {
        deg = min(deg, (int)myNeighbor.size());
      }
      #endif
      int noden = -1;
      for (int i = 0; i < deg; i++) {
        #ifndef ADVANCEDSEED
        noden = neighbor[i];
        #else
        if(!skip)
          noden = myNeighbor[i];
        else
          noden = neighbor[i];
        #endif
        frontcheck.insert(noden);
        searched++;
        float sim = getSimilarity(noden, current_query);
        bool cond1 = sim > curNeighborSim;
        bool cond2 = sim > secondLargestNeighborSim;
        if((cond1 || cond2) && !visited[noden]) {
          if(cond1) {
            secondLargestNeighborSim = curNeighborSim; 
            secondLargestNeighbor = curNeighborNode;
            curNeighborSim = sim;
            curNeighborNode = noden;
          }
          else if(cond2) {
            secondLargestNeighbor = noden;
            secondLargestNeighborSim = sim;
          }
        }
      }

      #ifdef LOG
      if (curNeighborNode != -1)
       cout << "Selected Curnode : " << curNeighborNode << " / rank : " << rank[curNeighborNode] << "\n";
      #endif

      if (curNeighborNode == -1)
        continue;

      if(curNeighborNode != -1)
        st.update(curNeighborNode, curNeighborSim);
      if(secondLargestNeighbor != -1)
        st.update(secondLargestNeighbor, secondLargestNeighborSim);
      
      if(skip) {
        skipwalk++;
        continue;
      }
      
      secondLargestNeighbor = -1;
      secondLargestNeighborSim = numeric_limits<float>::lowest();
      realwalk++;
      while(cont) {
        #ifdef LOG
        cout << "Curneighbor : " << curNeighborNode << " / rank : " << rank[curNeighborNode] << "\n";
        #endif
        cont = false;
        in_walk++;  
        vector<operand> mainneighbor;
        for(int i = 0; i < neighbor.size(); i++) {
          if (neighbor[i] != curNeighborNode)
            mainneighbor.push_back(make_pair(neighbor[i], 0));
        }
        #pragma omp parallel for num_threads(OMD)
        for(int i = 0; i < mainneighbor.size(); i++) {
          float sim = inner_product(data[curNeighborNode], data[mainneighbor[i].first], D);
          mainneighbor[i].second = sim;
        }
        int deg = (int)min(maxM, (int)mainneighbor.size());
        //int deg = mainneighbor.size();
        partial_sort(mainneighbor.begin(), mainneighbor.begin() + deg, mainneighbor.end(), compare_descending);        

        for(int i = 0; i < deg; i++) {
          int neighborNode = mainneighbor[i].first;
          if(frontcheck.count(neighborNode) != 0)
            continue;
          ct2++;
          searched++;
          float neighborSim = getSimilarity(neighborNode, current_query);
          bool cond1 = neighborSim > curNeighborSim;
          bool cond2 = neighborSim > secondLargestNeighborSim;
          if (cond1 || cond2) {
            if(cond1 && !visited[neighborNode]) {
              secondLargestNeighbor = curNeighborNode;
              secondLargestNeighborSim = curNeighborSim;
              curNeighborSim = neighborSim;
              curNeighborNode = neighborNode;
              cont = true;
            }
            else if(cond2 && !visited[neighborNode]) {
              secondLargestNeighbor = neighborNode;
              secondLargestNeighborSim = neighborSim;
            }
          }
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
   }
  void print_information() {
    cout << "[GraphSearchYJEmulated] \n";
  }
}; 
