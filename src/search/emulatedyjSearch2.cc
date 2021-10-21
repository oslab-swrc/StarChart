using phmap::flat_hash_map;
using phmap::flat_hash_set;
//#define SPARSIFYTYPE1
// #define SPARSIFYTYPE2
#define EXPAND
class GraphSearchYJEmulated2 : public BaseSearch {
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
  GraphSearchYJEmulated2(string name, int k, int nSearch, int mode_, int load_front_M_, int check_width_, int intersection_width_) : BaseSearch(name, k, nSearch) {
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
  ~GraphSearchYJEmulated2() {
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
    // vector<flat_hash_map<int, vector<int>>> additional_edges;
    // for (int i = 0; i < N; i++)
    //   additional_edges.push_back(flat_hash_map<int, vector<int>>());
    flat_hash_map<int, vector<int>> sparsified_graph;
    long long int tsize = 0;
    long long int t = 0;
    for (int i =0; i < N; i++) {
      vector<int>& neighbor = (*fg->graph)[i];
      #ifdef SPARSIFYTYPE1
      neighbor = GraphConstruct::sparsify(neighbor, i, true);
      tsize += neighbor.size();
      t++;
      #endif
      sparsified_graph.insert(make_pair(i, neighbor));
    }
    cout << "Average neighbor size : " << (float)tsize / t << "\n";
    set<int> reachable;
    set<int> unreachable;
    cout << "Reachable size : " << reachable.size() << " unreachable size : " << unreachable.size() << "\n";
    bfs(seed, reachable, unreachable, sparsified_graph);
    cout << "Preprocess done\n";
    timers[0] += compute_clock(t2, t1);
    ready = true;
  }
  float inline getSimilarity(int id, float *current_query) {
    float sim = inner_product(data[id], current_query, D);
    return sim;
  }
  // int get_start_ingraph_seed(vector<int>& neighbor) {
  //   float* centroid = get_centroid(neighbor);
  //   int start_seed = find_seed(centroid, neighbor);

  //   // flat_hash_map<int, vector<int>> ingraph = get_ingraph(seed, neighbor);
  //   // float[] centroid = get_centroid_with_ingraph(ingraph);
  //   // int start_seed = find_seed_with_ingraph(centroid, ingraph);
  //   // cout << "Starting Seed " << i << " : " << start_seeds[i] << "\n";
  //   // make_path(start_seed, ingraph, additional_edges);
  //   delete[] centroid;
  //   return start_seed;
  // }
  // flat_hash_map<int, vector<int>> get_ingraph(int node, vector<int>& neighbor) {
  //   int maxM = 32;
  //   flat_hash_map<int, vector<int>> ingraph;
  //   for (int j = 0; j < neighbor.size(); j++)
  //     ingraph.insert(make_pair(neighbor[j], vector<int>()));
    
  //   #pragma omp parallel for
  //   for (int j = 0; j < neighbor.size(); j++) {
  //     int curNeighborNode = neighbor[j];
  //     vector<operand> mainneighbor;
  //     for(int i = 0; i < neighbor.size(); i++) {
  //       if (neighbor[i] != curNeighborNode)
  //         mainneighbor.push_back(make_pair(neighbor[i], 0));
  //     }
  //     for(int i = 0; i < mainneighbor.size(); i++) {
  //       float sim = inner_product(data[curNeighborNode], data[mainneighbor[i].first], D);
  //       mainneighbor[i].second = sim;
  //     }

  //     int deg = (int)min(maxM, (int)mainneighbor.size());
  //     partial_sort(mainneighbor.begin(), mainneighbor.begin() + deg, mainneighbor.end(), compare_descending);
  //     // cout << "deg : " << deg << "\n";
  //     mainneighbor.resize(deg);
  //     for (int i = 0; i < mainneighbor.size(); i++)
  //       ingraph[curNeighborNode].push_back(mainneighbor[i].first);
  //   }
  //   return ingraph;
  // }
  float* get_centroid(vector<int>& neighbor) {
    float* centroid = new float[D];
    #pragma omp parallel for
    for (int j = 0; j < D; j++)
      centroid[j] = 0.0f;
    for (int i = 0; i < neighbor.size(); i++) {
      #pragma omp parallel for
      for (int j = 0; j < D; j++) 
        centroid[j] += data[neighbor[i]][j];
    }
    #pragma omp parallel for
    for (int j = 0; j < D; j++)
      centroid[j] /= neighbor.size();
    return centroid;
  }
  // float* get_centroid_with_ingraph(flat_hash_map<int, vector<int>>& ingraph) {
  //   float* centroid = new float[D];
  //   #pragma omp parallel for
  //   for (int j = 0; j < D; j++)
  //     centroid[j] = 0.0f;
  //   for (auto it : ingraph) {
  //     #pragma omp parallel for
  //     for (int j = 0; j < D; j++)
  //       centroid[j] += data[it.first][j];
  //   }
  //   #pragma omp parallel for
  //   for (int j = 0; j < D; j++)
  //     centroid[j] /= ingraph.size();
  //   return centroid;
  // }
  int find_seed(float* centroid, vector<int>& neighbors) {
    float maxSim = numeric_limits<float>::lowest();
    int maxID = -1;
    for (int i = 0; i < neighbors.size(); i++) {
      float sim = getSimilarity(neighbors[i], centroid);
      if (sim > maxSim) {
        // cout << "neighbor : " << neighbors[i] << " sim : " << sim << "\n";
        maxSim = sim;
        maxID = neighbors[i];
      }
    }
    if (maxID == -1) {
      cout << "neighbor size : " << neighbors.size() << "\n";
      for (int i = 0; i < D; i++)
        cout << "centroid : " << centroid[i] << "\n";
    }
    assert(maxID != -1);
    return maxID;
  }
  // int find_seed_with_ingraph(float* centroid, flat_hash_map<int, vector<int>>& ingraph) {
  //   float maxSim = numeric_limits<float>::lowest();
  //   int maxID = -1;

  //   for (auto it : ingraph) {
  //     float sim = getSimilarity(it.first, centroid);
  //     if (sim > maxSim) {
  //       maxSim = sim;
  //       maxID = it.first;
  //     }
  //   }
  //   return maxID;
  // }
  void make_path(int seed, flat_hash_map<int, vector<int>>& ingraph, vector<flat_hash_map<int, vector<int>>>& additional_edges) {
    set<int> reachable;
    set<int> unreachable;
    // cout << "BFS start\n";
    bfs(seed, reachable, unreachable, ingraph);
    // cout << "BFS end \n";
    cout << "reachable : " << reachable.size() << " unnrechable : " << unreachable.size() << "\n";
    for (auto it : unreachable) {
      float maxSim = numeric_limits<float>::lowest();
      int maxID = -1;
      for (auto itit : reachable) {
        float sim = inner_product(data[it], data[itit], D);
        if (sim > maxSim) {
          maxSim = sim;
          maxID = itit;
        }
      }
      additional_edges[seed][maxID].push_back(it);
      cout << "additional edges from " << maxID << " to " << it << "\n";
    }
  }

  void bfs(int seed, set<int>& reachable, set<int>& unreachable, flat_hash_map<int, vector<int>>& ingraph) {
    queue<int> q;
    vector<bool> visited(N, false);
    q.push(seed);
    reachable.insert(seed);
    visited[seed] = true;

    while(!q.empty()){
      int next_node = q.front();
      q.pop();
      // cout << "next node : " << next_node << "\n";
      assert(ingraph[next_node].size() == 32);
      for(int i = 0; i < ingraph[next_node].size(); i++){
        // cout << "nodenodenonde : " << ingraph[next_node][i] << " i : " << i << "\n";
        if(visited[ingraph[next_node][i]] == false) {
          q.push(ingraph[next_node][i]);
          // cout << "pushed : " << ingraph[next_node][i] << "\n";
          reachable.insert(ingraph[next_node][i]);
          visited[ingraph[next_node][i]] = true;
        }
      }
    }

    for (auto it : ingraph) {
      if (reachable.find(it.first) == reachable.end())
        unreachable.insert(it.first);
    }
    assert(reachable.size() + unreachable.size() == ingraph.size());
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
    for(int i = 0; i < N; i++)
      rank[answer[i].first] = i;
    #endif
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
    // cout << "id : " << id << "\n";
    if(!visited[id]) {
      visited[id] = true;
      top_k_candidates.emplace_back(id, sim);
    }
  }
  void walk_in_graph(int seed, float* current_query, vector<bool> &inserted, vector<bool>& visited) {
    #ifdef EXPAND
    StreamTracker st(8);
    #endif
    int curNode = seed;
    float curSim = getSimilarity(curNode, current_query);
    while(curNode != -1) {
      walk++;
      insertCandidate(inserted, curNode, curSim);
      visited[curNode] = true;
      st.deleteItem(curNode);
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
      for (int i = 0; i < min(load_front_M, (int)neighbor.size()); i++) {
        ct1++;
        float sim = getSimilarity(neighbor[i], current_query);
        if (sim > curNeighborSim) {
          if(!visited[neighbor[i]]) {
            searched++;
            curNeighborNode = neighbor[i];
            curNeighborSim = sim;
            //st.update(curNeighborNode, curNeighborSim);
            found = true;
          }
        }
        if(curNeighborSim > sim && sim > secondLargestNeighborSim) {
          if(!visited[neighbor[i]]) {
            secondLargestNeighbor = neighbor[i];
            secondLargestNeighborSim = sim;
          }
        }
        if(found && i >= 15)
          break;
      }
      if (curNeighborNode == -1) 
        cont = false;
      else
        st.update(curNeighborNode, curNeighborSim);
      if(secondLargestNeighbor != -1)
        st.update(secondLargestNeighbor, secondLargestNeighborSim);

      //curNeighborNode = get_start_ingraph_seed(curNode, neighbor);
      //curNeighborSim = getSimilarity(curNeighborNode, current_query);

      /* Start walk from NeighborSeed */
      #ifdef EXPAND
      int maxNeighbor = -1;
      float maxSim = numeric_limits<float>::lowest();
      #endif
      while(cont) {
        cont = false;
        in_walk++;  
        int maxM = 32;
        #ifdef EXPAND
        maxSim = numeric_limits<float>::lowest();
        #endif
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
        int reali = 0;
        #ifdef EXPAND
        int tmp = curNeighborNode;
        #endif
        #ifdef SPARSIFYTYPE2
        mainneighbor = GraphConstruct::sparsify(mainneighbor, curNeighborNode, maxM);
        #endif
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
      #ifdef EXPAND
      if (curNeighborNode !=-1) {
        st.update(curNeighborNode, curNeighborSim);
        insertCandidate(inserted, curNeighborNode, curNeighborSim);
      }
      if (maxNeighbor != -1)
        st.update(maxNeighbor, maxSim);
      operand nextNode = st.get();
      curNode = nextNode.first;
      assert(visited[curNode] == false);

      curSim = nextNode.second;
      #endif
      #ifndef EXPAND
      curNode = curNeighborNode;
      curSim = curNeighborSim;
      #endif
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
  }
}; 
