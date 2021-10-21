// SPDX-FileCopyrightText: Copyright Copyright 2021 Seoul National University
//
// SPDX-License-Identifier: MIT License

#ifdef FIGUREOUR
extern int nSearch_ = 1000;
extern int K_ = -1;
#endif

class GraphSearchYJBase : public BaseSearch {
public:
  vector<flat_hash_map<int, vector<int>>>* in_graph;
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
  int seed;
  int skipdegree;
  int skiplimit;
  int itc1;
  int itc2;

  GraphSearchYJBase(string name, int k, int nSearch, int load_front_M_, int skipdegree_, int skiplimit_) : BaseSearch(name, k, nSearch) {
    fg = new Graph();
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
    seed = 0;
    skipdegree = skipdegree_;
    skiplimit = skiplimit_;
  }
  ~GraphSearchYJBase() {
    delete fg;
  }
  void param_register(vector<int> params) {
  }
  void preprocess(){
    ready = true;
  }
  float inline getSimilarity(int id, float *current_query) {
    float sim = inner_product(data[id], current_query, D);
    return sim;
  }
  vector<operand> search(float* current_query){
    vector<bool> visited(N, false);
    vector<bool> inserted(N, false);
    StreamTracker st(16);
    //int magic = 4;
    // auto t1 = Clock::now();
    walk_in_graph(seed, current_query, inserted, visited, st);
    // auto t2 = Clock::now();
    partial_sort(top_k_candidates.begin(), top_k_candidates.begin() + K , top_k_candidates.end(), compare_descending);
    // auto t3 = Clock::now();
    // timers[1] += compute_clock(t2, t1);
    // timers[2] += compute_clock(t3, t2);
    return top_k_candidates;
  }
  void inline insertCandidate(vector<bool> &visited, int id, float sim) { ///, float* current_query) {
    if(!visited[id]) {
     visited[id] = true;
     top_k_candidates.emplace_back(id, sim);
    }
  }                                        
  void walk_in_graph(int seed, float* current_query, vector<bool> &inserted, vector<bool>& visited, StreamTracker& st) {
    int curNode = seed;
    float curSim = getSimilarity(curNode, current_query);
    st.update(curNode, curSim);
    while(st.size > 0 && searched < nSearch) {
      operand c = st.get();
      curNode = c.first;
      curSim = c.second;
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
      bool skip = false;
      int deg = min(load_front_M, neighborSize);
      if (neighbor.size() <= skiplimit) {
	      //deg =  neighbor.size();
        deg = min((int)neighbor.size(), skipdegree);
        skip = true;
      }
      for (int i = 0; i < deg; i++) {
        searched++;
        float sim = getSimilarity(neighbor[i], current_query);
        bool cond1 = sim > curNeighborSim;
        bool cond2 = sim > secondLargestNeighborSim;
        if((cond1 || cond2) && !visited[neighbor[i]]) {
          if(cond1) {
            secondLargestNeighborSim = curNeighborSim; 
            secondLargestNeighbor = curNeighborNode;
            curNeighborSim = sim;
            curNeighborNode = neighbor[i];
          }
          else if(cond2) {
            secondLargestNeighbor = neighbor[i];
            secondLargestNeighborSim = sim;
          }
        }
      }
      
      if(secondLargestNeighbor != -1)
        st.update(secondLargestNeighbor, secondLargestNeighborSim);

      if (curNeighborNode == -1)
        continue;

      st.update(curNeighborNode, curNeighborSim);
      skip = true;
      if(skip) {
        skipwalk++;
        continue;
      }
      realwalk++;
      secondLargestNeighbor = -1;
      secondLargestNeighborSim = numeric_limits<float>::lowest();
        
      while(cont) {
        cont = false;
        in_walk++;  
        vector<int> &neighborList = (*in_graph)[curNode][curNeighborNode];
        // deg = (int)min(64, (int)neighborList.size());
        deg = neighborList.size();
        for(int i = 0; i < deg; i++) {
          itc2++;
          int neighborNode = neighborList[i];
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
      st.update(curNeighborNode, curNeighborSim);
      if(secondLargestNeighbor != -1)
        st.update(secondLargestNeighbor, secondLargestNeighborSim);
      insertCandidate(inserted, curNeighborNode, curNeighborSim);
    }
  }

  void print_time() {
    // cout << "[" << name << "] " << (timers[1] + timers[2]) << "(ms)\n";
    // cout << "Preprocessing : " << timers[0] << "(ms)\n";
    // cout << "Search " << (timers[1] + timers[2]) << "(ms)\n";
    // cout << "---[1] Walk Front : " << timers[1] << "(ms)\n";
    // cout << "---[2] Walk : " << timers[2] << "(ms)\n";
    // cout << "Average Large walk : " << (float)walk << "\n";
    // cout << "Average Real walk : " << (float)realwalk << "\n";
    // cout << "Average steps in IG walk : " << (float)in_walk / realwalk << "\n";
    // cout << "Skip Walk : "<< (float)skipwalk / walk << "\n";
    // cout << "Per Node In-walk : "<< (float)ct2 / in_walk << "\n";
    // cout << "Average Walk Spent on Out : " << (float)ct1 << "\n";
    // cout << "Average Walk Spent on In : "  << (float)ct2 << "\n";
   }
  void print_information() {
    // cout << "[GraphSearchYJBase] \n";
  }
}; 
