using phmap::flat_hash_map;
using phmap::flat_hash_set;
class GraphSearchTJRand : public BaseSearch {
public:
  Graph* fg;
  int M;
  int walkc;
  vector<int> rank;
  GraphSearchTJRand(string name, int k, int nSearch) : BaseSearch(name, k, nSearch) {
    fg = new Graph();
    walkc = 0;
  }
  ~GraphSearchTJRand() {
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

  //#define LOG
  vector<operand> search(float* current_query){
    assert(ready == true);
    //#ifdef LOG
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
    //#endif
    vector<bool> visited(N, false);
    vector<bool> inserted(N, false);
    auto t1 = Clock::now();
    //int init = rand() % N;
    int init = 0;
    for(int i = 0; i < 2048; i++) {
      walk(init, visited, inserted, current_query);
      walkc++;
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
    int nextSeed = -1;
    //int obufSize = 0;
    //cout << "New Walk *** \n";
    while(found) {
      visited[curSeed] = true;
      //cout << "Node : " << rank[curSeed] << " / " << (int)visited[curSeed] << "\n";
      insertCandidate(inserted, curSeed, curSim);
      found = false;
      //vector<operand> obuf;
      if(rank[curSeed] < 10)
        searchEdge = 1024;
      else if(rank[curSeed] < 40)
        searchEdge = 512;
      //obuf.reserve(searchEdge);

      vector<int> &neighborList = (*fg->graph)[curSeed];
      float neighborMaxSim = numeric_limits<float>::lowest();
      int maxNeighbor = -1;

      for (int i=0; i<searchEdge; i++) {
        //int idx = i + ((int)visited[curSeed]-1) * searchEdge;
        int candidate = neighborList[i];
        float neighborSim = getSimilarity(candidate, current_query);
        searched++;
        //insertCandidate(inserted, candidate, neighborSim);
        //if(i < obufSize)
        //obuf.emplace_back(candidate, neighborSim);
        if(neighborSim > curSim) {
          if(visited[candidate])
            continue;
          curSim = neighborSim;
          nextSeed = candidate;
          found = true;
        }
        if(neighborSim > neighborMaxSim) {
          if(visited[candidate])
            continue;
          neighborMaxSim = neighborSim;
          maxNeighbor = candidate;
        }
        if (found && i>=16)
          break;
        //if(found && rank[curSeed] >= 25)
        //  break;
      }
      if(!found) {
        //for(int i=0; i<obuf.size(); i++) 
        //  insertCandidate(inserted, obuf[i].first, obuf[i].second);
        seed = maxNeighbor;
      }
      else {
        curSeed = nextSeed;
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
    cout << "---[*] Walk: " << walkc/queryCount << "\n";

    //cout << "---[2-2] Expand-Sort : " << timers[4] / queryCount << "(ms)\n";
   }
  void print_information() {
    cout << "[GraphSearchTJ] \n";
  }
}; 
