using phmap::flat_hash_map;
using phmap::flat_hash_set;
class GraphSearchYJ3 : public BaseSearch {
public:
  unordered_map<int, pair<bool, vector<int>>> bff;
  unordered_map<int, int> friends;
  map<int, vector<float>> mids;

  unordered_map<int, unordered_map<int, vector<int>>>* in_graph;
  static int cnt;
  static int cntcnt;
  static int walk;
  Graph* fg;
  int M;
  GraphSearchYJ3(string name, int k, int nSearch) : BaseSearch(name, k, nSearch) {
    fg = new Graph();
  }
  ~GraphSearchYJ3() {
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
    // vector<operand> answer;
    // vector<int> rank;
    // rank.reserve(N);
    // answer.reserve(N);
    // for(int i = 0; i < N; i++) {
    //   float sim = getSimilarity(i, current_query);
    //   answer.push_back(make_pair(i, sim));  
    //   rank.push_back(-1);
    // }
    // sort(answer.begin(), answer.end(), compare_descending);
    // for(int i = 0; i < N; i++)
    //   rank[answer[i].first] = i;
    // friendAnalysis();
    // assert(ready == true);
    flat_hash_set<int> visited;
    flat_hash_set<int> inserted;
    inserted.reserve(4096);
    visited.reserve(1024);
    auto t1 = Clock::now();
    // int init = rand() % N;
    // int walkcount = 0;
    // for(int i = 0; i < 128; i++) {
    //   walk_bff(init, visited, inserted, current_query);
    //   if(searched > nSearch) {
    //     walkcount = i;
    //     break;
    //   }
    // }
    int seed = 0;
    walk_in_graph(seed, current_query, inserted, visited);
    auto t2 = Clock::now();
    partial_sort(top_k_candidates.begin(), top_k_candidates.begin() + K, top_k_candidates.end(), compare_descending);
    auto t3 = Clock::now();
    
    // // for debugging
    // bool check = true;
    // for (int i = 0; i < top_k_candidates.size(); i++) {
    //   bool found = false;
    //   int index = top_k_candidates[i].first;
    //   if (bff[index].second.size() > 1024 || bff[index].second.size() == 1) {
    //     int bestFriend = friends[index];
    //     for (int j = 0; j < top_k_candidates.size(); j++) {
    //       if (top_k_candidates[j].first == bestFriend)
    //         found = true;
    //     }
    //     if (!found) {
    //       check = false;
    //       assert(false);
    //     }
    //   }
    // }
    timers[1] += compute_clock(t2, t1);
    timers[2] += compute_clock(t3, t2);
    return top_k_candidates;
  }
  float computeMids(int idx1, int idx2, float* current_query) {
    float* mid = new float[D];
    for (int i = 0; i < D; i++)
      mid[i] = (data[idx1][i] + data[idx2][i]) / 2;
    return inner_product(mid, current_query, D);
  }
  void inline insertCandidate(flat_hash_set<int> &visited, int id, float sim) { ///, float* current_query) {
    if(visited.find(id) == visited.end()) {
      visited.insert(id);
      top_k_candidates.emplace_back(id, sim);
    }
  }

  void friendAnalysis() {
    for (int i = 0; i < bff.size(); i++) {
      if (bff[i].first == false) {
        friends.insert(make_pair(i, bff[i].second[0]));
        friends.insert(make_pair(bff[i].second[0], i));
      }
    }
    cout << "Number of friends : " << friends.size() << "\n";
  }
  void walk_bff(int &seed, flat_hash_set<int> &visited, flat_hash_set<int> &inserted, float *current_query) {
    int curSeed = seed;
    float curSim;
    if (bff[curSeed].first == false) {
      assert(bff[curSeed].second.size() == 1);
      curSeed = bff[curSeed].second[0];
    }
    // if (friends.find(curSeed) == friends.end())
    //   curSim  = getSimilarity(curSeed, current_query);
    // else
    //   curSim = computeMids(curSeed, friends[curSeed], current_query);
    curSim = getSimilarity(curSeed, current_query);
    assert(curSeed != -1);
    bool found = true;
    int searchEdge = 512;
    int obufSize = 64;
    int pairInsertion = 0;
    while(found) {
      assert(bff[curSeed].first != false);
      visited.insert(curSeed);
      if (friends.find(curSeed) != friends.end()) {
        int bestFriend = friends[curSeed];
        float bestFriendSim = getSimilarity(bestFriend, current_query);
        // float bestFriendSim = curSim;
        assert(bff[bestFriend].second[0] == curSeed);
        insertCandidate(inserted, bestFriend, bestFriendSim);
        pairInsertion++;
        // cout << "Pair insertion1 : " << curSeed << " , "  << bestFriend << "\n";
      }
      insertCandidate(inserted, curSeed, curSim);
      found = false;
      vector<operand> obuf;
      obuf.reserve(obufSize);
      vector<int> &neighborList = bff[curSeed].second;
      assert(neighborList.size() != 1);
      float neighborMaxSim = numeric_limits<float>::lowest();
      int maxNeighbor = -1;
      for (int i = 0; i < searchEdge; i++) {
        int candidate = neighborList[i];
        if (bff[candidate].first == false) {
          assert(bff[candidate].second.size() == 1);
          candidate = bff[candidate].second[0];
        }
        float neighborSim = getSimilarity(candidate, current_query);
        // if (friends.find(candidate) == friends.end())
        //   neighborSim = getSimilarity(candidate, current_query);
        // else
        //   neighborSim = computeMids(candidate, friends[candidate], current_query);
        searched++;
        if(i < obufSize)
          obuf.emplace_back(candidate, neighborSim);
        if(neighborSim > curSim) {
          if(visited.find(candidate) != visited.end()) {
            continue;
          }
          curSim = neighborSim;
          curSeed = candidate;
          found = true;
        }
        if(neighborSim > neighborMaxSim) {
          if(visited.find(candidate) != visited.end())
            continue;
          neighborMaxSim = neighborSim;
          maxNeighbor = candidate;
        }
        if(found && i >= 8) {
          break;
        }
      }
      if(!found) {
        seed = maxNeighbor;
        for(int i = 0; i < obufSize; i++) {
          insertCandidate(inserted, obuf[i].first, obuf[i].second);
          if (friends.find(obuf[i].first) != friends.end()) {
            int bestFriend = friends[obuf[i].first];
            float bestFriendSim = getSimilarity(bestFriend, current_query);
            // float bestFriendSim = obuf[i].second;

            assert(bff[bestFriend].first == false);
            assert(bff[bestFriend].second[0] == obuf[i].first);
            insertCandidate(inserted, bestFriend, bestFriendSim);
            pairInsertion++;
            // cout << "Pair insertion2 : " << obuf[i].first << " ,  " << bestFriend << "\n";
          }
        }
      }
    }
    cout << "Total pair insertion : " << pairInsertion << " top_k_size : " << top_k_candidates.size() << "\n";
  }

  void walk_in_graph(int seed, float* current_query, flat_hash_set<int> &inserted, flat_hash_set<int>& visited) {
    int curSeed = seed;
    float curSim = getSimilarity(curSeed, current_query);
    visited.insert(curSeed);
    while(true) {
      walk++;
      insertCandidate(inserted, curSeed, curSim);
      visited.insert(curSeed);
      bool cont = true;
      int neighborSeed;

      auto t5 = Clock::now();

      vector<int>& neighbor = (*fg->graph)[curSeed];
      float neighborSim = numeric_limits<float>::lowest();

      for (int i = 0; i < 16; i++) {
        if (visited.find(neighbor[i]) == visited.end()) {
          searched++;
          float sim = getSimilarity(neighbor[i], current_query);
          if (sim > curSim && i > 7) {
            neighborSeed = neighbor[i];
            neighborSim = sim;
            break;
          }
          if (sim > neighborSim) {
            neighborSeed = neighbor[i];
            neighborSim = sim;
          }
        }
      }
      auto t6 = Clock::now();

      searched++;
      cntcnt++;
      while(cont) {
        cont = false;
        vector<int> &neighborList = (*in_graph)[curSeed][neighborSeed];
        insertCandidate(inserted, neighborSeed, neighborSim);
        cnt++;
        assert(neighborList.size() == 32);
        for (size_t i = 0; i <(int)neighborList.size(); i++) {
          if (visited.find(neighborList[i]) == visited.end()) {
            float tempSim = getSimilarity(neighborList[i], current_query);
            searched++;
            if (tempSim > curSim) {
              neighborSeed = neighborList[i];
              goto jump;
            }
            if(tempSim > neighborSim) {
              neighborSim = tempSim;
              neighborSeed = neighborList[i];
              cont = true;
              // break;
            }
            if (cont && i > 7) 
              break;
          }
        }
      }
      jump:
      auto t7 = Clock::now();

      if (searched >= nSearch)
        break;
      curSeed = neighborSeed;
      curSim = getSimilarity(curSeed, current_query);
      timers[3] += compute_clock(t6, t5);
      timers[4] += compute_clock(t7, t6);

    }
    // walk += top_k_candidates.size();

  }
  void print_time() {
    cout << "[GraphSearchYJ] " << (timers[1] + timers[2]) /queryCount << "(ms)\n";
    cout << "Preprocessing : " << timers[0] << "(ms)\n";
    cout << "Search " << (timers[1] + timers[2]) / queryCount << "(ms)\n";
    cout << "---[1] Walk Front : " << timers[1] / queryCount << "(ms)\n";
    cout << "---[2] Walk : " << timers[2] / queryCount << "(ms)\n";
    cout << "---[3] Loop 1 : " << timers[3] / queryCount << "(ms)\n";
    cout << "---[4] Loop 2: " << timers[4] / queryCount << "(ms)\n";
    cout << "Average in_graph walk : " << (float)cnt / cntcnt << "\n";
    cout << "Average out walk : " << (float)walk / queryCount << "\n";
   }
  void print_information() {
    cout << "[GraphSearchTJ] \n";
  }
}; 
int GraphSearchYJ3::cnt = 0;
int GraphSearchYJ3::cntcnt = 0;
int GraphSearchYJ3::walk = 0;
