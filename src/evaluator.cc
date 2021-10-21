class Evaluator {
public:
  string name;  
  float** data;
  float** query;
  int N, D, K;
  int nSearch;
  float timers[10];
  vector<float> recall;
  vector<int> queryCount;
  map<string, int> evaluator_map;
  vector<vector<operand>> true_top_k;
  Evaluator(string name_, float** data_, float** query_, int n, int d, int k, int nSearch_) : name(name_), data(data_), query(query_), N(n), D(d), K(k), nSearch(nSearch_) {
  }
  ~Evaluator() {

  }
  void run(string path, int qN) {
    for(int i = 0; i < qN; i++) {
      true_top_k.push_back(vector<operand>());
      true_top_k[i].reserve(N);
    }
    bool loadSuccess = false;
    if (path != "") { // load mode
      loadSuccess = load_groundtruth(qN, path);
    }
    auto t1 = Clock::now();
    if (path == "" || !loadSuccess) {   // generate mode
      #pragma omp parallel for num_threads(OMD)
      for (int j = 0; j < qN; j++) {
          float similarity;
          vector<operand> vec_private;
          vec_private.reserve(N);
          for(int i = 0; i < N; i++) {
            similarity = inner_product(data[i], query[j], D);
            vec_private.push_back(make_pair(i,similarity));
          }
          int sortrange = min(500, N);
          partial_sort(vec_private.begin(), vec_private.begin() + sortrange, vec_private.end(), compare_descending);
          for(int k = 0; k <sortrange; k++) {
            true_top_k[j].push_back(vec_private[k]);
          }
      }
      auto t2 = Clock::now();
      timers[0] += compute_clock(t2, t1);
    }
    else {
      auto t2 = Clock::now();
      timers[0] += compute_clock(t2, t1);
    }
  }
  bool load_groundtruth(int qN, string path) {
    fstream gtfile;
    gtfile.open(path, fstream::in);
    if (!gtfile.is_open() || gtfile.fail()) {
      cout << "Groundtruth file doesnt exist, it will be generated : " << path << endl;
      return false;
    }
    int id = 0;
    int idx;
    float sim;
    string line;
    for (int j = 0; j < qN; j++) {
      getline(gtfile, line);
      istringstream iss(line);
      for (int i = 0; i < K; i++){
        iss >> idx;
        iss >> sim;
        true_top_k[id].push_back(make_pair(idx, sim));
      }
      id++;
    }
    gtfile.close();
    return true;
  }
  void save_groundtruth(vector<vector<operand>> c, string path) {
    cout << "Saving groundtruth result to " << path << "\n";
    fstream gtfile;
    gtfile.open(path, fstream::out);
    for (auto it = c.begin(); it != c.end(); it++) {
      for (int i = 0; i < 500; i++) {
        operand f = (*it)[i];
        gtfile << f.first << " " << f.second << " ";
      }
      gtfile << "\n";
    }
    gtfile.close();
    cout << "Finished saving groudntruth result\n";
  }
  void set_evaluator(string name) {
    evaluator_map.insert(make_pair(name, evaluator_map.size()));
    recall.push_back(0.0f);
    queryCount.push_back(0);
    assert(evaluator_map.size() == recall.size() && recall.size() == queryCount.size());
  }
  int get_evaluator(string name) {
    assert(evaluator_map.find(name) != evaluator_map.end());
    return evaluator_map.find(name)->second;
  }
  //#define EXP
  #ifndef EXP
  float check_recall(vector<operand> true_top_k, vector<operand> top_k_candidates, bool silence=false, bool stat_mode=false) {
    int count = 0;
    size_t check = stat_mode ? top_k_candidates.size() : K;
    if (top_k_candidates.size() < (size_t)K ) {
      if (!silence) {
        cout << name  << " - ";
        cout << "[WARNING] Less than K (" << K << ") search results (" << top_k_candidates.size() << ")\n";
      }
      check = top_k_candidates.size();
    }
    // cout << "K : " << K << " check : " << check << " true.size : " << true_top_k.size() << " top_k : " << top_k_candidates.size() << "\n";
    for (int i = 0; i < K ; i++) {
      // cout << "truetopk : " << true_top_k[i].first << "\n";
      // cout << "top kcan : " << top_k_candidates[i].first << "\n";
      for (size_t j = 0; j < check; j++) {
        if (true_top_k[i].first == top_k_candidates[j].first) {
          count++;
          break;
        }
      }
    }
    // cout << "count : " << count << "\n";
    return (float)count / (float)K * 100.0;
  }
  #else
  float check_recall(vector<operand> true_top_k, vector<operand> top_k_candidates, bool silence=false, bool stat_mode=false) {
    int count = 0;
    size_t check = stat_mode ? top_k_candidates.size() : K;
    if (top_k_candidates.size() < (size_t)K ) {
      // if (!silence) {
      //   cout << name  << " - ";
      //   cout << "[WARNING] Less than K (" << K << ") search results (" << top_k_candidates.size() << ")\n";
      // }
      check = top_k_candidates.size();
    }
    // cout << "K : " << K << " check : " << check << " true.size : " << true_top_k.size() << " top_k : " << top_k_candidates.size() << "\n";
    for (int i = 0; i < K ; i++) {
      // cout << "truetopk : " << true_top_k[i].first << "\n";
      // cout << "top kcan : " << top_k_candidates[i].first << "\n";
      for (size_t j = 0; j < check; j++) {
        if (true_top_k[i].first == top_k_candidates[j].first) {
          count++;
          break;
        }
      }
    }
    if(count > 0)
      count = K;
    // cout << "count : " << count << "\n";
    return (float)count / (float)K * 100.0;
  }
  #endif


  float update(int evaluator_id, vector<operand> top_k_candidates) {
    float recall_tmp = check_recall(true_top_k[evaluator_id], top_k_candidates);
    // cout << "K " << K << " queryCount : " << queryCount[evaluator_id] << " recall : " << recall_tmp << "\n";
    recall[evaluator_id] += recall_tmp;
    queryCount[evaluator_id]++;
    return recall_tmp;
  }
  void clear() {
    recall.clear();
  }
  void print_result() {
    assert(recall.size() == queryCount.size());
    for (auto it = evaluator_map.begin(); it != evaluator_map.end(); it++) {
      int id = get_evaluator((*it).first);
      cout << (*it).first << " " << name << " : " << recall[id] / queryCount[id] << "(%)\n";
    }
  }
};
