#ifdef FIGURE
extern int ef_ = 100;
#endif
class HnswSearch : public BaseSearch {
public:
  hnswlib::HierarchicalNSW<float> *hnswg;
  int M;
  HnswSearch(string name, int k, int nSearch) : BaseSearch(name, k, nSearch) {
    hnsw_nSearch = nSearch; 
    #ifdef FIGURE
    hnsw_nSearch = N;
    #endif
  }
  ~HnswSearch() {
  }
  void param_register(vector<int> params) {
    cout << "[Parameters]\n--" << name << "\n";
    assert(params.size() == 1);   // there's only M
    for (size_t i = 0; i < params.size(); i++) {
      if (i == 0) {
        M = params[i];
        cout << "M : " << M << "\n";
      }
    }
  }
  void preprocess(){
    auto t1 = Clock::now();
    hnswg->setEf(200);
    hnswg->dist_calc=0;
    total_dist_calc=0;
    auto t2 = Clock::now();

    for (int i = 0; i < 10; i++)
      timers[i] = 0;
		searched = 0;
    totalSearched = 0;
    queryCount = 0;
    #ifdef FIGURE
    hnswg->setEf(ef_);
    #endif

    timers[0] += compute_clock(t2, t1);
    ready = true;
  }
  float inline getSimilarity(int id, float *current_query) {
    return inner_product(data[id], current_query, D);
  }

  vector<operand> search(float* current_query){
    assert(ready == true);
    auto t1 = Clock::now();
    priority_queue<pair<float, size_t>> gt = hnswg->searchKnn(current_query, K);
    while(!gt.empty()) {
      top_k_candidates.push_back(make_pair(gt.top().second, gt.top().first));
      gt.pop();
    }
    auto t2 = Clock::now();
    timers[1] += compute_clock(t2, t1);
    searched = 0;
    totalSearched = total_dist_calc ;
    return top_k_candidates;
  }

  void print_time() {
    cout << "[HnswSearch] " << (timers[1] ) /queryCount << "(ms)\n";
    cout << "Preprocessing : " << timers[0] << "(ms)\n";
    cout << "Search " << (timers[1]) / queryCount << "(ms)\n";
   }
  void print_information() {
    cout << "[HnswSearch] \n";
  }
}; 
