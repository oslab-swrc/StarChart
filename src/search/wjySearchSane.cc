using namespace kx;
class SearchWJYSane : public BaseSearch {
public:
  float** unitvec;
  unordered_map<int, vector<int>> gx;
  float cos_threshold;
  int candidate_size;
  int S; 

  SearchWJYSane(string name, int k, int nSearch) : BaseSearch(name, k, nSearch) {
    S = 1200; //Fix
    //Well work for 100000_100_cnormal
    cos_threshold=0.142;
    candidate_size=2000;
  }
  ~SearchWJYSane() {
    for (int i = 0; i < S; i++)
      delete[] unitvec[i];
    delete[] unitvec;
  }
  void param_register(vector<int> params) {
  }
  void preprocess(){
    auto t1 = Clock::now();

    random_device rd;
    mt19937 generator(rd());
    normal_distribution<> d1{0,1};
    unitvec = new float*[S];
    gx = unordered_map<int, vector<int>>();
    for (int i = 0; i<S; i++) {
      gx.insert(make_pair(i, vector<int>()));
      unitvec[i] = new float[D];
      for (int j = 0; j<D; j++) {
        unitvec[i][j] = d1(generator);
      }
      float norm = sqrt(squared_sum(unitvec[i], D));
      for (int j = 0; j<D; j++) {
        unitvec[i][j] /= norm;
      }
    }

    #pragma omp parallel for
    for (int i = 0; i<S; i++) {
      vector<operand> op;
      for (int j = 0; j<N; j++) {
        if (getSimilarity(j,unitvec[i]) > cos_threshold)
          gx[i].push_back(j);
      }
    }

    auto t2 = Clock::now();
    timers[0] += compute_clock(t2, t1);
    ready = true;
  }
  float inline getSimilarity(int id, float *current_query) {
    return inner_product(data[id], current_query, D);
  }

  vector<operand> search(float* current_query){
    assert(ready == true);
    int dirth = 92;
    flat_hash_map<int, int> score;
    vector<int> candidate;
    vector<operand> dirs;
    auto t1 = Clock::now();
    for (int i=0; i<S; i++) {
      dirs.push_back(make_pair(i, inner_product(current_query, unitvec[i], D)));
    }
    partial_sort(dirs.begin(), dirs.begin() + dirth, dirs.end(), compare_descending);

    for (int i=0; i<dirth; i++) {
      int dir = dirs[i].first;
      for (int n : gx[dir]) {
        if(score.find(n) == score.end())
          score.insert(make_pair(n, 1));
        else
          score[n]++;
      }
    }

    auto t2 = Clock::now();

    vector<intoperand> nums;
    for (auto it=score.begin(); it!=score.end(); ++it) {
      nums.push_back(make_pair(it->first, it->second));
    }
    //radix_sort(nums.begin(), nums.end(), radixsort_descending());
   
    auto t3 = Clock::now();
    int sizes = nums.size() < candidate_size ? nums.size() : candidate_size;
    partial_sort(nums.begin(), nums.begin() + sizes, nums.end(), int_compare_descending);
    
    //int reverse = nums.size()-1;
    //int until = reverse - sizes;
    // for (int i = reverse; i>until; i--) {
    //   top_k_candidates.push_back(make_pair(nums[i].first, getSimilarity(nums[i].first, current_query)));
    // }

    auto t4 = Clock::now();
    //cout << (int)nums[0].second << " " << (int)nums[sizes-1].second << endl;
    for(int i=0; i<sizes; i++)
      top_k_candidates.push_back(make_pair(nums[i].first, getSimilarity(nums[i].first, current_query)));
    partial_sort(top_k_candidates.begin(), top_k_candidates.begin()+K, top_k_candidates.end(), compare_descending);
    auto t5 = Clock::now();
    timers[1] += compute_clock(t2, t1);
    timers[2] += compute_clock(t3, t2);
    timers[3] += compute_clock(t4, t3);
    timers[4] += compute_clock(t5, t4);
    return top_k_candidates;
  }

  void print_time() {
    cout << "[SearchWJYSane] " << (timers[1] + timers[2] + timers[3]+timers[4]) /queryCount << "(ms)\n";
    cout << "Preprocessing : " << timers[0] << "(ms)\n";
    cout << "Search " << (timers[1] + timers[2] + timers[3]+timers[4]) / queryCount << "(ms)\n";
    cout << "---[1] S unitvec Sampling : " << timers[1] / queryCount << "(ms)\n";
    cout << "---[2] Sample Sorting : " << timers[2] / queryCount << "(ms)\n";
    cout << "---[3] Candidate innerproduct : " << timers[3] / queryCount << "(ms)\n";
    cout << "---[4] Candidate innerproduct : " << timers[4] / queryCount << "(ms)\n";
  }
  void print_information() {

  }
};


