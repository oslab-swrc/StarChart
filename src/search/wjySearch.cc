/*
using phmap::flat_hash_map;
using namespace kx;
typedef int diff_t;
typedef pair<int, pair<int, int>> ppi;
class SearchWJY : public BaseSearch {
public:
  float** unitvec;
  vector<diff_t*> gx;
  vector<int> gs;

  float cos_threshold;
  uint8_t visit_threshold;
  int candidate_size;
  int S; 

  SearchWJY(string name, int k, int nSearch) : BaseSearch(name, k, nSearch) {
    S = 1600; //Fix
    
    //Well work for 1000000_100_cnormal
    cos_threshold=0.19;
    candidate_size=1000;
    visit_threshold = 8;
  }
  ~SearchWJY() {
  }
  void param_register(vector<int> params) {
  }
  void preprocess(){
    auto t1 = Clock::now();
    int max = 0;

    random_device rd;
    mt19937 generator(rd());
    normal_distribution<> d1{0,1};
    unitvec = new float*[S];
    vector<int> duplicate(N);
    for (int i = 0; i<S; i++) {
      unitvec[i] = new float[D];
      gx.push_back((diff_t*)&i);
      gs.push_back(i);
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
      vector<int> op;
      for (int j = 0; j<N; j++) {
        if (getSimilarity(j,unitvec[i]) > 0.13)
          op.push_back(j);
      }
      sort(op.begin(), op.end());

      #pragma omp critical 
      {
        //int __attribute__((aligned(64))) x[op.size()];
        int opsize = op.size() * sizeof(diff_t); // *4
        diff_t *x = (diff_t*)malloc(opsize);
        int idx=0;

        for (int a : op) {
          x[idx++] = a;
        }
   
        gx[i] = x;
      }

      gs[i]=op.size();
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
    uint8_t visited[N] = {0};
    vector<int> candidate;
    vector<intoperand> nums;
    priority_queue<ppi, vector<ppi>, greater<ppi>> pq;

    auto t1 = Clock::now();
    int min = 1000000;
    for (int i=0; i<S; i++) {
      if (inner_product(current_query, unitvec[i],D) > cos_threshold) {
        pq.push({gx[i][0], {i, 0}});
      }
    }
    auto t2 = Clock::now();
    
    uint8_t cnt = 0;
    int curr_val = pq.top().first;
    while (!pq.empty()) {
      ppi curr = pq.top(); pq.pop();

      int i = curr.second.first;
      int j = curr.second.second;
      int now = curr.first;
    
      if (now == curr_val) {
        cnt++;
      }
      else {
        if (cnt > visit_threshold)
          nums.emplace_back(make_pair(curr_val, cnt));
        curr_val = now;
        cnt = 1;
      }

      int nextj = j+1;
      if (nextj < gs[i]) {
        pq.push({ gx[i][nextj], {i, nextj}});
      }
    }

    auto t3 = Clock::now();
    auto t4 = Clock::now();

    int sizes = nums.size() < candidate_size ? nums.size() : candidate_size;
    radix_sort(nums.begin(), nums.end(), radixsort_descending());
    auto t5 = Clock::now();
    int reverse = nums.size()-1;
    int until = reverse - sizes;
    for (int i = reverse; i>until; i--) {
      top_k_candidates.emplace_back(make_pair(nums[i].first, getSimilarity(nums[i].first, current_query)));
    }
    
    partial_sort(top_k_candidates.begin(), top_k_candidates.begin()+K, top_k_candidates.end(), compare_descending);

    auto t6 = Clock::now();
    timers[1] += compute_clock(t2, t1);
    timers[2] += compute_clock(t3, t2);
    timers[3] += compute_clock(t4, t3);
    timers[4] += compute_clock(t5, t4);
    timers[5] += compute_clock(t6, t5);
    return top_k_candidates;
  }

  void print_time() {
    cout << "[SearchWJY] " << (timers[1] + timers[2] + timers[3]+timers[4]+timers[5]) /queryCount << "(ms)\n";
    cout << "Preprocessing : " << timers[0] << "(ms)\n";
    cout << "Search " << (timers[1] + timers[2] + timers[3]+timers[4]+timers[5]) / queryCount << "(ms)\n";
    cout << "---[1] S unitvec Sampling : " << timers[1] / queryCount << "(ms)\n";
    cout << "---[2] Sample Sorting : " << timers[2] / queryCount << "(ms)\n";
    cout << "---[3] Candidate innerproduct : " << timers[3] / queryCount << "(ms)\n";
    cout << "---[4] Candidate innerproduct : " << timers[4] / queryCount << "(ms)\n";
    cout << "---[4] Candidate innerproduct : " << timers[5] / queryCount << "(ms)\n";
  }
  void print_information() {

  }
};

*/


using phmap::flat_hash_map;
using namespace kx;
typedef uint16_t diff_t; //1. uint8_t possible when normal 
                         //2. bitmap / invlist compression will much faster

class SearchWJY : public BaseSearch {
public:
  float** unitvec;
  vector<diff_t*> gx;
  vector<int> gs;
  vector<int> start;

  //It can be Further Optimize 
  const float cos_threshold = 0.245; //OPT1 : theta1 < theta2
  const uint8_t visit_threshold = 11;
  const int candidate_size = 3000;
  const int S = 8000; 

  int ctstat;

  SearchWJY(string name, int k, int nSearch) : BaseSearch(name, k, nSearch) {
    ctstat = 0;
  }
  ~SearchWJY() {
    for (int i = 0; i < S; i++)
      delete[] unitvec[i];
    delete[] unitvec;

  }
  void param_register(vector<int> params) {
  }
  void preprocess(){
    auto t1 = Clock::now();
    int max = 0;

    random_device rd;
    mt19937 generator(rd());
    normal_distribution<> d1{0,1};
    unitvec = new float*[S];
    vector<int> duplicate(N);
    for (int i = 0; i<S; i++) {
      unitvec[i] = new float[D];
      gx.push_back((diff_t*)&i);
      gs.push_back(i);
      start.push_back(i);

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
      vector<int> op;
      for (int j = 0; j<N; j++) {
        if (getSimilarity(j,unitvec[i]) > 0.14)
          op.push_back(j);
      }
      sort(op.begin(), op.end());

      int opsize = op.size() * sizeof(diff_t); // *4

      diff_t *x = (diff_t*)malloc(opsize);
      int idx=0;
      bool first = true;
      int prev = 0;

      for (int a : op) {
        if (first) {
          start[i] = a;
          first = false;
        }
        else {
          x[idx++] = a - prev;
          if (a-prev > max)
            max = a-prev;
        }
        prev = a;
      }
 
      gx[i] = x;

      gs[i]=op.size()-1;
    }

    cout << "MAX : " << max << endl;

    auto t2 = Clock::now();
    timers[0] += compute_clock(t2, t1);
    ready = true;
  }
  float inline getSimilarity(int id, float *current_query) {
    return inner_product(data[id], current_query, D);
  }
  vector<operand> search(float* current_query){
    assert(ready == true);
    uint8_t visited[N];
    memset(visited, 0, sizeof(uint8_t)*N);
    vector<int> candidate;
    vector<intoperand> nums(N);

    auto t1 = Clock::now();
    for (int i=0; i<S; i++) {
      if (inner_product(current_query, unitvec[i],D) > cos_threshold) {
        candidate.push_back(i);
      }
    }
    auto t2 = Clock::now();

    //Further optimization : 1. Bitmap un/compressed scancount
    //                       2. Bitmap un/compressed RBmrg
    for (int& i : candidate) {
      int size = gs[i];
      diff_t* g = gx[i];
      int first = start[i];

      for (int n =0; n<size; n++) {
        first += g[n];
        visited[first]++;
      }
    }

    auto t3 = Clock::now();

    int idx = 0;
    for (int i = 0; i<N; i++) {
      uint8_t& visits = visited[i];
      if (visits > visit_threshold)
        nums[idx++]=make_pair(i, visits);
    }
    nums.resize(idx);
    auto t4 = Clock::now();

    int sizes = nums.size() < candidate_size ? nums.size() : candidate_size;
    radix_sort(nums.begin(), nums.end(), radixsort_descending());
    auto t5 = Clock::now();
    int reverse = nums.size()-1;
    int until = reverse - sizes;
    for (int i = reverse; i>until; i--) {
      top_k_candidates.emplace_back(make_pair(nums[i].first, getSimilarity(nums[i].first, current_query)));
    }
     
    partial_sort(top_k_candidates.begin(), top_k_candidates.begin()+K, top_k_candidates.end(), compare_descending);

    auto t6 = Clock::now();
    timers[1] += compute_clock(t2, t1);
    timers[2] += compute_clock(t3, t2);
    timers[3] += compute_clock(t4, t3);
    timers[4] += compute_clock(t5, t4);
    timers[5] += compute_clock(t6, t5);
    return top_k_candidates;
  }

  void print_time() {
    cout << "[SearchWJY] " << (timers[1] + timers[2] + timers[3]+timers[4]+timers[5]) /queryCount << "(ms)\n";
    cout << "Preprocessing : " << timers[0] << "(ms)\n";
    cout << "Search " << (timers[1] + timers[2] + timers[3]+timers[4]+timers[5]) / queryCount << "(ms)\n";
    cout << "---[1] S unitvec Sampling : " << timers[1] / queryCount << "(ms)\n";
    cout << "---[2] Sample Sorting : " << timers[2] / queryCount << "(ms)\n";
    cout << "---[3] Candidate innerproduct : " << timers[3] / queryCount << "(ms)\n";
    cout << "---[4] Candidate innerproduct : " << timers[4] / queryCount << "(ms)\n";
    cout << "---[4] Candidate innerproduct : " << timers[5] / queryCount << "(ms)\n";
    //cout << "C : " << (float)sizeavg / queryCount << "\n";
  }
  void print_information() {

  }
};

