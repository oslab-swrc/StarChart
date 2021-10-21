class BaseSearch {
public:
  static float **data;
  static int N, D, K;
  string name;  
  int nSearch;
  float timers[12];
  int searched;
  int totalSearched;
  int queryCount;
  bool ready;
  vector<operand> top_k_candidates;
  float searchTime;
  vector<operand> buffer;
  BaseSearch(string name, int k, int nSearch) : name(name), nSearch(nSearch) {
    for (int i = 0; i < 12; i++)
      timers[i] = 0;
		searched = 0;
    totalSearched = 0;
    queryCount = 0;
    ready = false;
    searchTime = 0.0f;
    top_k_candidates.reserve(nSearch);
    buffer.reserve(nSearch);
  }
  ~BaseSearch() {
  }
  static void init(int n, int d, int k, float** data_in) {
    N = n;
    D = d;
    K = k;
    data = data_in;
  }
  virtual void param_register(vector<int>) = 0;
  virtual vector<operand> search(float *current_query) = 0;
  virtual void preprocess() = 0;
  virtual void update() {
    totalSearched += searched;
    queryCount++;
    top_k_candidates.clear();
    buffer.clear();
    searched = 0;
  }
  virtual void print_time() = 0;
  virtual void print_information() = 0;
  virtual void datasort() {}

  virtual void print_result() {
    cout << name << " : " << (float)totalSearched / queryCount << "\n";
  }
  void print_output() {
    datasort();
    int printK = K > 20 ? 20 : K;
    cout << name << " Result ( total : " << top_k_candidates.size() << ") \n";
    for (size_t i = 0; i < (size_t)printK && i < top_k_candidates.size(); i++)
      cout << "( " << top_k_candidates[i].first << " , " << top_k_candidates[i].second << " )\n";
    cout << "\n";
  }
  virtual void clear(){
    for (int i = 0; i < 10; i++)
      timers[i] = 0;
    searched = 0;
    totalSearched = 0;
    queryCount = 0;
  }
};
float** BaseSearch::data = NULL;
int BaseSearch::N = -1;
int BaseSearch::D = -1;
int BaseSearch::K = -1;