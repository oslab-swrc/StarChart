using Clock=::chrono::high_resolution_clock;

typedef pair<int, float> operand;
typedef pair<int, uint8_t> intoperand;

bool compare_ascending (operand first, operand second)
{ 
  if (first.second < second.second) 
    return true;
  else
    return false;
}

bool compare_descending_strict (operand first, operand second)
{ 
  if (first.second > second.second) 
    return true;
  else if(first.second == second.second && first.first > second.first) 
    return true;
  else
    return false;
}

bool compare_descending (operand first, operand second)
{ 
  if (first.second > second.second) 
    return true;
  else
    return false;
}

bool int_compare_descending (intoperand first, intoperand second)
{ 
  if (first.second > second.second) 
    return true;
  else
    return false;
}

bool rank_compare_ascending (pair<int, int> first, pair<int, int> second)
{ 
  if (first.second < second.second) 
    return true;
  else
    return false;
}

struct radixsort_descending {
  static const int nBytes = 1; //uint8_t 
  int kth_byte(const intoperand &x, int k) {
    return x.second;
  }
  bool compare(const intoperand &x, const intoperand &y) {
    return x.second < y.second;
  }
};

struct compare_priority_ascending {
    bool operator()(operand t, operand u){
        return t.second > u.second;
    }
};
struct compare_priority_descending {
    bool operator()(operand t, operand u){
        return t.second < u.second;
    }
};

typedef priority_queue <operand, vector<operand>, compare_priority_ascending> pq_a;
typedef priority_queue <operand, vector<operand>, compare_priority_descending> pq_d;

float inline compute_clock(chrono::system_clock::time_point t2, chrono::system_clock::time_point t1) {
  return (double)(chrono::duration_cast<chrono::microseconds>(t2 - t1).count()) / 1000;
}
void inline check_file(ofstream &ofile, string path) {
  if (!ofile.is_open()) {
    cout << "There is no such file ( " + path + " )\n";
    exit(0);
  }
}
void inline check_file(fstream &ofile, string path) {
  if (!ofile.is_open()) {
    cout << "There is no such file ( " + path + " )\n";
    exit(0);
  }
}

float inline squared_sum(float* __restrict__ in, int d) {
  float sum = 0.0f;
  #pragma omp simd reduction(+:sum)
  for (int i = 0; i < d; i++)
    sum += in[i] * in[i];
  return sum;
}
#define IP
float inline inner_product(float *__restrict__ a, float *__restrict__ b, int d, float anorm = 0.0f, float bnorm = 0.0f) {
  float similarity = 0.0f;
  #pragma omp simd reduction(+:similarity)
  #ifdef IP
  for (int m = 0; m < d; m++)
    similarity += a[m] * b[m];
  #else
  for (int m = 0; m < d; m++)
    similarity += (a[m] - b[m]) * (a[m] - b[m]);
  similarity = -sqrt(similarity);
  #endif
  return similarity;
}

float inline squared_euclidean_distance(float* __restrict__ a, float*__restrict__  b, int d) {
  float e_distance = 0.0f;
  #pragma omp simd reduction(+:e_distance)
  for (int i = 0; i < d; i++)
    e_distance += (a[i] - b[i]) * (a[i] - b[i]);
  return sqrt(e_distance);
}
float inline cosine_similarity(float *__restrict__ a, float *__restrict__ b, int d, float anorm, float bnorm) {
  float similarity = inner_product(a, b, d);
  float cossim = similarity / (anorm * bnorm);
  return cossim;
}
  
struct compare_priority_descending_tiebreak {
    bool operator()(operand t, operand u) {
        if(t.second < u.second)
          return true;
        else if(t.second == u.second && t.first < u.first)
          return true;
        else
          return false;
    }
};
struct dbl_pq_d { 
  set<operand, compare_priority_descending_tiebreak> s; 
  int size() {
    return s.size(); 
  }
  bool isEmpty() {
    return (s.size() == 0); 
  }
  void insert(operand x) { 
    s.insert(x); 
  }
  operand getMin() { 
    return *(s.begin()); 
  } 
  operand getMax() { 
    return *(s.rbegin()); 
  } 
  void deleteMin() { 
    if (s.size() == 0) 
      return; 
    s.erase(s.begin()); 
  } 
  void deleteMax() { 
    if (s.size() == 0) 
      return; 
    auto it = s.end(); 
    it--; 
    s.erase(it); 
  }
}; 
class DataManager {
public:
  string name;
  int N, D, qN;
  float **data;
  float **query;
  float timers[2];
  string data_path = "";
  string query_path= "";
  map<string, string> path_map;

  /* Data Processing */
  DataManager(string path, string data_name, string query_name, int q) : qN(q) {
    data_path = path + data_name + ".txt";
    query_path = path + query_name + "_query.txt";
    cout << "[Config] ./data/query path : " << data_path << " / " << query_path << "\n";
    //name = process_name(data_name);
    name = data_name;
    auto t1 = Clock::now();
    read_data(data_path, query_path);
    auto t2 = Clock::now();
    auto t3 = Clock::now();
    timers[0] = compute_clock(t2, t1);
    timers[1] = compute_clock(t3, t2);
  }
  ~DataManager() {
    #ifndef EXP
    cout << "[DataManager Information] \n";
    cout << "Data Reading : " << timers[0] << "(ms)\n";
    #endif
    for (int i=0; i < qN; i++)
      delete[] query[i];
    for (int i=0; i < N; i++)
      delete[] data[i];
    delete[] data;
    delete[] query;
  }
  string process_name (string name) {
    string processed_name = "";
    int ucount = 0;
    for (size_t i = 0; i < name.size(); i++) {
      if (ucount >= 2)
        processed_name += name[i];
      if (name[i] == '_')
        ucount++;
    }
    if (processed_name == "")
      processed_name = name;
    return processed_name;
  }
  string get_path(string name) {
    if(path_map.find(name) != path_map.end())
      return path_map.find(name)->second;
    else
      return "";
  }
  void set_paths(string name, string path) {
    path_map.insert(make_pair(name, path));
    cout << "[Config] " << name << " path : " << path << "\n";
  }
  void read_data(string data_path, string query_path) {
    fstream fs_data, fs_query;
    string line;
    string qline;
    cout << "Data path : " << data_path << " / query_path : " << query_path << endl;
    /* Open File */
    fs_data.open(data_path, fstream::in);
    fs_query.open(query_path, fstream::in);
    if (!fs_data.is_open() || !fs_query.is_open()){
      cout << "Please check data / query path.\n";
      exit(1);
    }
    getline(fs_data, line);
    getline(fs_query, qline);
    /* Read N & D */
    stringstream sts(line);
    stringstream qsts(qline);
    int q;
    sts >> N;
    sts >> D;
    qsts >> q;
    if (qN == -1)
      qN = q;
    int dummy = 0;
    qsts >> dummy;
    cout << "N : " << N << " / D : " << D << " / qN : " << qN << "\n";
    /* Initialize Structures */
    data = new float*[N];
    for (int i = 0; i < N; i++)
      data[i] = new float[D];
    query = new float*[q];
    for (int i = 0; i < q; i++)
      query[i] = new float[D];
    /* Read Data */
    for (int i = 0; i < N; i++) {
      getline(fs_data, line);
      stringstream sts(line);
      for (int j = 0; j < D; j++)
        sts >> data[i][j];
    }
    /* Read Query */
     for (int i = 0; i < q; i++) {
      getline(fs_query, qline);
      stringstream sts(qline);
      for (int j = 0; j < D; j++)
        sts >> query[i][j];
    }
    fs_data.close();
    fs_query.close();
  }
};

class IOManager {
public:
  static fstream load_head(string path, int N, int D, int argc, ...) {
    va_list vl;
    int n, d, arg1;
    fstream gfile;
    string line;
    gfile.open(path, fstream::in);
    check_file(gfile, path);
    getline(gfile, line);
    stringstream sts(line);
    cout << "Loading from " << path <<"\n";
    sts >> n;
    sts >> d;
    assert(n == N && d == D);
    cout << "N : " << N << " / D : " << D << " / ";
    va_start(vl, argc);
    for (int i = 0; i < argc; i+=2) {
      cout << va_arg(vl, char*) << " : ";
      int arg2 = va_arg(vl, int);
      cout << arg2 << " / ";
      sts >> arg1;
      // assert(arg1 == arg2);
    }
    cout << "\n";
    return gfile;
  }

  static ofstream save_head(string path, int N, int D, int argc, ...) {
    va_list vl;
    ofstream ofile;
    cout << "Saving to " << path << "\n";
    // cout << "Saving (size = " << N << " / MaxDegree : " << M << ") to " << path << "\n";
    ofile.open(path, ofstream::out | ofstream::trunc);
    check_file(ofile, path);
    ofile << N << " " << D << " "; 
    va_start(vl, argc);
    for (int i = 0; i < argc; i++)
      ofile << va_arg(vl, int) << " ";
    ofile << "\n";
    return ofile;
  }
};

string replaceAll(const string& str, const string& pattern, const string& replace) {
  string result = str;
  string::size_type pos = 0;
  string::size_type offset = 0;
  while ((pos = result.find(pattern, offset)) != string::npos) {
    result.replace(result.begin() + pos, result.begin() + pos + pattern.size(), replace);
    offset = pos + replace.size();
  }
  return result;
}


