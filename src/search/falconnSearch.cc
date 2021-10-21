#include <falconn/lsh_nn_table.h>

using falconn::construct_table;
using falconn::compute_number_of_hash_functions;
using falconn::DenseVector;
using falconn::DistanceFunction;
using falconn::LSHConstructionParameters;
using falconn::LSHFamily;
using falconn::LSHNearestNeighborTable;
using falconn::LSHNearestNeighborQuery;
using falconn::QueryStatistics;
using falconn::StorageHashTable;
using falconn::get_default_parameters;

typedef DenseVector<float> Point;

class FalconnSearch : public BaseSearch {
public:
  int M;
  const int NUM_QUERIES = 1000;
  const int SEED = 4057218;
  const int NUM_HASH_TABLES = 50;
  const int NUM_HASH_BITS = 18;
  const int NUM_ROTATIONS = 1;
  const float want_recall = 0.95;
  vector<Point> dataset;
  vector<Point> sample_dataset;
  vector<vector<int>> answers;
  Point center;
  unique_ptr<LSHNearestNeighborTable<Point>> table;
  unique_ptr<LSHNearestNeighborQuery<Point>> query_object;
  int num_probes;

  FalconnSearch(string name, int k, int nSearch) : BaseSearch(name, k, nSearch) {
  }
  ~FalconnSearch() {
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
    dataset.clear();
    cout << "[FALCONN-preprocess] re-align to dense vector" << endl;
    for(int i = 0; i<N; i++) {
      Point p;
      p.resize(D);
      for (int j = 0; j<D; j++){
        p[j] = data[i][j];
      }
      p.normalize();
      dataset.push_back(p);
    }

    // Find the center 
    center = dataset[0];
    for (size_t i = 1; i < dataset.size(); ++i) {
      center += dataset[i];
    }
    center /= dataset.size();

    // Recentering Data
    for (auto &datapoint : dataset) {
      datapoint -= center;
    }

    // selecting NUM_QUERIES data points as queries
    cout << "[FALCONN-preprocess] sample datapoints" << endl;
    gen_queries(&dataset, &sample_dataset);

    // running the linear scan
    cout << "[FALCONN-prerpocess] Find groundtruth for samples" << endl;
    gen_answers(dataset, sample_dataset, &answers);

    auto t1 = Clock::now();
    
    // Parameter Setting for Falconn 
    LSHConstructionParameters params;
    params.dimension = dataset[0].size();
    params.lsh_family = LSHFamily::CrossPolytope;
    params.l = NUM_HASH_TABLES;
    params.distance_function = DistanceFunction::NegativeInnerProduct;
    compute_number_of_hash_functions<Point>(NUM_HASH_BITS, &params);
    params.num_rotations = NUM_ROTATIONS;
    params.num_setup_threads = 0; //For multithread Falconn
    params.storage_hash_table = StorageHashTable::BitPackedFlatHashTable;

    cout << "[FALCONN-preprocess/main] Building the index based on the cross-polytope LSH" << endl;
    table = construct_table<Point>(dataset, params);

    // finding the number of probes via the binary search
    cout << "[FALCONN-prerpocess/main] finding the appropriate number of probes" << endl;
    num_probes = find_num_probes(&*table, sample_dataset, answers, params.l);
    query_object = table->construct_query_object(num_probes, N/2);

    auto t2 = Clock::now();
    timers[0] += compute_clock(t2, t1);
    ready = true;
  }
  float inline getSimilarity(int id, float *current_query) {
    return inner_product(data[id], current_query, D);
  }

  vector<operand> search(float* current_query){
    Point query;
    query.resize(D);

    for (int i = 0; i<D; i++) 
      query[i] = current_query[i];
    query.normalize();
    query -= center;
    
    assert(ready == true);
    auto t1 = Clock::now();
    vector<int32_t> result;
    query_object->find_k_nearest_neighbors(query, (int_fast64_t)K, &result);
    for (int r : result) {
      top_k_candidates.push_back(make_pair(r, 0));
    }
    auto t2 = Clock::now();
    timers[1] += compute_clock(t2, t1);
    return top_k_candidates;
  }

  void print_time() {
    cout << "[FalconnSearch] " << (timers[1] ) /queryCount << "(ms)\n";
    cout << "Preprocessing : " << timers[0] << "(ms)\n";
    cout << "Search " << (timers[1]) / queryCount << "(ms)\n";
   }
  void print_information() {
    cout << "[FalconnSearch] \n";
  }
  
  /*
   * Chooses a random subset of the dataset to be the queries. The queries are
   * taken out of the dataset.
   */
  void gen_queries(vector<Point> *dataset, vector<Point> *queries) {
    mt19937_64 gen(SEED);
    queries->clear();
    for (int i = 0; i < NUM_QUERIES; ++i) {
      uniform_int_distribution<> u(0, dataset->size() - 1);
      int ind = u(gen);
      queries->push_back((*dataset)[ind]);
      (*dataset)[ind] = dataset->back();
      dataset->pop_back();
    }
  }

  /*
   * Generates answers for the queries using the (optimized) linear scan.
   */
  void gen_answers(const vector<Point> &dataset, const vector<Point> &queries,
                   vector<vector<int>> *answers) {
    int outer_counter = 0;
    for (const auto &query : queries) {
      float best = -10.0;
      int inner_counter = 0;
      vector<operand> op;
      answers->push_back(vector<int>());
      for (const auto &datapoint : dataset) {
        float score = query.dot(datapoint);
        op.push_back(make_pair(inner_counter, score));
        ++inner_counter;
      }
      partial_sort(op.begin(), op.begin()+K, op.end(), compare_descending);
      for (int i = 0; i<K; i++) {
        (*answers)[outer_counter].push_back(op[i].first);
      }
      ++outer_counter;
    }
  }

  /*
   * Computes the probability of success using a given number of probes.
   */
  double evaluate_num_probes(LSHNearestNeighborTable<Point> *table,
                             const vector<Point> &queries,
                             const vector<vector<int>> &answers, int num_probes) {
    unique_ptr<LSHNearestNeighborQuery<Point>> query_object =
        table->construct_query_object(num_probes);
    int outer_counter = 0;
    int num_matches = 0;
    vector<int32_t> candidates;
    for (const auto &query : queries) {
      query_object->get_candidates_with_duplicates(query, &candidates);
      for (int answer_k : answers[outer_counter]){
        for (auto x : candidates) {
          if (answer_k == x) {
            ++num_matches;
            break;
          }
        }
      }
      ++outer_counter;
    }
    return (num_matches + 0.0) / ((queries.size() + 0.0) * K);
  }

  int find_num_probes(LSHNearestNeighborTable<Point> *table,
                      const vector<Point> &queries, const vector<vector<int>> &answers,
                      int start_num_probes) {
    int num_probes = start_num_probes;
    for (;;) {
      cout << "[FALCONN-preprocess] trying " << num_probes << " probes" << endl;
      double precision = evaluate_num_probes(table, queries, answers, num_probes);
      if (precision >= want_recall) {
        break;
    }
      num_probes *= 2;
    }

    int r = num_probes;
    int l = r / 2;

    while (r - l > 1) {
      int num_probes = (l + r) / 2;
      cout << "[FALCONN-preprocess] trying " << num_probes << " probes" << endl;
      double precision = evaluate_num_probes(table, queries, answers, num_probes);
      if (precision >= want_recall) {
        r = num_probes;
      } else {
        l = num_probes;
      }
    }

    return r;
  }
}; 
