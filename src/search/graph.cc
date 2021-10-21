using phmap::flat_hash_map;
using phmap::flat_hash_set;
#include <cstdio>

// #define SPARSIFY
// #define LOAD_AND_SPARSIFY

class Graph {
public:
  int N, D, M;
  vector<vector<int>>* graph;
  Graph() { 
  }
  ~Graph() {
  }
};
class GraphConstruct {
public:
  static int N;
  static int D;
  static float** data;
  static vector<float> norms;
  static float THRESHOLD;
  static int S_SIZE;

  GraphConstruct() {}
  ~GraphConstruct() {}

  static void init(int n, int d, float** data_in, int s_size, float threshold) {
    N = n;
    D = d;
    data = data_in;
    compute_norms();
    THRESHOLD = threshold;
    S_SIZE = s_size;
  }
  static void compute_norms() {
    for(int i = 0; i < N; i++)  {
      float v = sqrt(squared_sum(data[i], D));
      if(v == 0.0f)
        v = 1.0f;
      norms.push_back(v);
    }
  }
  static void bidirectionalize_light(flat_hash_map<int, vector<int>> &igraph) {
    flat_hash_map<int, flat_hash_set<int>> graph_set;
    flat_hash_map<int, operand> indeg;

    for(auto it : igraph) {
     flat_hash_set<int> s;
     graph_set.insert(make_pair(it.first, s));
     indeg.insert(make_pair(it.first, make_pair(it.first, 0)));
    }
    // #pragma omp parallel for
    for(auto it : igraph) {
     flat_hash_set<int> &s = graph_set[it.first];
     for(int j = 0; j < igraph[it.first].size(); j++) {
       int cur = igraph[it.first][j];
       indeg[cur].second++;
       s.insert(cur);
     }
    }
    vector<operand> sorted_indeg;
    for (auto it : indeg) 
      sorted_indeg.push_back(it.second);
    sort(sorted_indeg.begin(), sorted_indeg.end(), compare_descending);

    long long int total_added = 0;
    int maxAddID = -1;
    int maxAdd = -1;
    int cnt = 0;
    int threshold = sorted_indeg[(int)(3*igraph.size())/4].second;
    for(int i = sorted_indeg.size()-1; i >= 0; i--) {
     if(sorted_indeg[i].second < threshold) {
      cnt++;
       int cur = sorted_indeg[i].first;
       int added = 0;
       for(int j = 0; j < igraph[cur].size(); j++) {
         int target = igraph[cur][j];
         if(graph_set[target].count(cur) == 0) {
           graph_set[target].insert(cur);
           igraph[target].push_back(cur);
           added++;
         }
         if(added >= threshold - sorted_indeg[i].second)
           break;
       }
       if (added > maxAdd) {
        maxAdd = added;
        maxAddID = cur;
       }
       total_added += added;
      }
    }
  }

  static void bidirectionalize_light(vector<vector<int>> &graph) {
    vector<flat_hash_set<int>> graph_set;
    vector<operand> indeg;

    for(int i = 0; i < graph.size(); i++) {
     flat_hash_set<int> s;
     graph_set.push_back(s);
     indeg.push_back(make_pair(i, 0));
    }
    // #pragma omp parallel for
    for(int i = 0; i < graph.size(); i++) {
     flat_hash_set<int> &s = graph_set[i];
     for(int j = 0; j < graph[i].size(); j++) {
       int cur = graph[i][j];
       indeg[cur].second++;
       s.insert(cur);
     }
    }
    sort(indeg.begin(), indeg.end(), compare_descending);
    long long int total_added = 0;
    int maxAddID = -1;
    int maxAdd = -1;
    int cnt = 0;
    int threshold = indeg[(int)(3*N)/4].second;
    for(int i = indeg.size()-1; i >= 0; i--) {
     if(indeg[i].second < threshold) {
      cnt++;
       int cur = indeg[i].first;
       int added = 0;
       for(int j = 0; j < graph[cur].size(); j++) {
         int target = graph[cur][j];
         if(graph_set[target].count(cur) == 0) {
           graph_set[target].insert(cur);
           graph[target].push_back(cur);
           added++;
         }
         if(added >= threshold - indeg[i].second)
           break;
       }
       if (added > maxAdd) {
        maxAdd = added;
        maxAddID = cur;
       }
       total_added += added;
     }
    }
    cout << "Verification start...\n";
    vector<operand> new_indeg;
    int totalOut = 0;
    for(int i = 0; i < graph.size(); i++) {
     new_indeg.push_back(make_pair(i, 0));
    }

    for(int i = 0; i < graph.size(); i++) {
      totalOut += graph[i].size();
     for(int j = 0; j < graph[i].size(); j++) {
       int cur = graph[i][j];
       new_indeg[cur].second++;
     }
    }
    int count = 0;
    for (int i = 0; i < graph.size(); i++) {
      if (new_indeg[i].second < threshold) {
        count++;
        cout << i << "'s indegree is not enough : " << new_indeg[i].second << ", threshold : " << threshold << "\n";
      }
    }
    cout << "Max Added ID : " << maxAddID << " maxAdded : " << maxAdd << "\n";
    cout << "Total " << count << " nodes don't have enough edges\n";
    cout << "Total " << total_added << " edges are added to " << cnt << " nodes, average added edge : " << (float)total_added / cnt << "\n";
    cout << "Average out degree : " << (float)totalOut / graph.size() << "\n";
    cout << "Bidirectionalize light version finished\n";
  }

  static vector<int> sparsify(vector<int> &neighbor, float* p, bool parallelize=false) {
      vector<operand> res;
      res.reserve(2048);
      int insert = 1;
      for(int i=0; i<insert; i++)
        res.push_back(make_pair(neighbor[i], inner_product(p, data[neighbor[i]], D)));
      for(int i=insert; i<neighbor.size(); i++) {
        int q = neighbor[i];  
        float pq = acos(inner_product(p, data[q], D));
        //float pq = 1 - inner_product(data[p], data[q], D);
        bool flag = false;
        bool dopar = parallelize && (res.size() >= 16);
        #pragma omp parallel for if(dopar) num_threads(OMD)
        for(int j=0; j<res.size(); j++) {
          if(flag) 
            continue;
          int r = res[j].first;
          float pr = acos(res[j].second); 
          float rq = acos(inner_product(data[r], data[q], D));
          float cosres = (pr * pr + pq * pq - rq * rq)/(2*pr * pq);
          float threshold = THRESHOLD;
          if(cosres > threshold) {
            #pragma omp critical
            {
              flag = true;
            }
          }
        }
        if(!flag) {
          float sim = inner_product(p, data[q], D);
          res.push_back(make_pair(neighbor[i], sim));
        }
        if(res.size() >= 1024)
          break;
      }

      vector<int> resret;
      resret.reserve(res.size());
      for(int i=0; i<res.size(); i++) {
        resret.push_back(res[i].first);
      }
      return resret;
    }
  static vector<int> sparsify(vector<int> &neighbor, int p, bool parallelize=false) {
    vector<operand> res;
    res.reserve(2048);
    int insert = 1;
    for(int i=0; i<insert; i++)
      res.push_back(make_pair(neighbor[i], inner_product(data[p], data[neighbor[i]], D)));
    for(int i=insert; i<neighbor.size(); i++) {
      int q = neighbor[i];  
      float pq = acos(inner_product(data[p], data[q], D));
      //float pq = 1 - inner_product(data[p], data[q], D);
      bool flag = false;
      bool dopar = parallelize && (res.size() >= 16);
      #pragma omp parallel for if(dopar)  num_threads(OMD)
      for(int j=0; j<res.size(); j++) {
        if(flag) 
          continue;
        int r = res[j].first;
        float pr = acos(res[j].second); 
        float rq = acos(inner_product(data[r], data[q], D));
        float cosres = (pr * pr + pq * pq - rq * rq)/(2*pr * pq);
        float threshold = THRESHOLD;
        if(cosres > threshold) {
          #pragma omp critical
          {
            flag = true;
          }
        }
      }
      if(!flag) {
        float sim = inner_product(data[p], data[q], D);
        res.push_back(make_pair(neighbor[i], sim));
      }
      if(res.size() >= S_SIZE)
        break;
    }

    vector<int> resret;
    resret.reserve(res.size());
    for(int i=0; i<res.size(); i++) {
      resret.push_back(res[i].first);
    }
    return resret;
  }
  static vector<int> sparsify_seed(vector<int> &neighbor, int p, int deg, int minneighbor) {
    vector<operand> res;
    float th = 0.5;
    if(neighbor.size() <= minneighbor) {
      return neighbor;
    }
    int lastcut = neighbor.size();
    bool run = true;
    bool touch = false;
    int maxiter = 0;
    while(run) {
      res.clear();
      int insert = 1;
      maxiter++;
      for(int i=0; i<insert; i++)
        res.push_back(make_pair(neighbor[i], inner_product(data[p], data[neighbor[i]], D)));
      int rdeg = min(deg, (int)neighbor.size());
      for(int i=insert; i<rdeg; i++) {
        int q = neighbor[i];  
        float pq = acos(inner_product(data[p], data[q], D));
        //float pq = 1 - inner_product(data[p], data[q], D);
        bool flag = false;
        //#pragma omp parallel for if(dopar)  num_threads(OMD)
        for(int j=0; j<res.size(); j++) {
          int r = res[j].first;
          float pr = acos(res[j].second); 
          float rq = acos(inner_product(data[r], data[q], D));
          float cosres = (pr * pr + pq * pq - rq * rq)/(2*pr * pq);
          float threshold = th;
          if(cosres > threshold) {
            flag = true;
            break;
          }
        }
        if(!flag) {
          float sim = inner_product(data[p], data[q], D);
          res.push_back(make_pair(neighbor[i], sim));
        }
        if(res.size() > 80) {
          if(touch)
            run = false;
          th = th - 0.01;
          touch = false;
          break;
        }
      }
      if(res.size() < 64) {
        touch = true;
        th = th + 0.01;
      }
      if(maxiter == 10)
        break;
      if(res.size() < 80 && res.size() >= 64)
        break;
    }
    //cout << "Res : " <<res.size() <<"\n";
    vector<int> resret;
    resret.reserve(res.size());
    for(int i=0; i<min((int)res.size(),64); i++) {
      resret.push_back(res[i].first);
    }
    return resret;
  }
  static vector<operand> sparsify(vector<operand> &neighbor, int p, int deg) { 
    vector<operand> res;
    res.reserve(64);
    int insert = 1;
    for(int i = 0; i<insert; i++)
      res.push_back(make_pair(neighbor[0].first, neighbor[0].second));
    for(int i= insert; i<deg; i++) {
      int q = neighbor[i].first;
      float pq = acos(neighbor[i].second);
      bool flag = false;
      bool dopar = (res.size() >= 16);
      #pragma omp parallel for if(dopar)  num_threads(OMD)
      for(int j=0; j<res.size(); j++) {
        if(flag) 
          continue;
        int r = res[j].first;
        float pr = acos(res[j].second); 
        float rq = acos(inner_product(data[r], data[q], D));
        float cosres = (pr * pr + pq * pq - rq * rq)/(2*pr * pq);
        float threshold = THRESHOLD;
        if(cosres > threshold) {
          #pragma omp critical
          {
            flag = true;
          }
        }
      }
      if(!flag) {
        res.push_back(neighbor[i]);
      }
    }
    return res;
  }
  static vector<vector<int>> construct(int M, bool includeSelfConnect, float (*compute_similarity)(float*, float*, int, float, float)) {
    vector<vector<int>> graph;
    graph.reserve(N);
    for (int i = 0; i < N; i++)
      graph.push_back(vector<int>());
    bool print = N >= 100000 ? true : false;
    int progress = 0;
    cout << "Graph Generation Start.\n";
    long long int tote = 0;
    long long int tott = 0;
    float similarity = 0.0f;
    #pragma omp parallel for schedule(dynamic) num_threads(OMD)
    for (int i = 0; i < N; i++) {
      vector<operand> q;
      if (print && (i % (N / 100) == 0)) {
        #pragma omp critical
        {
          progress++;
          if (progress % 5 == 0)
            cout << "Graph Generation Progress [ " << progress << "% ]\n";
        } 
      }
      for (int j = 0; j < N; j++) { 
        if ((!includeSelfConnect && i !=j) || includeSelfConnect) {
          similarity = compute_similarity(data[i], data[j], D, norms[i], norms[j]);
          q.push_back(make_pair(j, similarity));
        }
      }
      size_t deg = (size_t)M > q.size() ? q.size() : M;
      partial_sort(q.begin(), q.begin() + deg, q.end(), compare_descending);
      graph[i].reserve(deg);
      for (int j = 0; j < deg; j++)
        graph[i].push_back(q[j].first);
      #pragma omp critical
      {
        tote += graph[i].size();
        tott++;
      }
    }
    graph.shrink_to_fit();
    cout << "Average Edge : " << (float)tote / tott << "\n";
    cout << "Graph Generation Complete.\n";
    return graph;
  }

  static vector<flat_hash_map<int, vector<int>>> construct_partial_graph_in_graph(int part_num, int per_part, int M, vector<vector<int>>& graph_in, float(*compute_similarity)(float*, float*, int, float, float)) {
    vector<flat_hash_map<int, vector<int>>> graph;
    graph.reserve(per_part);
    for (int i = 0; i < per_part; i++)
      graph.push_back(flat_hash_map<int, vector<int>>());
    assert(graph_in.size() == N);
    bool print = N >= 10000 ? true : false;
    int progress = 0;
    cout << "InGraph Generation Start.\n";
    long long int tote = 0;
    long long int tott = 0;
    #pragma omp parallel for schedule(dynamic) num_threads(OMD)
    for (int i = 0; i < per_part; i++) {
      int cte = 0;
      int ctt = 0;
      vector<int>& neighbors = graph_in[part_num * per_part + i];
      int inGraph_N = neighbors.size();
      float similarity = 0.0f;
      if (print && (i % (per_part / 100) == 0)) {
        #pragma omp critical
        {
          progress++;
          if (progress % 5 == 0)
            cout << "Graph Generation Progress [ " << progress << "% ]\n";
        } 
      }
      for (int j = 0; j < inGraph_N; j++) {
        vector<operand> q;
        q.reserve(inGraph_N);
        for (int k = 0; k < inGraph_N; k++) { 
          if (neighbors[j] != neighbors[k]) {
            similarity = compute_similarity(data[neighbors[j]], data[neighbors[k]], D, norms[neighbors[j]], norms[neighbors[k]]);
            q.push_back(make_pair(neighbors[k], similarity));
          }
        }
        // if (q.size() < M) {
        //   #pragma omp critical
        //   {
        //     cout << "Warning, q.size (" << q.size() << ") is not enough, inGraph_N : " << inGraph_N << "\n";
        //   }
        // }
        size_t deg = (size_t)M > q.size() ? q.size() : M;
        partial_sort(q.begin(), q.begin() + deg, q.end(), compare_descending);
        size_t id = 0;
        graph[i].insert(make_pair(neighbors[j], vector<int>()));
        graph[i].reserve(deg);
        while (id < q.size() && id < deg) {
          int cur = q[id].first;
          graph[i][neighbors[j]].push_back(cur);
          id++;
        }
        ctt++;
        cte += graph[i][neighbors[j]].size();
      }
      #pragma omp critical
      {
        tote += cte;
        tott += ctt;
      }
    }
    cout << "Average Edge : " << (float)tote / tott << "\n";
    graph.shrink_to_fit();
    cout << "InGraph Generation Complete.\n";
    return graph;
  }

  

  static hnswlib::HierarchicalNSW<float>* construct_hnsw() {
    cout << "HNSW Graph Generation Start.\n";
    hnswlib::L2Space l2space(D); //L2Space -> Actually it's Innerproduct
    hnswlib::HierarchicalNSW<float> *hnswg;
    hnswg = new hnswlib::HierarchicalNSW<float>(&l2space, N, 32, 1024);
    for (int i = 0; i<1; i++)
      hnswg->addPoint((void *)(data[i]), (size_t)i);
    #pragma omp parallel for
    for (int i = 1; i<N; i++)
      hnswg->addPoint((void *)(data[i]), (size_t)i);
    cout << "HNSW Graph Generation Complete.\n";
    return hnswg;
  }

  static void save_hnsw(hnswlib::HierarchicalNSW<float>* hnswg, string path) {
    hnswg->SaveIndex(path);
  }

  static hnswlib::HierarchicalNSW<float>* load_hnsw(string path) {
    hnswlib::L2Space l2space(D);
    hnswlib::HierarchicalNSW<float>* hnswg;
    hnswg = new hnswlib::HierarchicalNSW<float>(&l2space, path, false);
    #ifdef CACHE
    system("sync; echo 1 > /proc/sys/vm/drop_caches");
    cout << "Cache cleared.\n";
    #endif
    return hnswg;
  }

  static void save_graph (int M, vector<vector<int>>& graph, string path) {
    ofstream ofile = IOManager::save_head(path, graph.size(), D, 1, M);
    for(int i = 0; i < graph.size(); i++) {
      ofile << i << " " << graph[i].size() << "\n";
      for (auto itit = graph[i].begin(); itit != graph[i].end(); ++itit)
        ofile << *itit << " ";
      ofile << "\n";
    }
    ofile.close();
    cout << "Finished saving the graph.\n";
  }

  static void save_partial_in_graph (int part_num, int per_part, vector<flat_hash_map<int, vector<int>>>& graph, string path, int ingraph_num_node, int in_graph_M = -1) {
    ofstream ofile = IOManager::save_head(path, graph.size(), D, 1, in_graph_M);
    // assert(graph.size() == per_part);
    for(int i = 0; i < per_part; i++) {
      ofile << part_num * per_part + i << " " << graph[i].size() << "\n";
      if (graph[i].size() > ingraph_num_node)
        cout << "asdfadfsasdfd : " << graph[i].size() << "\n";
      assert(graph[i].size() <= ingraph_num_node);
      int count = 0;
      for (auto itit = graph[i].begin(); itit != graph[i].end(); ++itit) {
        count++;
        ofile << (*itit).first << " " << (*itit).second.size() << "\n";
        // if((*itit).second.size() < in_graph_M)
        //   cout << "Node" << part_num * per_part + i << "'s " << (*itit).first << "'s neighbor is not enough, neighbor size : " << (*itit).second.size() << " / ingraph_M : " << in_graph_M << "\n"; 
        for (auto ititit = itit->second.begin(); ititit != itit->second.end(); ++ititit)
          ofile << *ititit << " ";
        ofile << "\n";
      }
      for (int j = count; j < ingraph_num_node; j++)
        ofile << "\n\n";
    }
    ofile.close();
    cout << "Finished saving the graph.\n";
  }

  static vector<vector<int>> load_graph(string graph_path, int M = -1) {
    fstream gfile = IOManager::load_head(graph_path, N, D, 0);
    long long int tote = 0;    
    long long int tott = 0;
    int partition = 8;
    vector<vector<int>> graph;
    graph.reserve(N);
    for (int i = 0; i < N; i++)
      graph.push_back(vector<int>());
    int part = 0;
    int lastpart = 0;
    if(N % partition == 0) {
      part = N/partition;
      lastpart = part;
    }
    else {
      part = N/partition + 1;
      lastpart = N - part * (partition-1);
    }
    for(int p = 0; p < partition; p++) {
      cout << "Current partition : " << p << " out of " << partition << "\n";
      if (p == partition - 1)
        part = lastpart;
      vector<string> lines(part * 2);
      for(int i = 0; i < 2 * part; i++)
        getline(gfile, lines[i]);
      cout << "Reading from " << graph_path << " is done, data processing starts...\n";
      #pragma omp parallel for num_threads(OMD)
      /* Read Data */
      for (int i = 0; i < part; i++) {
        int nodeIdx = -1;
        int numNeighbor = -1;
        sscanf(lines[i * 2].c_str(), "%d %d", &nodeIdx, &numNeighbor);
        vector<int> neighbors;
        char *lstr = strdup(lines[i * 2 + 1].c_str());
        char* next_ptr;
        char* elem = strtok_r(lstr," ", &next_ptr);
        if (M != -1) {
          neighbors.reserve(M);
          for (int j = 0; j < M; j++) {
            neighbors.push_back(stoi(elem));
            elem = strtok_r(NULL, " ", &next_ptr);
      if (elem == NULL)
        break;
          }
        }
        else {
          neighbors.reserve(1024);
          while (elem != NULL) {
            neighbors.push_back(stoi(elem));
            elem = strtok_r(NULL, " ", &next_ptr);
          }
        }
        neighbors.shrink_to_fit();
        graph[nodeIdx] = neighbors;
        #pragma omp critical
        {
          tott++;
          tote += graph[nodeIdx].size();
        }
      }
    }
    cout << "Average Edge : " << (float)tote / tott << "\n";
    gfile.close();  
    #ifdef CACHE
    system("sync; echo 1 > /proc/sys/vm/drop_caches");
    cout << "Cache cleared.\n";
    #endif
    cout <<"Finished loading the graph.\n";
    graph.shrink_to_fit();
    int counter[21] = {0,};  // 100단위 0~2100
    cout << "Graph Analysis\n";
    for (int i = 0; i < graph.size(); i++) {
      counter[graph[i].size()/100]++;
    }
    for (int i = 0; i < 21; i++)
      cout << i*100 << " ~ " << (i+1)*100 << " : " << counter[i] << "\n";
    return graph;
  }

  static void load_partial_in_graph_for_bidirection(int part_num, int per_part, int graph_M, int in_graph_M, string graph_path, vector<flat_hash_map<int, vector<int>>>& graph, bool readM = false) {
    fstream gfile = IOManager::load_head(graph_path, per_part, D, 0);
    long long int numLine = (long long int)per_part + (long long int)per_part * 2 * graph_M;
    #ifdef CACHE
    system("sync; echo 1 > /proc/sys/vm/drop_caches");
    cout << "Cache cleared.\n";
    #endif
    cout << "Reading from " << graph_path << " is done, data processing starts...\n";
    for (int i = 0; i < per_part; i++)
      graph.push_back(flat_hash_map<int, vector<int>>());
    /* Read Data */
    long long int tote = 0;
    long long int tott = 0;
    // #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < per_part; i++) {
      int nodeIdx = -1;
      int numNeighbor = -1;
      long long int lnsize = 0;
      string curline;
      getline(gfile, curline);
      // cout << "curline : " << curline << "\n";
      sscanf(curline.c_str(), "%d %d", &nodeIdx, &numNeighbor);
      // sscanf(lines[(long long int)i * (2*graph_M +1)].c_str(), "%d %d", &nodeIdx, &numNeighbor);
      assert(nodeIdx != -1);
      assert(numNeighbor >= 0);
      vector<string> lines(graph_M * 2);
      for (int i = 0; i < graph_M * 2; i++)
        getline(gfile, lines[i]);
      graph[nodeIdx - part_num * per_part].reserve(numNeighbor);
      #pragma omp parallel for num_threads(OMD)
      for (int j = 0; j < numNeighbor; j++) {
        int nodeNodeIdx = -1;
        int numNeighborNeighbor = -1;
        sscanf(lines[(long long int)j * 2].c_str(), "%d %d", &nodeNodeIdx, &numNeighborNeighbor);
        vector<int> neighbors;
        if(in_graph_M != -1)
          neighbors.reserve(in_graph_M);
        else
          neighbors.reserve(128);
        char *lstr = strdup(lines[(long long int)j * 2 + 1].c_str());
        char* next_ptr;
        char* elem = strtok_r(lstr," ", &next_ptr);
        while(elem != NULL){
          neighbors.push_back(stoi(elem));
          elem = strtok_r(NULL, " ", &next_ptr);
          if (readM && neighbors.size() == in_graph_M)
            break;
        }
        #pragma omp critical
        {
          lnsize += neighbors.size();
          graph[nodeIdx - part_num * per_part].insert(make_pair(nodeNodeIdx, neighbors));
          graph[nodeIdx - part_num * per_part][nodeNodeIdx].shrink_to_fit();
        }
      }
      graph[nodeIdx - part_num * per_part].rehash(graph[nodeIdx - part_num * per_part].size());
      #pragma omp critical
      {
        tote += lnsize;
        tott += numNeighbor;
      }
    }
    gfile.close();
    cout << "Average Edges : " << (float) tote / tott << "\n";
    cout <<"Finished loading the graph.\n"; 
  } 
  static void load_partial_in_graph_for_bidirection_and_cut(int part_num, int per_part, int graph_M, int load_front_M, string graph_path, vector<flat_hash_map<int, vector<int>>>& graph, vector<vector<int>>& sfgbl) {
    fstream gfile = IOManager::load_head(graph_path, per_part, D, 0);
    long long int numLine = (long long int)per_part + (long long int)per_part * 2 * graph_M;
    cout << "size : " << numLine << "\n";

    #ifdef CACHE
    system("sync; echo 1 > /proc/sys/vm/drop_caches");
    cout << "Cache cleared.\n";
    #endif
    cout << "Reading from " << graph_path << " is done, data processing starts...\n";
    for (int i = 0; i < per_part; i++)
      graph.push_back(flat_hash_map<int, vector<int>>());
    /* Read Data */
    long long int tote = 0;
    long long int tott = 0;
    // #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < per_part; i++) {
      int nodeIdx = -1;
      int numNeighbor = -1;
      long long int lnsize = 0;
      string curline;
      getline(gfile, curline);
      // cout << "curline : " << curline << "\n";
      sscanf(curline.c_str(), "%d %d", &nodeIdx, &numNeighbor);
      // sscanf(lines[(long long int)i * (2*graph_M +1)].c_str(), "%d %d", &nodeIdx, &numNeighbor);
      assert(nodeIdx != -1);
      assert(numNeighbor >= 0);
      vector<string> lines(graph_M * 2);
      for (int i = 0; i < graph_M * 2; i++)
        getline(gfile, lines[i]);
      flat_hash_set<int> front_M_nodes;
      int deg = min(load_front_M, (int)sfgbl[nodeIdx].size());
      front_M_nodes.reserve(deg);
      for (int j = 0; j < deg; j++) {
        front_M_nodes.insert(sfgbl[nodeIdx][j]);
      }
      graph[nodeIdx - part_num * per_part].reserve(numNeighbor);
      #pragma omp parallel for num_threads(OMD)
      for (int j = 0; j < numNeighbor; j++) {
        int nodeNodeIdx = -1;
        int numNeighborNeighbor = -1;
        sscanf(lines[(long long int)j * 2].c_str(), "%d %d", &nodeNodeIdx, &numNeighborNeighbor);
        vector<int> neighbors;
        neighbors.reserve(128);
        char *lstr = strdup(lines[(long long int)j * 2 + 1].c_str());
        char* next_ptr;
        char* elem = strtok_r(lstr," ", &next_ptr);
        while(elem != NULL){
          int node = stoi(elem);
          if (front_M_nodes.find(node) == front_M_nodes.end())
            neighbors.push_back(node);
          elem = strtok_r(NULL, " ", &next_ptr);
        }
        #pragma omp critical
        {
          lnsize += neighbors.size();
          graph[nodeIdx - part_num * per_part].insert(make_pair(nodeNodeIdx, neighbors));
          graph[nodeIdx - part_num * per_part][nodeNodeIdx].shrink_to_fit();
        }
      }
      graph[nodeIdx - part_num * per_part].rehash(graph[nodeIdx - part_num * per_part].size());
      #pragma omp critical
      {
        tote += lnsize;
        tott += numNeighbor;
      }
    }
    gfile.close();
    cout << "Average Edges : " << (float) tote / tott << "\n";
    cout <<"Finished loading the graph.\n"; 
  } 
  static void load_partial_in_graph(int part_num, int per_part, int graph_M, int in_graph_M, string graph_path, vector<flat_hash_map<int, vector<int>>>& graph, bool readM = false) {
    fstream gfile = IOManager::load_head(graph_path, per_part, D, 0);
    long long int numLine = (long long int)per_part + (long long int)per_part * 2 * graph_M;
    cout << "size : " << numLine << "\n";

    vector<string> lines(numLine);
    for(long long int i = 0; i < numLine; i++)
      getline(gfile, lines[i]);
    gfile.close();
    #ifdef CACHE
    system("sync; echo 1 > /proc/sys/vm/drop_caches");
    cout << "Cache cleared.\n";
    #endif
    cout << "Reading from " << graph_path << " is done, data processing starts...\n";
    for (int i = 0; i < per_part; i++)
      graph.push_back(flat_hash_map<int, vector<int>>());
    /* Read Data */
    long long int tote = 0;
    long long int tott = 0;
    #pragma omp parallel for schedule(dynamic) num_threads(OMD)
    for (int i = 0; i < per_part; i++) {
      int nodeIdx = -1;
      int numNeighbor = -1;
      long long int lnsize = 0;
      sscanf(lines[(long long int)i * (2*graph_M +1)].c_str(), "%d %d", &nodeIdx, &numNeighbor);
      assert(nodeIdx != -1);
      assert(numNeighbor >= 0);
      graph[nodeIdx].reserve(numNeighbor);
      for (int j = 0; j < numNeighbor; j++) {
        int nodeNodeIdx = -1;
        int numNeighborNeighbor = -1;
        sscanf(lines[(long long int)i * (2*graph_M +1) + 2 * j +1].c_str(), "%d %d", &nodeNodeIdx, &numNeighborNeighbor);
        vector<int> neighbors;
        if(in_graph_M != -1)
          neighbors.reserve(in_graph_M);
        else
          neighbors.reserve(128);
        char *lstr = strdup(lines[(long long int)i * (2 * graph_M + 1) + 2 * j + 2].c_str());
        char* next_ptr;
        char* elem = strtok_r(lstr," ", &next_ptr);
        while(elem != NULL){
          neighbors.push_back(stoi(elem));
          elem = strtok_r(NULL, " ", &next_ptr);
          if (readM && neighbors.size() == in_graph_M)
            break;
        }
        lnsize += neighbors.size();
        graph[nodeIdx].insert(make_pair(nodeNodeIdx, neighbors));
        graph[nodeIdx][nodeNodeIdx].shrink_to_fit();
      }
      graph[nodeIdx].rehash(graph[nodeIdx].size());
      #pragma omp critical
      {
        tote += lnsize;
        tott += numNeighbor;
      }
    }
    cout << "Average Edges : " << (float) tote / tott << "\n";
    cout <<"Finished loading the graph.\n"; 
  } 

};
int GraphConstruct::N = -1;
int GraphConstruct::D = -1;
int GraphConstruct::S_SIZE = -1;
float** GraphConstruct::data = NULL;
float GraphConstruct::THRESHOLD = 0.0f;
vector<float> GraphConstruct::norms = vector<float>();
