class GraphSearchWJY : public BaseSearch {
public:
  Graph* cg;
  int M;
  float** unitvec;
  vector<operand> selectvec;
  unordered_map<int, vector<int>> gx;
  vector<set<int>> cnt;
  vector<int> cnt_S;
  float cos_threshold = 0.5;
  int candidate_size = 1200;
  int S; 

  GraphSearchWJY(string name, int k, int nSearch) : BaseSearch(name, k, nSearch) {
    cg = new Graph();
  }
  ~GraphSearchWJY() {
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
/*
    #pragma omp parallel for
    for (int i = 0; i <N; i++) {
      gx.insert(make_pair(i, vector<int>()));
      for (int j = 0; j<N; j++) {
        if (i!=j)
          gx[i].push_back(j);
        //int rand_ = rand()%N;
        //if (find(gx[i].begin(), gx[i].end(), rand_) == gx[i].end()) {
        //  gx[i].push_back(rand_);
        //}
      }
    }
*/
    //int scccnt;
    //map<int,int> x = SCCAnalysis::tarjanSCC(cg->graph, scccnt);
    //cout << "SCC " <<  scccnt << endl;


    S=1200;
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
        //op.push_back(make_pair(j, getSimilarity(j, unitvec[i])));
      }
      /*O
      partial_sort(op.begin(), op.begin()+1024, op.end(), compare_descending);
      for (int j = 0; j<1024; j++) {
        gOx[i].push_back(op[j].first);
      }
      */
    }



/*
    //vector<vector<int>> cnt(S, vector<int>(N));
    cnt = vector<set<int>>(S, set<int>());
    cnt_S = vector<int>(S);

    int K_ = 1000;
    int progress = 0;

    #pragma omp parallel for
    for (int i = 0; i<N; i++) {
      vector<operand> unit_selector;
      for (int j = 0; j<S; j++) {
        unit_selector.push_back(make_pair(j, getSimilarity(i, unitvec[j])));
      }
      partial_sort(unit_selector.begin(), unit_selector.begin()+3, unit_selector.end(), compare_descending);

      if ((i % (N / 100) == 0)) {
       #pragma omp critical
        {
          progress++;
          if (progress % 5 == 0)
            cout << "Generation Progress [ " << progress << "% ]\n";
        } 
      }
 

      
      vector<int>& neighbor_i = cg->graph[i];
      

      #pragma omp critical
      {
        for (int s_=0; s_<1; s_++) {
          int s = unit_selector[s_].first;
          cnt_S[s]++;
          for (int j = 0; j<K_; j++) {
            int neighbor = neighbor_i[j];
            //cnt[s][neighbor]++;
            cnt[s].insert(neighbor);
          }
        }
      }
    }

    progress = 0;
    int extra = 2*N;
    extra = N;
    #pragma omp parallel for
    for (int i = 0; i<extra; i++) {
      float* unit;
      unit = new float[D];
      for (int j = 0; j<D; j++) {
        unit[j] = d1(generator);
      }
      float norm = sqrt(squared_sum(unit, D));
      for (int j = 0; j<D; j++) {
        unit[j] /= norm;
      }
 
      if ((i % (extra / 100) == 0)) {
        #pragma omp critical
        {
          progress++;
          if (progress % 5 == 0)
            cout << "Generation Progress [ " << progress << "% ]\n";
        } 
      }
 
      vector<operand> unit_selector;
      for (int j = 0; j<S; j++) {
        unit_selector.push_back(make_pair(j, inner_product(unitvec[j], unit, D)));
      }
      partial_sort(unit_selector.begin(), unit_selector.begin()+1, unit_selector.end(), compare_descending);

      
      vector<operand> true_top;
      for (int j = 0; j<N; j++) {
        true_top.push_back(make_pair(j, getSimilarity(j, unit)));
      }
      partial_sort(true_top.begin(), true_top.begin()+K_, true_top.end(), compare_descending);

      #pragma omp critical
      {
        
        for (int s_=0; s_<1; s_++) {
          int s = unit_selector[s_].first;
          cnt_S[s]++;
          for (int j = 0; j<K_; j++) {
            int neighbor = true_top[j].first;
            //cnt[s][neighbor]++;
            cnt[s].insert(neighbor);
          }
        }
      }

      delete[] unit;
    }

 
    for (int i = 0; i<S; i++) {
      cout << i << "th unitvec : " << cnt_S[i] << " / " << 11*N << " --- ";
      int zero = 0;
      cout << "Nonzero / Zero : " << cnt[i].size() << " / " << N-cnt[i].size() << endl;
      cout << endl;
    }
*/   


    auto t2 = Clock::now();
    timers[0] += compute_clock(t2, t1);
    ready = true;
  }
  float inline getSimilarity(int id, float *current_query) {
    return inner_product(data[id], current_query, D);
  }

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
    cout << " < " ;
    for (int i = 0; i < K ; i++) {
      for (size_t j = 0; j < check; j++) {
        if (true_top_k[i].first == top_k_candidates[j].first) {
          cout << "(" << i << "," << j << ") ";
          count++;
          break;
        }
      }
    }
    cout << " > " ;
    return (float)count / (float)K * 100.0;
  }

  int localcount = 0;
  float localrecall = 0;
  float neighborrecall = 0;
  int neighborcount = 0;

  float diversity(vector<int> R, float* current_query) {
    float** diff = new float*[R.size()];
    for (int i = 0; i < R.size(); i++) {
      diff[i] = new float[D];
      int p = R[i];
      for (int j = 0; j<D; j++){
        diff[i][j] = unitvec[p][j] - current_query[j];
      }
      float norm = sqrt(squared_sum(diff[i], D));
      for (int j = 0; j<D; j++) {
        diff[i][j] /= norm;
      } 
    }

    float* diffsum = new float[D];
    for (int i = 0; i<D; i++) {
      diffsum[i] = 0;
      for (int j = 0; j<R.size(); j++) {
        diffsum[i] += diff[j][i];
      }
    }

    float norm = sqrt(squared_sum(diffsum, D));

    return 1.0 - norm/(float)R.size();
  }


  vector<operand> search(float* current_query){
    assert(ready == true);
    vector<uint8_t> visited(N);
    vector<int> candidate;
    auto t1 = Clock::now();
    for (int i=0; i<S; i++) {
      if (inner_product(current_query, unitvec[i],D) > cos_threshold) {
        candidate.push_back(i);
      }
    }

    cout << candidate.size() << endl;

    for (int i : candidate) {
      for (int n : gx[i]) {
        visited[n]++;
      }
    }

    vector<intoperand> nums;
    auto t2 = Clock::now();
    for (int i = 0; i<N; i++) {
      uint8_t visits = visited[i];
      if (visits > 13)
        nums.push_back(make_pair(i, visits));
    }

    int sizes = nums.size() < 3000 ? nums.size() : 3000;
    partial_sort(nums.begin(), nums.begin()+sizes, nums.end(), int_compare_descending);

    auto t3 = Clock::now();
    for (int i = 0; i<sizes; i++) {
      top_k_candidates.push_back(make_pair(nums[i].first, getSimilarity(nums[i].first, current_query)));
    }
    
    partial_sort(top_k_candidates.begin(), top_k_candidates.begin()+K, top_k_candidates.end(), compare_descending);
     
    cout << (int)nums[0].second << " " << (int)nums[sizes-1].second << endl;


    /*
    vector<operand> true_top;
    for (int i = 0; i<N; i++) {
      true_top.push_back(make_pair(i, getSimilarity(i,current_query)));
    }
    
    partial_sort(true_top.begin(), true_top.begin()+K, true_top.end(), compare_descending);
    */



/*
    int e = 64;
    unordered_map<int, vector<int>> ingraph;

    vector<int>& neighbors = (*cg->graph)[rand()&N];
    
    for (int i : neighbors) {
      ingraph.insert(make_pair(i, vector<int>()));
      vector<operand> op;
      for (int j : neighbors) {
        if (i==j) continue;
        op.push_back(make_pair(j, getSimilarity(i, data[j])));
      }
      partial_sort(op.begin(), op.begin()+e, op.end(), compare_descending);
      for (int j=0; j<e; j++) {
        ingraph[i].push_back(op[j].first);
      }
    }


    for (int i = 0; i<100; i++) {
      int n = neighbors[i];
      vector<int>& neighbors_n = (*cg->graph)[n];

      vector<int> intersect;
      int cnt = 0;
      bool flag = false;
      for (int x : neighbors_n) {
        if (flag) {
          break;
        }
        for (int j : neighbors) {
          if (j==x) {
            intersect.push_back(j);
            cnt++;
            if (cnt == e) {
              cnt = 0;
              flag = true;
            }
            break;
          }
        }
      }
      
      cout << "intersect size : " << intersect.size() << " / " << flush;
      flag = false;
      for (int j=0; j<intersect.size(); j++) {
        if (ingraph[n][j] != intersect[j]) {
          flag = true;
          break;
        }
      }
      if (flag)
        cout << "wrong" << endl;
      else
        cout << "same" << endl;
    }

*/  
   /*
    for (int i=0; i<S; i++) {
      candidate.push_back(make_pair(i, inner_product(current_query, unitvec[i], D)));
    }
    
    partial_sort(candidate.begin(), candidate.begin()+1000, candidate.end(), compare_descending);

    int best;
    int max = 0;

    for (int i=0; i<200; i++) {
      int s = candidate[i].first;
      vector<int>& neighbor = gx[s];
      vector<operand> op;
      for (int n : neighbor) {
        op.push_back(make_pair(n, getSimilarity(n, current_query)));
      }
      partial_sort(op.begin(), op.begin()+K, op.end(), compare_descending);

      float recall = check_recall(true_top, op);
      if (max < recall) {
        best = s;
        max = recall;
      }

      if (i<200)
        cout << "S unitvec Top" << i +1 << " : " << candidate[i].first << " / " << candidate[i].second << " / " << recall << endl;
    }

    vector<int> h;
    float maxdiversity = 0;

    for (int i =0; i<200; i++) {
      if (candidate[i].first == best) continue;
      for (int j = 0; j<200; j++) {
        if (i==j || candidate[j].first==best) continue;
        vector<int> x{ best, candidate[i].first, candidate[j].first};
        float d = diversity(x, current_query);
        if (d > maxdiversity) {
          h = x;
          maxdiversity = d;
        }
      }
    }

    cout << "Max diversity : " << maxdiversity << " / Member : ";
    
    vector<operand> op;
    for (int i : h) {
      cout << i << " ";
      vector<int>& neighbor = gx[i];
      for (int n : neighbor) {
        if (visited[n]) continue;
        visited[n] = 1;
        op.push_back(make_pair(n, getSimilarity(n, current_query)));
      }
    }
    cout << endl;
    
    partial_sort(op.begin(), op.begin()+K, op.end(), compare_descending);
    float recall = check_recall(true_top, op);
  
    cout << "Final Size : " << op.size() << " / Recall : " << recall << endl;
    cout << endl;
*/
/*
    vector<operand> maxima;
    vector<operand> maximaneighbor;
    for (int i = 0; i<N; i++) {
      float sim = getSimilarity(i, current_query);
      bool localmaxima = true;
      for (int j = 0; j<1024; j++) {
        if (getSimilarity((*cg->graph)[i][j], current_query) > sim) {
          localmaxima = false;
          break;
        }
      }
      if (localmaxima) {
        maxima.push_back(make_pair(i, getSimilarity(i, current_query)));
        for (int j = 0; j<1024; j++) {
          if (visited[(*cg->graph)[i][j]]) continue;
          visited[(*cg->graph)[i][j]] = 1;
          maximaneighbor.push_back(make_pair((*cg->graph)[i][j], getSimilarity((*cg->graph)[i][j], current_query)));
        }
      }
    }
    localcount += maxima.size();
    
   
    vector<operand> true_top;
    for (int i = 0; i<N; i++) {
      true_top.push_back(make_pair(i, getSimilarity(i,current_query)));
    }
    
    partial_sort(true_top.begin(), true_top.begin()+K, true_top.end(), compare_descending);
    partial_sort(maximaneighbor.begin(), maximaneighbor.begin()+K, maximaneighbor.end(), compare_descending);
    sort(maxima.begin(), maxima.end(), compare_descending);

    neighborcount += maximaneighbor.size();
    localrecall += check_recall(true_top, maxima, true, false);
    neighborrecall += check_recall(true_top, maximaneighbor, true, false);
*/

/*
    for (int z = 0; z<6; z++) {
      vector<int>& top1_neighbor = cg->graph[true_top[z].first];

      float* virtual_top1 = new float[D];
      for (int i = 0; i<D; i++) virtual_top1[i] = 2*current_query[i]*true_top[z].second - data[true_top[z].first][i];

      cout << true_top[z].second << " vs " << inner_product(virtual_top1, current_query, D)  << endl;
      

      vector<operand> virtual_top;
      for (int i = 0; i<N; i++) {
        virtual_top.push_back(make_pair(i, inner_product(virtual_top1, data[i], D)));
      }
      sort(virtual_top.begin(), virtual_top.end(), compare_descending);

      for (int i = 0; i<1024; i++){
        int n1 = top1_neighbor[i];
        if (visited[n1] == 0) {
          candidate.push_back(make_pair(n1, getSimilarity(n1, current_query)));
          visited[n1] = 1;
        }
        int n2 = virtual_top[i].first;
        if (visited[n2] == 0) {
          candidate.push_back(make_pair(n2, getSimilarity(n2, current_query)));
          visited[n2] = 1;
        }
      }
    }
    cout << candidate.size() << endl;
    partial_sort(candidate.begin(), candidate.begin() + K, candidate.end(), compare_descending);
    cout << check_recall(true_top, candidate) << endl;
*/


/*
    vector<operand> unit_selector;
    for (int j = 0; j<S; j++) {
      unit_selector.push_back(make_pair(j, inner_product(current_query, unitvec[j], D)));
    }
    partial_sort(unit_selector.begin(), unit_selector.begin()+3, unit_selector.end(), compare_descending);

    int s = unit_selector[0].first;
    for (int x : cnt[s]) {
      visited[x] = 1;
      top_k_candidates.push_back(make_pair(x, getSimilarity(x, current_query)));
    }
  
    s = unit_selector[1].first;
    for (int x : cnt[s]) {
      if (visited[x]) continue;
      visited[x] = 1;
      top_k_candidates.push_back(make_pair(x, getSimilarity(x, current_query)));
    }
    
    s = unit_selector[2].first;
    for (int x : cnt[s]) {
      if (visited[x]) continue;
      visited[x] = 1;
      top_k_candidates.push_back(make_pair(x, getSimilarity(x, current_query)));
    }
    
    cout << top_k_candidates.size() << endl; 
    partial_sort(top_k_candidates.begin(), top_k_candidates.begin()+K, top_k_candidates.end(), compare_descending);
*/




/*
    int s=10;
    int rD = D/s;
    float* query_reduce = new float[rD];
    for (int i = 0; i<rD; i++) {
      query_reduce[i] = 0;
      for (int j = 0; j<s; j++) {
        query_reduce[i] += current_query[j+i*s];
      }
      query_reduce[i] /= s;
    }

    float** data_reduce = new float*[N];
    for(int i =0; i<N; i++) {
      data_reduce[i] = new float[rD];
      for (int j = 0; j<rD; j++) {
        data_reduce[i][j] = 0;
        for (int k = 0; k<s; k++) {
          data_reduce[i][j] += data[i][k+j*s];
        }
        data_reduce[i][j] /= s;
      }
    }

    vector<operand> candidate;
    for (int i = 0; i<N; i++) {
      float norm_data = sqrt(squared_sum(data_reduce[i], D));
      float norm_query = sqrt(squared_sum(query_reduce, D));
      candidate.push_back(make_pair(i, cosine_similarity(data_reduce[i], query_reduce, rD, norm_data, norm_query)));
    }

    partial_sort(candidate.begin(), candidate.begin() + K, candidate.end(), compare_descending);

    cout << check_recall(true_top, candidate, false, true) << endl;
*/



/*    vector<operand> true_top;
    vector<int> true_top_idx;
    auto t1 = Clock::now();
    for (int i = 0; i<N; i++) {
      true_top.push_back(make_pair(i, getSimilarity(i,current_query)));
    }
    //partial_sort(true_top.begin(), true_top.begin()+K, true_top.end(), compare_descending);
    sort(true_top.begin(), true_top.end(), compare_descending);
    for (int i = 0; i<1024; i++) {
      true_top_idx.push_back(true_top[i].first);
    }
  */
  /*
    float **x;
    x = simplex(D);
    float* query = new float[D];
    for (int i = 0; i<D; i++) {
      query[i] = 0;
      for (int j = 0; j<D; j++) {
        query[i] += x[j][i];
      }
    }
    float norm = sqrt(squared_sum(query, D));
    for (int i = 0; i<D; i++) {
      query[i] /= norm;
    }

    vector<operand> true_top;
    for (int i = 0; i<N; i++) {
      true_top.push_back(make_pair(i, getSimilarity(i,query)));
    }
    
    partial_sort(true_top.begin(), true_top.begin()+N, true_top.end(), compare_descending);
 
    for (int i = 0; i<D; i+=10) {
      vector<operand> op;
      for (int j = 0; j<N; j++) {
        op.push_back(make_pair(j, inner_product(data[j], x[i], D)));
      }
      partial_sort(op.begin(), op.begin()+100, op.end(), compare_descending);
      
      for (int j = 0; j<100; j++) {
        int idx = walk(x[i]);
        if (visited[idx]) continue;
        visited[idx] = 1;
        searched++;
        candidate.push_back(make_pair(idx, getSimilarity(idx, query)));
        

        //for (int k = 0; k<64; k++) {
        //  int n = neighborList[k];
        //  if (visited[n]) continue;
        //  visited[n] = 1;
        //  searched++;
        //  candidate.push_back(make_pair(n, getSimilarity(n, query)));    
        //}
      }
      
      
      //for (int j = 0; j<100; j++) {
      //  int idx = op[j].first;
      //  if (visited[idx]) continue;
      //  visited[idx] = 1;
      //  candidate.push_back(make_pair(idx, getSimilarity(idx, query)));
      //}
      
    }

    partial_sort(candidate.begin(), candidate.begin()+K, candidate.end(), compare_descending);
    cout << candidate.size() << " : " ;
    cout << check_recall(true_top, candidate) << endl;
*/

/*
    vector<operand> simxop;
    for (int i =0; i<D; i++) {
      simxop.push_back(make_pair(2*i, inner_product(x[i], current_query, D)));
      simxop.push_back(make_pair(2*i+1, -1.0*inner_product(x[i], current_query, D)));
    }
    partial_sort(simxop.begin(), simxop.begin() + 10, simxop.end(), compare_descending);
   
    for (int i = 0; i<10; i++) {
      cout << "Top" << i+1 << " : " << simxop[i].first << " / " << simxop[i].second << endl; 
    }
    cout << endl;

    vector<int> visited(N);
    vector<operand> candidate;

    for (int z=0; z<5; z++) {
      int top1 = (simxop[z].first)/2;
      float sign = (simxop[z].first)%2==1 ? -1.0 : 1.0;
      float **simx = new float*[D];
      int id = 0;
      for (int i =0; i<D+1; i++) {
        if (top1 == i) continue;
        simx[id] = new float[D];
        for (int j=0; j<D; j++) {
          simx[id][j] = sign * (x[top1][j] - x[i][j]);
        }
        id++;
      }
      
      float z_sim = sign * inner_product(x[top1], current_query, D);
      cout << "zsim" << " " << z_sim << endl;
      vector<operand> simplex_top;
      int bigcount = 0;
      for (int i = 0; i<D; i++) {
        float z_simplex_sim = cosine_similarity(simx[i], current_query, D, sqrt(squared_sum(simx[i],D)), 1);
        if (z_simplex_sim > z_sim) {
          bigcount++;
          
          vector<operand> simplex_top;
          for (int j=0; j<N; j++) {
            simplex_top.push_back(make_pair(j, cosine_similarity(simx[i], data[j], D, sqrt(squared_sum(simx[i], D)), 1)));
          }
          
          sort(simplex_top.begin(), simplex_top.end(), compare_descending);
          for (int j = 0; j<100; j++) {
            int idx = simplex_top[j].first;
            if (visited[idx] == 1) continue;
            visited[idx] = 1;
            candidate.push_back(make_pair(idx, getSimilarity(idx, current_query)));  
          }
        }
      } 
      cout << "big " << bigcount << endl;
      //simplex_top.push_back(make_pair(j, inner_product(simx[i], data[j], D)));
      //sort(simplex_top.begin(), simplex_top.end(), compare_descending);

      
      for (int i = 0; i<D; i++) {
        vector<operand> simplex_top;
        for (int j=0; j<N; j++) {
          simplex_top.push_back(make_pair(j, inner_product(simx[i], data[j], D)));
        }
        sort(simplex_top.begin(), simplex_top.end(), compare_descending);
        for (int j = 0; j<1024; j++) {
          int idx = simplex_top[j].first;
          if (visited[idx] == 1) continue;
          visited[idx] = 1;
          candidate.push_back(make_pair(idx, getSimilarity(idx, current_query)));  
        }
      }

      for (int i =0; i<D; i++) delete[] simx[i];
      delete[] simx;
    }
    cout << "Csize " << candidate.size() << endl; 
    partial_sort(candidate.begin(), candidate.begin()+K, candidate.end(), compare_descending);
    
    vector<operand> true_top;
    for (int i = 0; i<N; i++) {
      true_top.push_back(make_pair(i, getSimilarity(i,current_query)));
    }
    sort(true_top.begin(), true_top.end(), compare_descending);
    cout << check_recall(true_top, candidate) << endl;
     
*/






/*
    float* query = new float[D];
    for (int i = 0; i<D; i++) query[i] = 1;
    float norm = sqrt(squared_sum(query, D));
    for (int i = 0; i<D; i++) {
      query[i] /= norm;
      cout << query[i] << " ";
    }
    cout << endl;
*/
/*
    for (int i = 0; i<D; i++) query[i] = current_query[i];
    auto t1 = Clock::now();


    vector<int> visited(N);
    vector<operand> candidate;

    for (int i = 0; i<D; i++) {
      vector<operand> simplex_top;
      cout << simx[i] << endl;
      for (int j=0; j<N; j++) {
        simplex_top.push_back(make_pair(j, inner_product(x[i], data[j], D)));
      }
      sort(simplex_top.begin(), simplex_top.end(), compare_descending);
      for (int j = 0; j<100; j++) {
        int idx = simplex_top[j].first;
        if (visited[idx] == 1) continue;
        visited[idx] = 1;
        candidate.push_back(make_pair(idx, getSimilarity(idx, query)));  
      }
    }
    cout << candidate.size() << endl; 
    partial_sort(candidate.begin(), candidate.begin()+K, candidate.end(), compare_descending);

    vector<operand> true_top;
    for (int i = 0; i<N; i++) {
      true_top.push_back(make_pair(i, getSimilarity(i,query)));
    }
    sort(true_top.begin(), true_top.end(), compare_descending);
    cout << check_recall(true_top, candidate) << endl;
    



    //double** x = new double*[D];
    //for (int i = 0; i<D; i++) x[i] = new double[D+1];
   */ 

/*
    vector<operand> maxima;
    for (int i = 0; i<N; i++) {
      float sim = getSimilarity(i, current_query);
      bool localmaxima = true;
      for (int j = 0; j<256; j++) {
        if (getSimilarity(cg->graph[i][j], current_query) > sim) {
          localmaxima = false;
          break;
        }
      }
      if (localmaxima) {
        maxima.push_back(make_pair(i, getSimilarity(i, current_query)));
      }
    }
   
    for (int i = 0; i<10; i++) {
      cout << "Top" << i << " : " << true_top[i].first << " / " << true_top[i].second << endl;
    }
   
    float max = -1;
    for (int i = 0; i<500; i++) {
      int query_topi = true_top[i].first;
      vector<int>& query_topi_neighbor = cg->graph[query_topi];
      vector<operand> topi_operand;
      for (int qn : query_topi_neighbor) {
        topi_operand.push_back(make_pair(qn, 0));
      }
      float x = check_recall(true_top, topi_operand, false, true); 
      if (max < x) { max = x;} 
      cout << "Top" << i+1 << " " << getSimilarity(query_topi, current_query) << " : " << x << "%" << endl;
    }

  cout << max << endl;
*/
/*
    for (int k = 0; k<5; k++) {
      int source = maxima[k].first; //should be local maxima
      
      priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>> > pq;
      vector<int> dist(N, 1000000);
      vector<vector<int>> prev(N, vector<int>());

      pq.push(make_pair(0, source));
      dist[source] = 0;
      int edge = 32;
      while (!pq.empty()) {
        int u = pq.top().second;
        pq.pop();

        vector<int> neighborList = cg->graph[u];
        
        vector<operand> queryoperand;
        vector<int> queryidx;
        for(int i = 0; i<edge; i++) {
          queryoperand.push_back(make_pair(neighborList[i], getSimilarity(neighborList[i], current_query)));
        }
        sort(queryoperand.begin(), queryoperand.end(), compare_descending);
        for (auto x : queryoperand) queryidx.push_back(x.first);


        for (int i = 0; i<edge; i++) {
          int v = neighborList[i];
          ptrdiff_t pos = distance(true_top_idx.begin(), find(true_top_idx.begin(), true_top_idx.end(), v));
          ptrdiff_t querypos = distance(queryidx.begin(), find(queryidx.begin(), queryidx.end(), v));
          if(querypos > 9) {
            continue;
          }
          //if(pos > 80) {
          //  continue;
          //}
          int weight = 1;
          if (dist[v]>=dist[u]+weight) {
            if(dist[v]>dist[u]+weight) {
              dist[v] = dist[u]+weight;
              prev[v].clear();
            }
            prev[v].push_back(u);
            pq.push(make_pair(dist[v], v));
          }
        }
      }
      
      for (int i = 0; i<K; i++) {
        int u = true_top[i].first;
        cout << "Top" << i+1 << " : ";
        stack<pair<int, pair<int, int>>> S;
        if (prev[u].size()!=0 || u==source) {
          while (prev[u].size() != 0) {
            int prevu = -1;
            int ranku = 10000000;
            int edgeu = 100000;
            for (auto i : prev[u]) {
              vector<int> &neighborList = cg->graph[i];
              
              ptrdiff_t pos = distance(true_top_idx.begin(), find(true_top_idx.begin(), true_top_idx.end(), i));
              ptrdiff_t edgepos = distance(neighborList.begin(), find(neighborList.begin(), neighborList.end(), u));
              if (edgepos < edgeu) {
                ranku = pos;
                prevu = i;
                edgeu = edgepos;
              }
            }

            S.push(make_pair(prevu, make_pair(ranku, edgeu)));
            u = prevu;
          }

          while(!S.empty()) {
            int v = S.top().first;
            int localrank = S.top().second.second;
            vector<int>& neighborList = cg->graph[v];
            
            vector<operand> queryoperand;
            vector<float> querysim;
            for(int i = 0; i<edge; i++) {
              queryoperand.push_back(make_pair(neighborList[i], getSimilarity(neighborList[i], current_query)));
            }
            sort(queryoperand.begin(), queryoperand.end(), compare_descending);

            for(int i = 0; i<edge; i++) {
              querysim.push_back(queryoperand[i].first);
            }
            ptrdiff_t localrank_query = distance(querysim.begin(), find(querysim.begin(), querysim.end(), neighborList[localrank]));


            cout << "(" << S.top().first << ", " << (S.top().second.first)+1 << "th/" << (S.top().second.second) + 1 << "th/" << localrank_query + 1 << "th, " <<  getSimilarity(S.top().first, current_query) << " / " << getSimilarity(S.top().first, data[true_top[i].first]) << ") --> ";

            if(S.size() == 1) {
              cout << "(" << true_top[i].first << ", " << i+1 << "th, " << true_top[i].second << " / " << getSimilarity(true_top[i].first, data[true_top[i].first]) << ")";
            }
            S.pop();
          }

          cout << endl;
        }
       else {cout<<endl;}
      }
      cout << endl;
    }

    //walk(gx, current_query);
    //partial_sort(top_k_candidates.begin(), top_k_candidates.begin()+K, top_k_candidates.end(), compare_descending);
    //sort(top_k_candidates.begin(), top_k_candidates.end(), compare_descending);
    //for (int i = 0; i<10; i++) {
    //  cout << "MyTop" << i << " : " << top_k_candidates[i].first << " / " << top_k_candidates[i].second << endl;
    //}
*/
    auto t4 = Clock::now();
    timers[1] += compute_clock(t2, t1);
    timers[2] += compute_clock(t3, t2);
    timers[3] += compute_clock(t4, t3);
    return top_k_candidates;
  }
 /* 
  int walk(unordered_map<int, vector<int>> &g, float* current_query) {
    int seed = rand()%N;
    float seed_similarity = getSimilarity(seed, current_query);
    cout << "Seed : " << seed << " / Sim : " << seed_similarity << endl;
    int step = 40;

    while (step--) {
      vector<int>& neighborList = g[seed];
      vector<operand> cos_diff;

      vector<operand> new_sim;
      float** diff = new float*[neighborList.size()];
      for (int i = 0; i<neighborList.size(); i++) {
        diff[i] = new float[D];
        for (int j = 0; j <D; j++) {
          diff[i][j] = data[neighborList[i]][j] - data[seed][j];
        }
      }

      float* q_x = new float[D];
      for (int i = 0; i<D; i++) q_x[i] = current_query[i] - data[seed][i];



      for (int i =0; i<neighborList.size(); i++) {
        int neighbor = neighborList[i];
        float neighbor_similarity = getSimilarity(neighbor, current_query);
        cos_diff.push_back(make_pair(neighbor, neighbor_similarity));


        new_sim.push_back(make_pair(neighbor, cosine_similarity(q_x, diff[i], D, sqrt(squared_sum(q_x,D)), sqrt(squared_sum(diff[i],D)))));
      }

      partial_sort(cos_diff.begin(), cos_diff.begin() + 1, cos_diff.end(), compare_descending);
      partial_sort(new_sim.begin(), new_sim.begin() + 5, new_sim.end(), compare_descending);

      //seed = cos_diff[0].first;
      //cout << "Seed : " << seed << " / Sim : " << getSimilarity(seed, current_query) << endl;
      for (int i = 0; i<5; i++)
        cout << "NewSeed : " << new_sim[i].first << " / RelativeSim : " << new_sim[i].second << " / RealSim : " << getSimilarity(new_sim[i].first, current_query)  << endl;
      cout << endl;
      seed = new_sim[0].first;
    }
    return seed;
  }
*/
/*
  int walk(unordered_map<int, vector<int>> &g, float* current_query) {
    int seed = rand()%N; //originally random
    float seed_similarity = getSimilarity(seed, current_query);
    cout << "Seed : " << seed << " / Sim : " << seed_similarity << endl;
    int step = 40;
    vector<int> visited(N);

    while (step--) {
      vector<int>& neighborList = g[seed];
      vector<operand> cos_diff;

      for (int i =0; i<neighborList.size(); i++) {
        int neighbor = neighborList[i];
        float neighbor_similarity = getSimilarity(neighbor, current_query);
        cos_diff.push_back(make_pair(neighbor, neighbor_similarity));
      }

      partial_sort(cos_diff.begin(), cos_diff.begin() + (visited[seed]+1), cos_diff.end(), compare_descending);
      seed = cos_diff[visited[seed]].first;
      cout << "Seed : " << seed << " / Sim : " << getSimilarity(seed, current_query) << endl;
      visited[seed]++;
    }
    return seed;
  }
  */
  int walk(float* current_query) {
    int maxCheck = 64;
    int curSeed = rand()%N;
    float curSim = getSimilarity(curSeed, current_query);
    bool cont = true;
    int cnt = 0;
    while(cont) {
      cont = false;
      vector<int> &neighborList = (*cg->graph)[curSeed];
      for (size_t i = 0; i < min(maxCheck, (int)neighborList.size()); i++) {
        float neighborSim = getSimilarity(neighborList[i], current_query);
        searched++;
        if(neighborSim > curSim) {
          curSim = neighborSim;
          curSeed = neighborList[i];
          cont = true;
        }
      }
    }
    //cout << "walk : " << cnt << endl;
    return curSeed;
  }

  void print_time() {
    cout << "[GraphSearchWJY] " << (timers[1] + timers[2]) /queryCount << "(ms)\n";
    cout << "Preprocessing : " << timers[0] << "(ms)\n";
    cout << "Search " << (timers[1] + timers[2]) / queryCount << "(ms)\n";
    cout << "---[1] Walk : " << timers[1] / queryCount << "(ms)\n";
    cout << "---[2] Expand : " << timers[2] / queryCount << "(ms)\n";
    cout << "---[2-1] Expand-Main : " << timers[3] / queryCount << "(ms)\n";
    cout << "---[2-2] Expand-Sort : " << timers[4] / queryCount << "(ms)\n";
    cout << "local count : " << localcount / (float)queryCount << endl;
    cout << "local recall : " << localrecall / (float)queryCount << endl;
    cout << "localneighbor count : " << neighborcount / (float)queryCount << endl;
    cout << "localneighbor recall : " << neighborrecall / (float)queryCount << endl;

   }
  void print_information() {
    cout << "[GraphSearch] \n";
  }

float** simplex(int D) {
    double *x = new double[D*(D+1)];
    for (int i = 0; i<D*(D+1); i++) { x[i] = 0; }
    
    double a;
    double c;
    int i;
    int j;
    double s;
    int n = D;

    for ( i = 0; i < n; i++ )
    {
      x[i+i*n] = 1.0;
    }

    a = ( 1.0 - sqrt ( 1.0 + ( double ) ( n ) ) ) / ( double ) ( n );

    for ( i = 0; i < n; i++ )
    {
      x[i+n*n] = a;
    }
    for ( i = 0; i < n; i++ )
    {
      c = 0.0;
      for ( j = 0; j < n + 1; j++ )
      {
        c = c + x[i+j*n];
      }
      c = c / ( double ) ( n + 1 );
      for ( j = 0; j < n + 1; j++ )
      {
        x[i+j*n] = x[i+j*n] - c;
      }
    }
    s = 0.0;
    for ( i = 0; i < n; i++ )
    {
      s = s + x[i+0*n] * x[i+0*n];
    }
    s = sqrt ( s );

    for ( j = 0; j < n + 1; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        x[i+j*n] = x[i+j*n] / s;
      }
    }
    
    float **newx = new float*[D+1];
    for (int i = 0; i<D+1; i++) newx[i] = new float[D];

    for (int i = 0; i<D+1; i++) {
      for (int j =0; j<D; j++) {
        newx[i][j] = (float)((double)x[i*D+j] - (double)x[D*D+j]);
        //newx[i][j] = (float)((double)x[i*D+j]);
      }
    }
    delete[] x;

    for (int i = 0; i<D; i++) {
      float sqsum = sqrt(squared_sum(newx[i],D));
      for (int j =0; j<D; j++) {
        newx[i][j] /= sqsum;
      }
    }
    return newx;
}

};



/*
int source = maxima[0].first; //should be local maxima
    
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>> > pq;
    vector<int> dist(N, 1000000);
    vector<vector<int>> prev(N, vector<int>());

    pq.push(make_pair(0, source));
    dist[source] = 0;
    int edge = 512;
    while (!pq.empty()) {
      int u = pq.top().second;
      pq.pop();

      vector<int> neighborList = cg->graph[u];
      for (int i = 0; i<edge; i++) {
        int v = neighborList[i];
        ptrdiff_t pos = distance(true_top_idx.begin(), find(true_top_idx.begin(), true_top_idx.end(), v));
        if(pos > 80) {
          continue;
        }
        int weight = 1;
        if (dist[v]>=dist[u]+weight) {
          if(dist[v]>dist[u]+weight) {
            dist[v] = dist[u]+weight;
            prev[v].clear();
          }
          prev[v].push_back(u);
          pq.push(make_pair(dist[v], v));
        }
      }
    }
    
    for (int i = 0; i<K; i++) {
      int u = true_top[i].first;
      cout << "Top" << i+1 << " : ";
      stack<pair<int, pair<int, int>>> S;
      if (prev[u].size()!=0 || u==source) {
        while (prev[u].size() != 0) {
          int prevu = -1;
          int ranku = 10000000;
          int edgeu = -1;
          for (auto i : prev[u]) {
            vector<int> &neighborList = cg->graph[i];
            
            ptrdiff_t pos = distance(true_top_idx.begin(), find(true_top_idx.begin(), true_top_idx.end(), i));
            ptrdiff_t edgepos = distance(neighborList.begin(), find(neighborList.begin(), neighborList.end(), u));
            if (pos < ranku) {
              ranku = pos;
              prevu = i;
              edgeu = edgepos;
            }
          }

          S.push(make_pair(prevu, make_pair(ranku, edgeu)));
          u = prevu;
        }

        while(!S.empty()) {
          int v = S.top().first;
          int localrank = S.top().second.second;
          vector<int>& neighborList = cg->graph[v];
          
          vector<operand> queryoperand;
          vector<float> querysim;
          for(int i = 0; i<edge; i++) {
            queryoperand.push_back(make_pair(neighborList[i], getSimilarity(neighborList[i], current_query)));
          }
          sort(queryoperand.begin(), queryoperand.end(), compare_descending);

          for(int i = 0; i<edge; i++) {
            querysim.push_back(queryoperand[i].first);
          }
          ptrdiff_t localrank_query = distance(querysim.begin(), find(querysim.begin(), querysim.end(), neighborList[localrank]));


          cout << "(" << S.top().first << ", " << (S.top().second.first)+1 << "th/" << (S.top().second.second) + 1 << "th/" << localrank_query + 1 << "th, " <<  getSimilarity(S.top().first, current_query) << ") --> ";

          if(S.size() == 1) {
            cout << "(" << true_top[i].first << ", " << i+1 << "th, " << true_top[i].second << ")";
          }
          S.pop();
        }

        cout << endl;
      }
     else {cout<<endl;}
    }
    cout << endl;
*/






class SCCAnalysis {
public:
  struct SCC_struct {
    vector<vector<int>> adj; 
    vector<int> discovered, finished;
    stack<int> st;
    vector<int> sccId;
    int sccCounter, vertexCounter;
  };
  static int scc(SCC_struct& scs, int here) { 
    int ret = scs.discovered[here] = scs.vertexCounter++;
    scs.st.push(here);  
    for (size_t i = 0; i < scs.adj[here].size(); ++i) { 
      int there = scs.adj[here][i];
      if (scs.discovered[there] == -1) 
        ret = min(ret, scc(scs, there));
      else if (scs.discovered[there] < scs.discovered[here] && scs.sccId[there] == -1) 
        ret = min(ret, scs.discovered[there]);
    }  
    if (ret == scs.discovered[here]) {
      while(true) {
        int t = scs.st.top();
        scs.st.pop();
        scs.sccId[t] = scs.sccCounter;
        if (t == here) 
          break;
      }
      ++scs.sccCounter;
    }
    scs.finished[here] = 1;
    return ret;
  }
  static map<int,int> tarjanSCC(unordered_map<int, vector<int>> &graph, int &sccCount) {
    SCC_struct scs;
    scs.adj = vector<vector<int>>(graph.size());
    map<int, int> iid;
    map<int, int> iid_rev;
    int cnt = 0;
    for (auto itr = graph.begin(); itr != graph.end(); itr++) {
      int i = itr->first;
      iid_rev.insert(make_pair(cnt, i));
      iid.insert(make_pair(i, cnt++));
    }
    for (auto itr = graph.begin(); itr != graph.end(); itr++) {
      int i = itr->first;
      vector<int> neighbor = itr->second;
      for (size_t j = 0; j < neighbor.size(); j++)
        scs.adj[iid[i]].push_back(iid[neighbor[j]]);
    }
    scs.sccId = scs.discovered = vector<int>(scs.adj.size(), -1);
    scs.finished = vector<int>(scs.adj.size(), 0);
    scs.vertexCounter = 0;
    scs.sccCounter = 0;
    for (size_t i = 0; i < scs.adj.size(); i++) {
      if (scs.discovered[i] == -1)
        scc(scs, i);
    }
    sccCount = scs.sccCounter;
    map<int, int> sccId_decoded;
    for (int i = 0; i < cnt; i++) {
      int org_idx = iid_rev[i];
      sccId_decoded.insert(make_pair(org_idx, scs.sccId[i]));
    }
    if(sccCount == 0)
      assert(false);
    return sccId_decoded;
  }
};



