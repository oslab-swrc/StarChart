// SPDX-FileCopyrightText: Copyright Copyright 2021 Seoul National University
//
// SPDX-License-Identifier: MIT License
class Manager {
public:
  vector<BaseSearch*> objs;
  DataManager* dm;
  Evaluator* eval;
  map<string, vector<int>> param_map;
  int queryTimeCounter[100];
  Manager() { }
  ~Manager() { }
  void init(string data_path, string data_name, string query_name, int K, int nSearch, int qN, int s_size, float threshold) {
    dm = new DataManager(data_path, data_name, query_name, qN);
    eval = new Evaluator("Evaluator", dm->data, dm->query, dm->N, dm->D, K, nSearch);
    GraphConstruct::init(dm->N, dm->D, dm->data, s_size, threshold);
    BaseSearch::init(dm->N, dm->D, K, dm->data);
    for(int i = 0; i < 100; i++)
      queryTimeCounter[i] = 0;
  }
  void parse_param(string process_mode, int& K, int& nSearch, int& qN, int& construct_fg_M, int& load_fg_M, int& load_front_M, int& ingraph_num_node, int& ingraph_M, float& threshold, int& s_size) {
    if (process_mode == "-d") {
      K = 100;
      nSearch = K * 10;
      qN = 100;
    }
    else if (process_mode == "-c") {
      fstream cfile;
      cfile.open("./config.txt", fstream::in);
      check_file(cfile, "./config.txt");
      string line;
      for (int i = 0; i < 10; i++) {
        getline(cfile, line);
        stringstream sts(line);
        if (i == 0) 
          sts >> K;
        else if (i == 1)
          sts >> nSearch;
        else if (i == 2)
          sts >> qN;
        else if (i == 3)
          sts >> construct_fg_M;
        else if (i == 4)
          sts >> load_fg_M;
        else if (i == 5)
          sts >> load_front_M;
        else if (i == 6)
          sts >> ingraph_num_node;
        else if (i == 7)
          sts >> ingraph_M;
        else if (i == 8)
          sts >> threshold;
        else if (i == 9)
          sts >> s_size;
        else
          assert(false);
      }
      cfile.close();
    }
    else
      assert(false);  
  }
  void path_register(string name, string path) {
    dm->set_paths(name, path);
  }
  void objs_register(int args, ...) {
    va_list argc;
    va_start(argc, args);
    for (int i = 0; i < args; i++) {
      BaseSearch* obj = va_arg(argc, BaseSearch*);
      objs.push_back(obj);
      if (param_map.find(obj->name) != param_map.end()) {
        obj->param_register(param_map.find(obj->name)->second);
      }
    }
    va_end(argc);
  }

  void objs_register_for_query_level_par(BaseSearch* object) {
    // va_list argc;
    // va_start(argc, args);
    // BaseSearch* obj = va_arg(argc, BaseSearch*);
    // for (int i = 0; i < qN; i++) {
      objs.push_back(object);
    // }
    // va_end(argc);
  }
  void run(int qN) {
    vector<vector<operand>> resArray(objs.size(), vector<operand>());
    #ifndef FIGURE
    #ifndef FIGUREOUR
    eval->run(dm->get_path("Groundtruth"), qN);
    for (size_t i = 0; i < objs.size(); i++)
      objs[i]->preprocess();
    float **myQuery = dm->query;
    cout << "obj size : " << objs.size() << "\n";
    float recall = 0.0f;
    for (size_t i = 0; i < objs.size(); i ++)
      eval->set_evaluator(objs[i]->name);
    auto t1 = Clock::now();
    #pragma omp parallel for num_threads(6)
    for (size_t i = 0; i < qN; i ++) {
      // cout << i << "th query\n";
      resArray[i] = vector<operand>(objs[i]->search(myQuery[i]));
      // for (int j = 0; j < 5; j++)
      //   cout << "(" << resArray[i][j].first << ", " << resArray[i][j].second << ") ";
      // cout << "\n";
    }
    auto t2 = Clock::now();
    for (size_t i = 0; i < objs.size(); i ++)
      objs[i]->update();
    for (size_t i = 0; i < objs.size(); i ++)
      recall += eval->update(eval->get_evaluator(objs[i]->name), resArray[i]);
    float timer = compute_clock(t2, t1);
    cout << "Query level parallelism search time : " << timer << "\n";
    cout <<  "Query level parallelism recall : " << recall / qN << "\n";
    return;
    #endif
    #endif

    #ifdef FIGURE
    #ifndef OURFIGURE
    eval->run(dm->get_path("Groundtruth"), qN);
    int iter = 1400;
    int cnt = 0;
    for (size_t ef = 50; ef < iter; ef+=50) {
      ef_ = ef;
      for (size_t i = 0; i < objs.size(); i++)
        objs[i]->preprocess();
      for (size_t i = 0; i < objs.size(); i ++) {
        string hnswname = objs[i]->name + to_string(ef);
        eval->set_evaluator(hnswname);
        float **myQuery = dm->query;
        for (int j = 0; j < qN; j++) {
          resArray[i][j] = vector<operand>(objs[i]->search(myQuery[j]));
          objs[i]->update();
          eval->update(eval->get_evaluator(hnswname), resArray[i][j]);
        }
        cout << eval->recall[cnt]/eval->queryCount[cnt] << " " << objs[i]->timers[1] / eval->queryCount[cnt] << " " << objs[i]->totalSearched/objs[i]->queryCount  <<  endl;
      }
      cnt++;
    }
    return;
    #endif
    #endif

    #ifdef FIGUREOUR
    int ITER = 5;
    eval->recall.clear();
    eval->queryCount.clear();
    vector<int> K_list = {1, 5};
    vector<int> nSearchList = {1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 12000, 13000, 14000, 15000, 16000, 17000, 18000, 20000, 22000, 24000, 25000, 26000, 28000, 30000, 32000, 34000, 36000, 38000, 40000};
    for (size_t m = 0; m < K_list.size(); m++) {
      int topK = K_list[m];
      eval->K = topK;
      eval->run(dm->get_path("Groundtruth"), qN);
      cout << "Top " << topK << " report\n";
      for (int l = 0; l < ITER; l++) {
        for (size_t k = 0; k < nSearchList.size(); k++) {
          nSearch_ = nSearchList[k];
          K_ = topK;
          for (size_t i = 0; i < 1; i++)
            objs[i]->preprocess();
          for (size_t i = 0; i < 1; i ++) {
            string hnswname = objs[i]->name + "top" + to_string(topK) + "_nSearch" + to_string(nSearch_) + "_" + to_string(l);
            eval->set_evaluator(hnswname);
            float **myQuery = dm->query;
            for (int j = 0; j < qN; j++) {
              resArray[i][j] = vector<operand>(objs[i]->search(myQuery[j]));
              objs[i]->update();
              eval->update(eval->get_evaluator(hnswname), resArray[i][j]);
            }
            int idx = nSearchList.size() * ITER * m + nSearchList.size() * l + k;
            cout << eval->recall[idx]/eval->queryCount[idx] << " / " << objs[i]->timers[1] / eval->queryCount[idx] << " / " << objs[i]->totalSearched/objs[i]->queryCount  <<  endl;
          }
        }
        cout << "\n";
      }
    }
    #endif

  }
  void print_result() {
    cout << "\n" << "[Search count] \n";
    for (size_t i = 0; i < objs.size(); i++)
      objs[i]->print_result();
    cout << "\n" << "[Timing] \n";
    // for (size_t i = 0; i < objs.size(); i++){
    //   objs[i]->print_time();
    //   cout << "\n";
    // }
    // cout << "\n" << "[Statistics] \n";
    // for (size_t i = 0; i < objs.size(); i++)
    //   objs[i]->print_information();
    cout << "\n" << "[Evaluations] \n";
    eval->print_result();
  }
};
