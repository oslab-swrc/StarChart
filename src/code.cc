#define OMD (56)
#include "libs.hh"
#include "evaluator.cc"
#include "search/graph.cc"
#include "search/base.cc"
#include "search/graphSearch.cc"
#include "search/tjSearch.cc"
#include "search/streamTracker.cc"
#include "search/yjSearch.cc"
#include "search/emulatedNIPS.cc"
#include "manager.cc"
#ifdef FALCONN
#include "search/falconnSearch.cc"
#endif

void parse_and_save(string data_path, int N, int D, int partitionNum) {
  fstream gfile;
  gfile.open(data_path + ".txt", fstream::in);
  check_file(gfile, data_path);
  uint32_t count = 0;
  uint32_t per_partition = N / partitionNum;
  // assert(N % partitionNum == 0);
  cout << "perPartiton : " << per_partition << "\n";
  string line;
  vector<string> temp_line;
  getline(gfile, line);
  while(getline(gfile, line)) {
    count++;
    temp_line.push_back(line);
    if (count % per_partition == 0) {
      ofstream ofile;
      ofile.open(data_path + "_" + to_string(count/per_partition) + ".txt");
      ofile << per_partition << " " << D << "\n";
      cout << "Writing to file : " << data_path + "_" + to_string(count/per_partition) + ".txt" << "\n";
      for (auto it : temp_line) {
        ofile << it << "\n";
      }
      temp_line.clear();
      ofile.close();
    }
  }
  gfile.close();
}
int main (int argc, char** argv) {
  int K, nSearch, qN;
  int construct_fg_M, load_fg_M, load_front_M, ingraph_num_node, ingraph_M, s_size;
  float threshold;
  int numPart = 8;
  Manager manager;

  /* Wrong Usage */
  if (argc != 4) {
    cout << "Usage:\n";
    cout << "./bin/ssa -s -c dataName -> save data" << endl;
    cout << "./bin/ssa -l (-c || -d || -i) dataName (queryname) -> load data and process\n";
    cout << "./bin/ssa -g (-c || -d || -i) dataName (queryname) -> generate data and process\n";
    cout << "./bin/ssa -parse #part dataName -> save data" << endl;
    exit(0);
  }
  cout << "Check the Current Core Count : " << OMD << "\n";
  cout << "Make sure it is the half of the CPUs (check with lscpu) \n";

  string data_mode = argv[1];
  string process_mode = argv[2];
  string data_name = argv[3];
  string query_name = argv[3];
    
  cout << "numPart : " << numPart << "\n";
  string path = "./data/";

  if (data_mode == "-parse") {
    // parse_and_save(path + data_name, 1120000, 100, stoi(process_mode));
    parse_and_save(path + data_name, 1917488, 300, stoi(process_mode));
    exit(0);
  }

  manager.parse_param(process_mode, K, nSearch, qN, construct_fg_M, load_fg_M, load_front_M, ingraph_num_node, ingraph_M, threshold, s_size);
  manager.init(path, data_name, query_name, K, nSearch, qN, s_size, threshold);
  assert(manager.dm->N % numPart == 0);
  
  /* Data Save & Load */
  if (data_mode == "-s" || data_mode == "-fg2sfg" || data_mode == "-sfg2sfgbl" || data_mode == "-sfgbl2sigbl" || data_mode == "-sigbl2sbligbl" || data_mode == "-sbligbl2sbligblc"
    || data_mode == "-elsfgbl"
    || data_mode == "-olblc"
    || data_mode == "-hl") {   
    manager.path_register("Groundtruth",path + "groundtruth/" + manager.dm->name + "_gt.txt");  
    path += "preprocess/";
    manager.path_register("HNSW", path + manager.dm->name + ".hnsw");
    manager.path_register("SFG", path + manager.dm->name + "_fg_sparsify.txt");
    manager.path_register("SFGBL", path + manager.dm->name + "_fg_sparsify_bidirectional_light.txt");
    manager.path_register("FG", path + manager.dm->name + "_fg.txt");
    for (int i = 0; i < numPart; i++)
      manager.path_register("SIGBL" + to_string(i), path + manager.dm->name + "_ig" + to_string(i) + "_w_sparsified_bidirectional_light_fg.txt");
    for (int i = 0; i < numPart; i++)
      manager.path_register("SBLIGBL" + to_string(i), path + manager.dm->name + "_bidirectionalized_light_ig" + to_string(i) + "_w_sparsified_bidirectional_light_fg.txt");
    for (int i = 0; i < numPart; i++)
      manager.path_register("SBLIGBLC" + to_string(i), path + manager.dm->name + "_bidirectionalized_light_cut_ig" + to_string(i) + "_w_sparsified_bidirectional_light_fg.txt");
  }
  else {
    assert(false);
  }

  GraphSearch gs("GS", K, nSearch);
  // HnswSearch hnswgs("HNSW", K, nSearch);
  // GraphSearchYJBase gsyjb1("GSYJBase1", K, nSearch, load_front_M, 192, 256);
  GraphSearchFinal gsf1("GSF1", K, nSearch, load_front_M);
  GraphSearchYJBase* gsyjb1[qN];
  for (int i = 0; i < qN; i++) 
    gsyjb1[i] = new GraphSearchYJBase("GSYJBase" + to_string(i), K, nSearch, load_front_M, 192, 256);
  #ifdef FALCONN
  FalconnSearch falconn("Falconn", K, nSearch);
  #endif

  cout << "[Parameters]\n";
  cout << "---Search\n";
  cout << "Search for " << K << " items out of " << manager.dm->N << " total items each having " << manager.dm->D << " dimensions. \n";
  cout << "Checks up to " << nSearch << " items \n";
  cout << "construct_fg_M : " << construct_fg_M << " / load_fg_M : " << load_fg_M << " / load_front_M : " << load_front_M << " / ingraph_num_node : " << ingraph_num_node << " / ingraph_M : " << ingraph_M << " / threashold : " << threshold << " / s_size : " << s_size << "\n";

 if(data_mode.find("-el") != string::npos)
    manager.objs_register(1, &gsf1);
  else if(data_mode.find("-ol") != string::npos) {
    // manager.objs_register(1, &gsyjb1);
    for (int i = 0; i < qN; i++) {
      manager.objs_register_for_query_level_par(gsyjb1[i]);
    }
  }
  // else if(data_mode.find("-hl") != string::npos)
  //   manager.objs_register(1, &hnswgs);

  vector<vector<int>> fg;
  vector<vector<int>> graph;
  vector<flat_hash_map<int, vector<int>>> ig;

  if (data_mode == "-s") {
    manager.eval->run("", qN);
    manager.eval->save_groundtruth(manager.eval->true_top_k, manager.dm->get_path("Groundtruth"));
    // vector<vector<int>> graph = GraphConstruct::construct(construct_fg_M, false, inner_product);
    // hnswlib::HierarchicalNSW<float> *hnsw_graph = GraphConstruct::construct_hnsw();    
    // GraphConstruct::save_graph(construct_fg_M, graph, manager.dm->get_path("FG"));
    // GraphConstruct::save_hnsw(hnsw_graph, manager.dm->get_path("HNSW"));
    exit(0);
  }
  else if (data_mode == "-fg2sfg") {
    long long int tote = 0;
    long long int tott = 0;
    auto fg = GraphConstruct::load_graph(manager.dm->get_path("FG"), load_fg_M);
    #pragma omp parallel for schedule(dynamic) num_threads(OMD)
    for (int i = 0; i < fg.size(); i++) {
      if (i % 1000 == 0) {
      	#pragma omp critical
      	{
              cout << "Sparsify progress : " << (float)i / fg.size() * 100 << "%, average edge so far : " << (float)tote/tott << "\n";
      	}
      }
      fg[i] = GraphConstruct::sparsify(fg[i], i);
      tote += fg[i].size();
      tott++;
    }
    GraphConstruct::save_graph(-1, fg, manager.dm->get_path("SFG"));
    cout << "Average edge : " << (float)tote / tott << "\n";
    exit(0);
  }
  else if (data_mode == "-sfg2sfgbl") {
    auto fg = GraphConstruct::load_graph(manager.dm->get_path("SFG"));
    GraphConstruct::bidirectionalize_light(fg);
    GraphConstruct::save_graph(-1, fg, manager.dm->get_path("SFGBL"));
    exit(0);
  }
  else if (data_mode == "-sfgbl2sigbl") {
    auto sfg = GraphConstruct::load_graph(manager.dm->get_path("SFGBL"));
    for (int i = 0; i < numPart; i++) {
      auto ig = GraphConstruct::construct_partial_graph_in_graph(i, manager.dm->N/numPart, ingraph_M, sfg, inner_product);
      GraphConstruct::save_partial_in_graph(i, manager.dm->N/numPart, ig, manager.dm->get_path("SIGBL" + to_string(i)), ingraph_num_node, ingraph_M);
    }  
    exit(0);
  }
  else if (data_mode == "-sigbl2sbligbl") {
    for (int i = 0; i < numPart; i++) {
      GraphConstruct::load_partial_in_graph_for_bidirection(i, manager.dm->N / numPart, ingraph_num_node, ingraph_M, manager.dm->get_path("SIGBL" + to_string(i)), ig, true);
      #pragma omp parallel for num_threads(OMD)
      for (int i = 0; i < ig.size(); i++) {
        if (i % 1000 == 0) {
          #pragma omp critical
          {
          cout << "Ingraph " << i << "'s Bidirectionalize progress : " << (float)i / ig.size() * 100 << "%\n";
          }
        }
        GraphConstruct::bidirectionalize_light(ig[i]);
      }
      GraphConstruct::save_partial_in_graph(i, manager.dm->N/numPart, ig, manager.dm->get_path("SBLIGBL" + to_string(i)), ingraph_num_node); 
      ig.clear();
    }
    exit(0);
  }
  else if (data_mode == "-sbligbl2sbligblc") {
    auto sfgbl = GraphConstruct::load_graph(manager.dm->get_path("SFGBL"));

    for (int i = 0; i < numPart; i++) {
      GraphConstruct::load_partial_in_graph_for_bidirection_and_cut(i, manager.dm->N / numPart, ingraph_num_node, load_front_M, manager.dm->get_path("SBLIGBL" + to_string(i)), ig, sfgbl);
      GraphConstruct::save_partial_in_graph(i, manager.dm->N/numPart, ig, manager.dm->get_path("SBLIGBLC" + to_string(i)), ingraph_num_node); 
      ig.clear();
    }
    exit(0);
  }
  
  if (data_mode != "-hl"){    // -l, -el, -ol
    if (data_mode == "-elsfgbl" || data_mode == "-olblc")
      fg = GraphConstruct::load_graph(manager.dm->get_path("SFGBL"));
    // ig.reserve(manager.dm->N);
    // if (data_mode == "-olblc") {
    //   for (int i = 0; i < numPart; i++)
    //     GraphConstruct::load_partial_in_graph(i, manager.dm->N / numPart, ingraph_num_node, -1, manager.dm->get_path("SBLIGBLC" + to_string(i)), ig);
    // }
  }

  if (data_mode.find("-el") != string::npos) {
    gsf1.fg->graph = &fg;
  }
  else if (data_mode.find("-ol") != string::npos) {
    for (int i = 0; i < qN; i++) {
      gsyjb1[i]->fg->graph = &fg;
      gsyjb1[i]->in_graph = &ig;
    }
  }
  // else if (data_mode.find("-hl") != string::npos) {
  //   hnswgs.hnswg = GraphConstruct::load_hnsw(manager.dm->get_path("HNSW"));
  // }
  #ifdef CACHE
  system("sync; echo 1 > /proc/sys/vm/drop_caches");
  cout << "Cache cleared.\n";
  #endif
  manager.run(qN);
  #ifdef CACHE
  system("sync; echo 1 > /proc/sys/vm/drop_caches");
  cout << "Cache cleared.\n";
  #endif

  // manager.print_result();
  return 0;
}
