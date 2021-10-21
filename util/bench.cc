#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include "../src/libs.hh"

using phmap::flat_hash_map;
using phmap::flat_hash_set;
#include "../src/search/streamTracker.cc"

using namespace std;

int main(int argc, char *argv[]) {
  StreamTracker st(128);
  vector<float> vals;
  vector<operand> ops;
  vector<int> ids;
  srand(7777);
  for(int i=0; i<500000; i++) {
    vals.push_back(rand() % 1000000 * 1.124f);
  }
  cout << "Mid \n";
  auto t1 = Clock::now();
  for(int i=0; i<100000; i++) {
    int id = rand()%10;
    st.update(id, vals[i*2]);
    id = rand()%10;
    st.update(id, vals[i*2+1]);
    operand k = st.get();
    ops.push_back(k);
  }
  auto t3 = Clock::now();

  cout << "---[1] Update & Get : " << compute_clock(t3, t1) << "(ms)\n";

  // for(int i=0; i<100; i++) {
  //   cout << ops[i].first << "\n";
  //   cout << ops[i].second << "\n";
  // }
  return 0;
}
