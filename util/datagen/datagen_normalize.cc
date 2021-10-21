#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <cstring>
#include <sstream>

using namespace std;

/* Computation */
float inline squared_sum(float* __restrict__ in, int d) {
  float sum = 0.0f;
  #pragma omp simd reduction(+:sum)
  for (int i = 0; i < d; i++)
    sum += in[i] * in[i];
  return sum;
}
int main(int argc, char *argv[]) {
    
    if(argc != 2) {
      cout << "Usage: ./datagen_normalize dataname \n";
      return 0;
    }
    string data_path = "./data/" + string(argv[1]) + ".txt";
    string query_path = "./data/" + string(argv[1]) + "_query.txt";
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
    int N, D, qsize, qN;
    float** data;
    float** query;
    int q;
    sts >> N;
    sts >> D;
    qsts >> q;
    int dummy = 0;
    qsts >> dummy;
    qN = q;
    cout << "N : " << N << " / D : " << D << " / qN : " << q << "\n";
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
      float norm = sqrt(squared_sum(data[i], D));
      for (int j = 0; j < D; j++)
        data[i][j] /= norm;
    }
    /* Read Query */
     for (int i = 0; i < q; i++) {
      getline(fs_query, qline);
      stringstream sts(qline);
      for (int j = 0; j < D; j++)
        sts >> query[i][j];
      
      float norm = sqrt(squared_sum(query[i], D));

      for (int j = 0; j < D; j++)
        query[i][j] /= norm;
    }
    fs_data.close();
    fs_query.close();

    ofstream ofile, qfile;
    ofile.open("./data/c"+ string(argv[1]) +".txt", ofstream::out | ofstream::trunc);
    qfile.open("./data/c" + string(argv[1]) + "_query.txt", ofstream::out | ofstream::trunc);
    if(!ofile.is_open() || !qfile.is_open())
      exit(1);

    // Print
    int n = N;
    int d = D;
    ofile << n << " " << d << "\n";
    for(int i=0; i<n; i++) {
      for (int j=0; j<d; j++) {
        ofile << data[i][j];
        if(j != d-1)
          ofile << " ";
      }
      if(i != n-1)
        ofile << "\n";
    }
    qsize = q;
    qfile << qsize << " " << d << "\n";
    for(int i=0; i<qsize; i++) {
      for (int j=0; j<d; j++) {
        qfile << query[i][j];
        if(j != d-1)
          qfile << " ";
      }
      if(i != qsize-1)
        qfile << "\n";
    }
}
