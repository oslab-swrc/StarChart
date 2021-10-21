#include <iostream>
#include <fstream>
#include <string>
#include <random>

using namespace std;

int main(int argc, char *argv[]) {
    
    if(argc != 3 && argc != 4) {
      cout << "Check the command line inputs. 3 inputs are required but " << argc << " inputs are received \n";
      cout << "Usage: ./datagen n d query#\n";
      return 0;
    }
    ofstream ofile;
    ofstream qfile;

    random_device rd;
    mt19937 generator(rd());
    int n = stoi(argv[1]);
    int d = stoi(argv[2]);
		int qsize = stoi(argv[3]);
    //string dist = argv[3];
    //int qsize = 1;
    if (qsize > n)
      qsize = n;

    // Create Distribution
    normal_distribution<> dq{0,1};
    normal_distribution<> d1{2,1};
    normal_distribution<> d2{-2,1};

    string path = "data/";
    string p = "";
    p += to_string(n) + "_" + to_string(d) + "_" + "hetero";
    cout << "Generating Distribution " << p << "\n";
    ofile.open(path + p + ".txt", ofstream::out | ofstream::trunc);
    qfile.open(path + p + "_query.txt", ofstream::out | ofstream::trunc);
    if(!ofile.is_open() || !qfile.is_open())
      exit(1);
    float **data = new float*[n];
    for (int i=0; i<n; i++) {
      data[i] = new float[d];
    }
    float **query = new float*[n];
    for (int i=0; i<qsize; i++) {
      query[i] = new float[d];
    }
    // Sample Distribution
    for(int i=0; i<n; i++) {
      for (int j=0; j<d; j++) {
        if(i % 2 == 0)
          data[i][j] = d1(generator);
        else if(i % 2 == 1)
          data[i][j] = d2(generator);
      }
    }
    for(int i=0; i<qsize; i++) {
      for (int j=0; j<d; j++) {
        query[i][j] = dq(generator);
      }
    }
    // Print
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

    for(int i=0; i<qsize; i++)
      delete[] query[i];
    for(int i=0; i<n; i++)
      delete[] data[i];
    delete[] data;
    delete[] query;
    qfile.close();
    ofile.close();
    return 0;
}
