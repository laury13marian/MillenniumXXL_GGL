#ifndef GENERAL_H
#define GENERAL_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include <ctime>
#include <new>
#include <algorithm>


using namespace std;

inline void check_stream(ifstream &in, string file_name) 
  
{ 
  if (!in) {
    cerr << "Cannot open file " + file_name << endl;
    exit(1);
  }
  else 
    cerr<<endl<<"reading file " + file_name<<endl; 

}

#endif
