//Resources and the Dynamic-programming Dynamic-array Dynamic solution to the subset_sum problem
//Jacob Liddy
//
#pragma once
#include <map>
#include <vector>
#include <NTL/ZZ.h>
#include "tools.hpp"

using std::vector;
using std::map;
using NTL::ZZ;

struct BoolVector{
  BoolVector(){}

  bool defined(ZZ pos){
    //std::cout << "Defined? "  << (lookup.find(pos) != lookup.end()) << '\n';
    return lookup.find(pos) != lookup.end();
  }

  bool& operator[](ZZ pos){
    if (lookup.find(pos) == lookup.end()) {
      lookup[pos] = false;
    }
    return lookup[pos];
  }

  map<ZZ, bool> lookup; 
};


struct BoolMatrix{

  BoolMatrix(size_t len): rows(len) {}
  
  bool defined(size_t i, ZZ j){return rows[i].defined(j);}

  BoolVector& operator[](size_t i){ return rows[i]; }

  vector<BoolVector> rows;
};

bool isSubsetSum2(vector<ZZ>& set, ZZ sum);
