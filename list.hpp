//To house the list structure
#pragma once

#define TEST

#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <vector>
#include <iostream>
#include <utility>

//#include "tools.hpp"

using std::vector;
using NTL::RR;
using NTL::ZZ;
using std::pair;
using std::make_pair;

struct List{

  List(vector<ZZ> primes);

  int operator[](int index){ return storage[index].second; }
  RR  delta(int index){ return storage[index].first; }
  RR  delta(){return storage[position].first;}
  int index(){return storage[position].second;}
  int size(){return storage.size();}
  int at(){return position;}
  List& operator++(){ ++position; return *this;}

  void resort(ZZ& p, ZZ& q);

  int position;
  vector<pair<RR, int> > storage;    
};
void display(List& l);
