//To house the list structure
#include "list.hpp"
#define TEST

#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <vector>
#include <iostream>
#include <utility>

#include "tools.hpp"

using std::vector;
using NTL::RR;
using NTL::ZZ;
using std::pair;
using std::make_pair;

List::List(vector<ZZ> primes){
    int i = 0;
    position = 0;

    for (p : primes){
      storage.push_back(make_pair(Delta(p), i));
      ++i;
    }
  }

void List::resort(ZZ& p, ZZ& e){

  RR d = Delta(p, e);

  int index = storage[position].second;
  storage[position].first = d;

  if (position == size() - 1){ return; }
  if (storage[position + 1].first < d){ return; }

  storage.erase(storage.begin() + position);
  
  for (int i = position; i < storage.size(); ++i){
    if (storage[i].first < d){
      storage.insert(storage.begin() + i, make_pair(d, index));
      return;
    }
  }

  storage.push_back(make_pair(d, index));
  return;
}

void display(List& l){

  std::cout << "Position: " << l.position << '\n';
  std::cout << "Elements: \n";
  for (auto pair : l.storage){
    std::cout << pair.first;
    std::cout << '\t' << pair.second;
    std::cout << '\n';
  }
  std::cout << std::endl;

}
