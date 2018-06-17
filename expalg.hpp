// The exponent algorithms.
#pragma once

#include <vector>
#include <NTL/ZZ.h>
#include <algorithm> // std::fill

//#include "tools.hpp"
#include "list.hpp" // Doubly linked list, sorted, of integers.
using std::vector; // vector<>
using NTL::ZZ; // ZZ

/*
  expAlg: Returns the exponents found when looking for all
primitive abundant numbers whose primes compose P(primes)
*/
vector<vector<ZZ> > expAlg(vector<ZZ>& primes);
void initiate(vector<ZZ>& primes, 
              vector<ZZ>& expos, 
              vector<vector<ZZ> >& exp_seqs);
bool orderDel(vector<ZZ>& primes,
              vector<ZZ>& expos,
              vector<vector<ZZ> >& exp_seqs,
              List& G);
bool expIncAlg(vector<ZZ>& primes,
               vector<ZZ>& expos,
               vector<vector<ZZ> >& exp_seqs,
               List& G);

vector<vector<ZZ> > expAlg(vector<ZZ>& primes);

