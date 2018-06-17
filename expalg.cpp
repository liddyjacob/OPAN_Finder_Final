// The exponent algorithms.
#include "expalg.hpp"
#include <vector>
#include <NTL/ZZ.h>
#include <algorithm> // std::fill

#include "tools.hpp"
#include "list.hpp" // Doubly linked list, sorted, of integers.
using std::vector; // vector<>
using NTL::ZZ; // ZZ

/*
  expAlg: Returns the exponents found when looking for all
primitive abundant numbers whose primes compose P(primes)
*/
vector<vector<ZZ> > expAlg(vector<ZZ>& primes){

  int divisors = primes.size();
  vector<ZZ> expos(divisors);
  std::fill (expos.begin(), expos.end(), ZZ(1));
  
  vector<vector<ZZ> > exp_seqs;

  if (b_inf(primes) < RR(2)){ return exp_seqs; }

  initiate(primes, expos, exp_seqs);

  return exp_seqs;
}

void initiate(vector<ZZ>& primes, 
              vector<ZZ>& expos, 
              vector<vector<ZZ> >& exp_seqs){

  #ifdef DEBUG
    std::cout << "In initiate with expos ";
    printvect(expos);
    std::cout << std::endl;
  #endif

  if (b(primes, expos) >= RR(2)){
    if (b(primes, expos) / b(primes.back()) < RR(2)){
      exp_seqs.push_back(expos);
    }
    return;
  }
  
  List G(primes); // Initialize a doubly linked list as 0, 1, ..., p.

  orderDel(primes, expos, exp_seqs, G);

  return;
}

bool orderDel(vector<ZZ>& primes,
              vector<ZZ>& expos,
              vector<vector<ZZ> >& exp_seqs,
              List& G){

  #ifdef DEBUG
    std::cout << "In orderDel with expos ";
    printvect(expos);
    std::cout << std::endl;
    display(G);
  #endif

  bool found = false;
  int i = 0;
  RR precalc_b = b(primes, expos);

  while (G.at() < G.size()){

    if (precalc_b * G.delta() >= RR(2)){
      //Set new exponents:
      vector<ZZ> expos_new = expos;
        expos_new[G.index()] += ZZ(1);

      if (primitive(primes, expos_new, G)){
        exp_seqs.push_back(expos_new);
        found = true;
      }

      ++G;
    } else {
      break;
    }
  }

  if (G.at() != G.size()){
    if (expIncAlg(primes, expos, exp_seqs, G)){ 
      found = true;
    }
  }

  return found;
}
      

bool expIncAlg(vector<ZZ>& primes,
               vector<ZZ>& expos,
               vector<vector<ZZ> >& exp_seqs,
               List& G){

  #ifdef DEBUG
    std::cout << "In expIncAlg with expos ";
    printvect(expos);
    std::cout << std::endl;
    display(G); std::cout << '\n';
  #endif
  bool found = false;

  while (G.at() < G.size()){

    
    if (mb(primes, expos, G) > RR(2)){

      vector<ZZ> expos_new = expos;
        expos_new[G.index()] += ZZ(1);        

      List Gnew(G);
      //Does NOT modify iter - modifies contents of G.

      Gnew.resort(primes[G[G.at()]], expos_new[G[G.at()]]);
      if (orderDel(primes, expos_new, exp_seqs, Gnew)){
        found = true;
      }

      ++G;

    } else { break; }
  }
  return found;
}
