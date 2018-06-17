//Tools.
#pragma once

#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <iostream>
#include <vector>
 
#include "list.hpp"

using NTL::RR;
using NTL::ZZ;
using std::vector;

using NTL::conv;
using NTL::NextPrime;
//#define DEBUG
void printwithprime(vector<ZZ>& primes, ZZ& prime, std::ostream& stream);
void printwithgeneric(vector<ZZ>& primes, std::ostream& stream);
void printexponents(vector<vector<ZZ> >& exp_sets, std::ostream& stream);
template<typename T>
void printvectos(vector<T>& vect, std::ostream& stream){
  for (auto v : vect){
    stream << v << ", ";
  } 
}

template<typename T>
void printvect(vector<T> vect){
  for (auto v : vect){
    std::cout << v << ", ";
  }
}

RR Delta(ZZ& p, ZZ& e);
RR Delta(ZZ& p);

RR Delta(ZZ& p, ZZ& e);

RR Delta(ZZ& p);
RR del_neg(ZZ& p, ZZ& e);
RR b_inf(ZZ& p); 
RR b_inf(vector<ZZ>& primes);
RR b(ZZ& p, ZZ& e);
RR b_1(ZZ& p);
RR b_1(vector<ZZ>& prime);
RR b(ZZ& p);
RR b(vector<ZZ>& primes, vector<ZZ>& expos);
bool primitive(vector<ZZ>& primes, vector<ZZ>& expos, List& l);
//REQUIRES AT LEAST TWO PRIMES.
bool primitive(vector<ZZ>& primes, vector<ZZ>& expos);
RR mb(vector<ZZ>& primes,vector<ZZ>& expos, List& l);
ZZ min_deficient(vector<ZZ>& primes);


ZZ product(vector<ZZ>& primes, vector<ZZ>& exps);
ZZ abundance(vector<ZZ>& primes, vector<ZZ>& exps);
vector<ZZ> div_below_a(vector<ZZ>& primes, vector<ZZ>& exps);

