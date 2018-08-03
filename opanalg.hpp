//This file will house the opn algorithm.
#pragma once
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <vector>
#include <string>
#include <iostream> 
#include <fstream>
#include <omp.h>
#include <fstream>
#include <primecount.hpp>
#include "tools.hpp"
#include "tree.hpp"
#include "expalg.hpp"
using std::vector;
using std::string;

using NTL::ZZ;
using NTL::RR;
using NTL::pow;

//#define MAX 0


/*
 *
 * Stats:
 *   Helps display and collect information
 *   about the results of the algorithm
 *  
 *   Core primes are all but the last prime of a prime sequence
 *   tail prime is the last prime of a prime sequence
 */
struct Stats{

  Stats(string fn = "-", ZZ updatefreq = ZZ(10007))
    : fn(fn), out(std::cout.rdbuf())// 100003
  {
    freq = updatefreq; 
    number_found = 0;
    std::streambuf* buf = std::cout.rdbuf();
    if (fn != "-") { 
      of.open(fn.c_str(), std::ios::out);
      buf = of.rdbuf();
    }

    out.rdbuf(buf);
  }

  ~Stats(){
   freq = 0; 
   if (fn != "-") { 
     of.close();
   }

  }


  std::ofstream of;
  string fn;
  ZZ number_found;
  ZZ freq;
  ZZ product;

  std::ostream out;

  vector<ZZ> prev_core;
  vector<ZZ> init_core;
  vector<vector<ZZ> > prev_exps;
  ZZ init_tail;
  ZZ prev_tail;
};


void efficiency(Stats& s, vector<ZZ>& primes, vector<vector<ZZ> >& exp_seqs);
void Reset(Stats& s);
void Update(Stats& s, Tree& t);
void Show(Stats& s, Tree& t);
void dump_primes(Stats& s);

void Write(Stats& s, vector<ZZ>& primes, vector<vector<ZZ> >& exp_sets);
//void modify(Stats& s, vector<ZZ>& primes, vector<vector<ZZ> >& exp_seqs);

bool cap_check(vector<ZZ>& primes, 
               vector<Tree>& factor_trees, 
               int& factors);

bool exp_find(vector<ZZ>& primes, 
              vector<vector<ZZ> >& exp_seqs);

void OPAN(int d, string fname);
void printline(int line);


//bool isWeird(vector<ZZ>& primes, vector<ZZ>& exps);
//void check_weird(vector<ZZ>& primes,vector<vector<ZZ> >& exp_seqs);
//void check_weirds(ZZ& start, vector<ZZ>& primes,vector<vector<ZZ> >& exp_seqs);
