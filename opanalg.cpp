/* This file contains the main algorithm and its components */
#include <iostream>
#include <limits>
#include <vector>
  using std::string;
  using std::vector;
#include "opanalg.hpp"
//#include "subset_sum.hpp"
#include <primecount.hpp>
#include <algorithm>


void OPAN(int d, string fname){

  std::vector<Tree> record_forest;
  Stats stats(fname);

  for (int factors = 3; factors <= d; ++factors){
    
    record_forest.push_back(Tree());
    Tree& record_tree = record_forest.back();

    /* This is where the algorithm starts in the paper */
    vector<ZZ> primes;
    bool running = true;
 
    vector<vector<ZZ> > exp_seqs; /* Exponent sequences for the prime sequence
                                     (The prime sequence is kept track of by the tree) */  
    while (running){

      if (primes.size() != factors){
        if (b_1(primes) >= RR(2)){
          primes.pop_back();
          ZZ q = min_deficient(primes); /* find minimum prime so that Y(primes + r, {1,1,1,1...}) is deficient */
          primes.push_back(q);
        }
        
        ZZ s = find_s(primes, record_tree);
        primes.push_back(s);
        continue;
      }

      if (b_inf(primes) <= RR(2)){
        running = fail(primes, record_tree);
        continue;
      }

      if (exp_find(primes, exp_seqs)){
        
        /* Issues here */
       
        //Applying Efficiency Theorem...
        Write(stats, primes, exp_seqs); 
        efficiency(stats, primes, exp_seqs);
        Write(stats, primes, exp_seqs);
       
        replace(record_tree, primes.back()); 
        /* ^ Result of the efficiency theorem ^ */
        success(primes, record_tree); /* Store in record tree for later ref*/
        //Write(stats, primes, exp_seqs);

        continue;
      
        /* End issues section */
      }
      
      ZZ s = backup(primes, record_tree);
      primes = strip_primes(record_tree);

      if (primes.size() == factors - 1){ /* Follows from the continuity of OPANS theorem */
        running = fail(primes, record_tree);
        continue;
      }

      /* Find cap_d(P) */
      ZZ cap_d = findmax(primes, factors, record_forest);
      
      /* Check theorem concering capd(P) */
      if (s <= cap_d) {
        grow(record_tree, conv<RR>(s) + RR(.5));        
        primes.push_back(NextPrime(s + ZZ(1)));
        continue;
      }

    /* Theorem tells us that we should modify the primes and cannot look here anymore */
      running = fail(primes, record_tree);
      }
    
    /* Issues in set_max of tree.cpp. Need to Properly find the cap */
    Write(stats, primes, exp_seqs);
    std::cerr << "> Number found with " << factors << " prime factors: " 
              << stats.number_found << '\n';
    //stats.number_found = 0;
    Reset(stats);   
  }
  return;   
}

bool exp_find(vector<ZZ>& primes, vector<vector<ZZ> >& exp_seqs){
  /* First, test all exponents from previous primes,
   * taking advantage of the Efficiency Lemma */

  if (!exp_seqs.empty()){
    int i;

    for (i = 0; i < exp_seqs.size(); ++i){
      vector<ZZ>& exps = exp_seqs[i];
      if (b(primes, exps) < RR(2)){ break; }            
      if (!primitive(primes, exps)) { break; }
    }
    if (i != exp_seqs.size()) {

      /* Efficiency Lemma had a false hypothesis, find exponents manually */
      exp_seqs = expAlg(primes);
      return (!exp_seqs.empty());

    } else {
      /* Efficiecy Lemma had a true hypothesis */
      return true;
    }
  }
  /* There were no given exponents to check. Find exponents manually*/
  exp_seqs = expAlg(primes);
  return (!exp_seqs.empty());
}

ZZ countprimes(ZZ p1, ZZ p2){

  int64_t number = primecount::pi(conv<int64_t>(p2))
                 - primecount::pi(conv<int64_t>(p1));
  return conv<ZZ>(number);
}



void efficiency(Stats& s, vector<ZZ>& primes, vector<vector<ZZ> >& exp_seqs){
  /* Assumes that input primes work with exp_seqs. */
  
  ZZ prime_init = primes.back();

  ZZ increment(2);
  ZZ bad(0);
  bool comedown = false;

  while (increment >= ZZ(2)){//(last_exps == exp_seqs){        
    vector<ZZ> new_primes = primes;
    vector<vector<ZZ> > last_exps = exp_seqs;

    ZZ prev = primes[primes.size() - 1];
    ZZ next = NextPrime(prev + increment);
    new_primes.pop_back();
    new_primes.push_back(next);
 
    if (bad != next){ 
      exp_find(new_primes, exp_seqs);
    } else {
      exp_seqs.clear();
    }

    if (last_exps == exp_seqs){

      if (!comedown){ increment *= ZZ(4); }
      
      primes = new_primes; // Move on to new primes

    } else {
             
      comedown = true;
      increment /= ZZ(4);
      bad = next;
      exp_seqs = last_exps;
    }      
  }

  ZZ prime_final = primes.back();

  s.number_found += countprimes(prime_init - ZZ(1), prime_final + ZZ(1)) * exp_seqs.size();

}


void Reset(Stats& s){

  //std::cout << "Reset\n";

  s.number_found = 0;
  s.product = 0;
  s.prev_core.clear();
  s.init_core.clear();
  s.prev_exps.clear();
  s.init_tail == 0;
  s.prev_tail == 0;
}

void Update(Stats& s, Tree& t){

  //std::cout << "Update\n";

}

void Show(Stats& s, Tree& t){
  vector<ZZ> primes = strip_primes(t);
  std::cout << "P = {";
  printvectos(primes, std::cout);
  std::cout << "}\n";
}

void dump_primes(Stats& s){
  printvectos(s.prev_core, s.out);//s.out);
}

void Write(Stats& s, vector<ZZ>& primes, vector<vector<ZZ> >& exp_sets){
 
  
  //std::cout << "Write\n";

  bool echange = false;
  vector<ZZ> core_primes = primes;
  //std::cout << "Prime length in write: " << primes.size() << '\n'; 
  if (primes.size() != 0){ core_primes.pop_back(); }

  if (exp_sets != s.prev_exps){ 
    echange = true;
  } else if (core_primes != s.prev_core){
    echange = true; /* Dump primes in form p1^e1 p2^e2 p3^e3 p4^e4 q^e5 for e's : and q's: */ 
  }

  
 
  /* If there was a change in exponents, 
   * we need to display all prime sequences
   * with old exponents */

  if (echange) {
    /* Just print one prime */
    if (s.prev_tail == s.init_tail){
      printwithprime(s.prev_core, s.prev_tail, s.out);
      printexponents(s.prev_exps, s.out);
      s.out << '\n';
    } else {
      printwithgeneric(s.prev_core, s.out);
      s.out << "q: " << s.init_tail << " -> " << s.prev_tail << '\n';
      printexponents(s.prev_exps, s.out);
      s.out << '\n';
    }

    if (primes.size() != 0){s.init_tail = primes.back();}

    s.init_core = core_primes;

  }

  if (primes.size() != 0){s.prev_tail = primes.back();}
  s.prev_core = core_primes;
  s.prev_exps = exp_sets;
  
}


