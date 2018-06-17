//Tools.
#include <tools.hpp>

#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <iostream>

#include <vector>

using NTL::RR;
using NTL::ZZ;
using std::vector;
//using NTL::power

using NTL::conv;
using NTL::NextPrime;
//#define DEBUG


RR Delta(ZZ& p, ZZ& e);
RR Delta(ZZ& p);

RR Delta(ZZ& p, ZZ& e){

  RR pr = conv<RR>(p); // convPrec(p, 53)
  long er = conv<long>(e);

  return (power(pr, er + 2.0) - RR(1)) /
/*      -----------------------------------       */
         (pr * (power(pr, er + 1.0) - RR(1)));
}

RR Delta(ZZ& p){
  ZZ e(1);
  return Delta(p, e);
}

RR del_neg(ZZ& p, ZZ& e){

  RR pr = conv<RR>(p);
  long er = conv<long>(e);

  return (pr * (power(pr, er) - RR(1))) /
/*       --------------------------      */
         (power(pr, er + 1.0) - RR(1));
}

#include "list.hpp"

RR b_inf(ZZ& p) {
  RR pr = conv<RR>(p);
  return (pr)/(pr - RR(1));
}

RR b_inf(vector<ZZ>& primes){

  RR product(1);

  for (prime : primes){ product *= b_inf(prime); }

  return product;
}

RR b(ZZ& p, ZZ& e){
  RR pr = conv<RR>(p);
  long er = conv<long>(e);
  return (power(pr, er + 1.0) - RR(1)) /
/*       ----------------------------    */
         ((pr - RR(1)) * power(pr, er));

}

RR b_1(ZZ& p){

  RR pr = conv<RR>(p);
  return (pr + RR(1)) / (pr);
}

RR b_1(vector<ZZ>& prime){  

  RR product(1);
  for (auto& p : prime){product*=b_1(p);}

  return product;
}

RR b(ZZ& p){
  ZZ e(1);
  return b(p, e);
}

RR b(vector<ZZ>& primes, vector<ZZ>& expos){

  #ifdef DEBUG
    std::cout << "In b with: ";
    printvect(primes);
    std::cout << "\nAnd exponents: ";
    printvect(expos);
    std::cout << std::endl;
  #endif

  if (primes.size() != expos.size()){
    std::cerr << "from b(P, E): Exponent and prime sequence sizes are not the same!\n";
    return RR(-1);
  }

  RR product(1);

  for (int i = 0; i < primes.size(); ++i){ product *= b(primes[i], expos[i]); }

  return product;
}
//REQUIRES AT LEAST TWO PRIMES.
bool primitive(vector<ZZ>& primes, vector<ZZ>& expos, List& l){

  #ifdef DEBUG
    std::cout << "In Primitive with: ";
    printvect(primes);
    std::cout << "\nAnd exponents: ";
    printvect(expos);
    display(l);
  #endif

  for (int i = l.size() - 1; i > 0; --i){
    ZZ p = primes[l[i]];
    ZZ e = expos[l[i]];
 
    if (b(primes, expos) * del_neg(p, e) >= RR(2))
      return false;
  }

  #ifdef DEBUG
    std::cout << "Primitive Succeeded\n";
  #endif

  return true;
}


//REQUIRES AT LEAST TWO PRIMES.
bool primitive(vector<ZZ>& primes, vector<ZZ>& expos){

  #ifdef DEBUG
    std::cout << "In Primitive(nolist) with: ";
    printvect(primes);
    std::cout << "\nAnd exponents: ";
    printvect(expos);
  #endif

  for (int i = primes.size() - 1; i >= 0; --i){
    ZZ p = primes[i];
    ZZ e = expos[i];
 
    if (b(primes, expos) * del_neg(p, e) >= RR(2))
      return false;
  }

  #ifdef DEBUG
    std::cout << "Primitive Succeeded\n";
  #endif

  return true;
}

RR mb(vector<ZZ>& primes,vector<ZZ>& expos, List& l){
  #ifdef DEBUG
    std::cout << "In mb with: ";
    printvect(primes);
    std::cout << "\nAnd exponents: ";
    printvect(expos);
    std::cout << std::endl;
  #endif

  RR product(1);

  for (int i = 0; i < l.at(); ++i){
    product *= b(primes[l[i]], expos[l[i]]);
  }
  
  #ifdef DEBUG
    std::cout << "Finite Product: " << product << std::endl;
  #endif

  for (int i = l.at(); i < l.size(); ++i){
    product *= b_inf(primes[l[i]]);
  }

  #ifdef DEBUG
    std::cout << "Total Product: " << product << std::endl;
  #endif
 
  return product;
}

ZZ min_deficient(vector<ZZ>& primes){

  RR b1 = b_1(primes);

  RR bound = b1 / (RR(2) - b1);

  ZZ prime = NextPrime(NTL::CeilToZZ(bound));

  return prime;
}

ZZ product(vector<ZZ>& primes, vector<ZZ>& exps){
  //assert(primes.size() == exps.size());
  ZZ product(1);
  for (int i= 0; i < primes.size(); ++i){
    product *= NTL::power(primes[i], NTL::conv<long>(exps[i]));
  }
  return product;
}

ZZ a(vector<ZZ>& primes, vector<ZZ>& exps){
  ZZ prod(1);
  for (int i = 0; i < primes.size(); ++i){
    ZZ& p = primes[i];
    long e = conv<long>(exps[i]);
    prod *= (power(p, e + 1.0) - 1.0) / (p - 1.0);
  }

  return prod;
}

ZZ abundance(vector<ZZ>& primes, vector<ZZ>& exps){
  return a(primes, exps) - ZZ(2) * product(primes, exps);
}

vector<ZZ> div_below_a(vector<ZZ>& primes, vector<ZZ>& exps){

  vector<ZZ> divisors;
    divisors.push_back(ZZ(1));

  ZZ ab = abundance(primes, exps);

  for (int i = 0; i < primes.size(); ++i){

    ZZ ppower(1);
    size_t oldsize = divisors.size();

    for (ZZ exp(1); exp <= exps[i]; ++exp){

      mul(ppower, ppower, primes[i]);      
      if (ppower > ab){break;}

      for (int di = 0; di < oldsize; ++di){
        ZZ newdivisor = ppower * divisors[di];

        if (newdivisor <= ab){divisors.push_back(newdivisor);}
        else break;
      }
    }

    if (ppower == 1){ break; }
  }
  return divisors;
}


void printwithgeneric(vector<ZZ>& primes, std::ostream& stream){


  if (primes.size() == 0) {return;}

  int i = 0;
  
  for (i = 0; i < primes.size(); ++i){
    stream << primes[i] << "^" << "e" << i << ' ';
  }
  stream << "q^e" << i << '\n';
}

void printwithprime(vector<ZZ>& primes, ZZ& prime, std::ostream& stream){

  if (primes.size() == 0){ return; }

  int i = 0;
  
  for (i = 0; i < primes.size(); ++i){
    stream << primes[i] << "^" << "e" << i << ' ';
  }
  stream << prime << "^e" << i << '\n';

}

void printexponents(vector<vector<ZZ> >& exp_sets, std::ostream& stream){

  if (exp_sets.size() == 0) return;

  stream << "With {";
  for (int i = 0; i < exp_sets[0].size() ; ++i){
    stream << "e" << i << ", ";
  }
  stream << "}\n";

  for (vector<ZZ>& exps : exp_sets){
    stream << "     {";
    for (ZZ& e : exps){
      stream << e << ", ";
    }
    stream << "}\n";
  }
}


