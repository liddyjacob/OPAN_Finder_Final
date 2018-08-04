//Tree.hpp

#include "tree.hpp"

#include <string>
#include <vector>
#include <algorithm>

#include <NTL/ZZ.h>
#include <NTL/RR.h>

using std::vector;
using std::max;
using std::string;

using NTL::ZZ;
using NTL::RR;
using NTL::NextPrime;
using NTL::conv;


void grow(Tree& t, RR bound){

  std::cout << "In grow\n";

  bound = max(bound, conv<RR>(t.curr->p + ZZ(1)));

  // If there are any children: Need to worry about them.
  if (!t.curr->children.empty()){

    RR maxchild = conv<RR>(t.curr->children.back()->p);
    bound = max(maxchild + RR(1), bound);

  }

  ZZ nextprime = NextPrime(NTL::CeilToZZ(bound));

  Node* extend = new Node(nextprime, t.curr);
  
  t.curr = extend;
}

void strip_helper(Node* n, vector<ZZ>& primes){
  if (n->p == 1) { return; }
  strip_helper(n->parent, primes);
  primes.push_back(n->p);
}

vector<ZZ> strip_primes(Tree& t){
  vector<ZZ> primes;
  strip_helper(t.curr, primes);
  return primes;
}

void replace(Tree& t, ZZ rep){
  t.curr->p = rep;
}
bool fail(vector<ZZ>& primes, Tree& t){

  while (!t.curr->success){

    Node* failed = t.curr;
    t.curr = t.curr->parent;
    delete failed;
    t.curr->children.pop_back();
    primes.pop_back();

    if (t.curr->p == ZZ(1) && (!t.curr->success))
      
      return false;
  }

  t.curr->tried = true;

  while (t.curr->tried){
    t.curr = t.curr->parent;
    primes.pop_back();
  }
  return true;
}

void displaynode(Node* n, int depth){
  for (int i = 0; i < depth; ++i){ std::cout << "---";} 
  std::cout << '+';
  std::cout << n->p << '(' << n->max << ')';
  return;
}

void displayhelper(Node* n, int depth = 0){
  displaynode(n, depth);
  for (int i = 0; i < n->children.size(); ++i){
    std::cout << std::endl;
    for (int i = 0; i < depth; ++i){std::cout << "   ";}
    std::cout << '|'; 
    displayhelper(n->children[i], depth + 1);
  }
}

void display(Tree& t){
  Node* curr = t.root;
  displayhelper(curr);
  std::cout << '\n';
}

ZZ backup(vector<ZZ>& primes, Tree& t){

  ZZ lastprime(-1);
  while (!t.curr->success){

    lastprime = t.curr->p;

    Node* failed = t.curr;
    t.curr = t.curr->parent;
    delete failed;
    t.curr->children.pop_back();
    primes.pop_back();


    if (t.curr->p == ZZ(1) && (!t.curr->success))
      return ZZ(1);
  }

  //primes = strip_primes(t);
  return lastprime;
}

void success(vector<ZZ>& primes, Tree& t){
  Node* node = t.curr;
  node = node->parent;

  if (node->children.size() == 2){
    delete node->children[0];
    node->children.erase(node->children.begin());
  }

  while (node->p != ZZ(1)){
    node->success = true;
    node = node->parent; 

    // A shortcut:
    if (node->success == true) {break;}
    //primes.pop_back(); 
  }

  t.curr = t.curr->parent;
  primes = strip_primes(t);
}


ZZ findmax(vector<ZZ>& primes, int factors, vector<Tree>& Trees){

  int divisors = Trees.size() + 2;
  ZZ max(0);  

  for (int d = divisors; d > primes.size(); d--){
    Tree& t = Trees[d - 3]; // 
    max = std::max(max_branch(primes, t), max);
  }

  return max;
}

ZZ max_helper(Node* n, vector<ZZ>& primes, int i = 0){
  if (i == primes.size()){
    return n->max;
  }

  for(int j = 0; j < n->children.size(); ++j){
    //Replace with binary search.
    if (n->children[j]->p == primes[i]){
      return max_helper(n->children[j], primes, i + 1);
    }
  }

  return ZZ(0); // Incomplete;
}

ZZ max_branch(vector<ZZ>& primes, Tree& t){

  Node* curr = t.root;
  int i = 0;

  return max_helper(curr, primes);
}



void set_max(Node* n, ZZ& p){

  n->p = p;
  Node* curr = n;  

  while (curr->p != ZZ(1)){
    if (p > curr->max){
      curr->max = p;
      curr = curr->parent;
    } else {return;}
  }
  curr->max = p;
}
ZZ find_s(vector<ZZ>& primes, Tree& t){
  grow(t, RR(3));
  ZZ s = t.curr->p;
  return s;
}

void replace_next(vector<ZZ>& primes, Tree& t){ return ;}


