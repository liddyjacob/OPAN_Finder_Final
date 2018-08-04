//Tree.hpp
#pragma once


#include <vector>
#include <algorithm>

#include <NTL/ZZ.h>
#include <NTL/RR.h>

using std::vector;
using std::max;

using NTL::ZZ;
using NTL::RR;
using NTL::NextPrime;
using NTL::conv;

struct Node{

  Node(ZZ _p, Node* _parent)
  {

    tried = false;
    success = false;
    parent = _parent;
    if (parent != NULL)
      parent->children.push_back(this);
    p = _p; 

  }

  ZZ p;
  ZZ max;
  bool tried;
  bool success;

  Node* parent;
  std::vector<Node*> children;
};

struct Tree{

  Tree(){
    root = new Node(ZZ(1), NULL);
    height = 0;
    curr = root;
  }

  int height;
  Node* root;
  Node* curr;

};
void grow(Tree& t, RR bound);
void strip_helper(Node* n, vector<ZZ>& primes);

vector<ZZ> strip_primes(Tree& t);

void replace(Tree& t, ZZ rep);

bool fail(vector<ZZ>& primes, Tree& t);
void display(Tree& t);
ZZ backup(vector<ZZ>& primes, Tree& t);
void success(vector<ZZ>& primes, Tree& t);

ZZ find_s(vector<ZZ>& primes, Tree& t);
void replace_next(vector<ZZ>& primes, Tree& t);


ZZ max_branch(vector<ZZ>& primes,Tree& t);

void set_max(Node* n, ZZ& p);

ZZ findmax(vector<ZZ>& primes, int factors, vector<Tree>& trees);
