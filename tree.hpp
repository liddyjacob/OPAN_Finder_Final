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

void failline(int line);


void grow(Tree& t, RR bound, bool verbose);
void strip_helper(Node* n, vector<ZZ>& primes);

vector<ZZ> strip_primes(Tree& t);

void replace(Tree& t, ZZ rep);

bool fail(Tree& t, bool verbose);
void display(Tree& t);
ZZ backup(Tree& t);
void success(Tree& t);

ZZ max_branch(vector<ZZ>& primes,Tree& t);

void set_max(Node* n, ZZ& p);

ZZ findmax(vector<ZZ>& primes, vector<Tree>& trees);
