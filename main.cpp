
#include <iostream>
#include <string>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "opanalg.hpp"
using std::string;

bool parse(int& divisors,bool& verbose, bool& steps, string& fname, int argc, char** argv);


int main(int argc, char** argv){
  int divisors;
  bool verbose;
  bool steps;
  string fname;

  if (parse(divisors, verbose, steps, fname, argc, argv)){
    OPAN(divisors, verbose, steps, fname);
    return 0;
  }

  return -1;
}

bool
parse(int& divisors,bool& verbose, bool& steps, string& fname, int argc, char** argv){

  verbose = false;
  steps = false;
  divisors = -1;
  fname = '-';

  if (argc == 1){ 
    std::cerr << "For options, type -help\n";
  }

  for (int i = 1; i < argc; ++i){
    if (strcmp(argv[i], "-help") == 0 ){
      std::cout << "Options:\n\n"
                << "\t-n <Number of Divisors>: Specify number of divisors\n"
                << "\t-v : Verbose\n"
                << "\t-s : Display each step, then wait for user input\n"
                << "\t-f <filename> : Save all integers found to filename\n";
      return false;
    } else 
    if (strcmp(argv[i], "-n") == 0){
      ++i;
      if (i >= argc){
        std::cerr << "Type a positive integer after -n\n";
        return false;
      }
      divisors = atoi(argv[i]);
    } else 
    if (strcmp(argv[i], "-s") == 0){
      steps = true;
    } else 
    if (strcmp(argv[i], "-v") == 0){
      verbose = true;
    } else
    if (strcmp(argv[i], "-f") == 0){
      ++i;
      if (i >= argc){
        std::cerr << "Type a filename after -f\n";
        return false;
      }
      fname = string(argv[i]);

    } else {
      std::cerr << argv[i] << ": Invalid Argument\n";
      return false;
    }

  }

  if (divisors <= -1){ 
    std::cerr << "Please specify the number of prime divisors\n";
    return false;
  }
  return true;

}
