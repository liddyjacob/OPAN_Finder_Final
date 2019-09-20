#include "subset_sum.hpp"
#include "tools.hpp"

void determine(const vector<ZZ>& set, size_t si, ZZ sum, BoolMatrix& bm){
  //Need a base case.i
  //
  /*std::cout << "INFO: \n" 
            << "\tset index...\t" << si << '\n'
            << "\t      sum...\t" << sum << '\n';
  */
  if (sum == 0) { bm[si][sum] = true; return; }
  if (si == 0) { 
    //std::cout << "Checking if sum is element; si = " << si << "\n";
    bm[si][sum] = (sum == set[si]);     
    //std::cout << "\tSums are checked\n"; 
    return; 
  }

  //std::cout << "Eliminating high elements, si = " << si << "\n";
  size_t t = si;
  while (sum < set[t]){
    if (t == 0){ break;}
    t--;
  }

  if (t == 0) {
    bool b = (sum == set[t]);
    for (t; t <= si; ++t){
      bm[t][sum] = b;
    }
    return;
  }
  //std::cout << "\tDone: t = " << t << "\n";
 
  //In most cases, we will be using as many elements as possible,
  if (!bm.defined(t - 1, sum - set[t])){ 
    //std::cout << "Previous case was not defined, si = " << si << "\n";
    determine(set, t - 1, sum - set[t], bm);
    //std::cout << "Result from defined: " << bm[t - 1][ sum - set[t]] << '\n';
  }

  if (bm[t - 1][sum - set[t]]){ bm[si][sum] = true; return;}
    //std::cout << "Previous case was false, si = " << si << '\n';
  if (!bm.defined(t - 1, sum)){
    //std::cout << "Latter case was not defined, si = " << si << "\n";
    determine(set, t - 1, sum, bm);
  }
  //std::cout << "Latter case now defined, si = " << si << "\n";

  bm[si][sum] = bm[t - 1][sum];
  return;
}
bool isSubsetSum2(vector<ZZ>& set, ZZ sum){

  BoolMatrix bm(set.size());
  determine(set, set.size() - 1, sum, bm); 
  return bm[set.size() - 1][sum];
}


