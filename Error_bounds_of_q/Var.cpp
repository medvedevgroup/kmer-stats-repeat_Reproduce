#include <iostream>
#include <random>
#include <string>
#include <set>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <map>
#include <unordered_map>
using namespace std;


#include "mpreal.h"
using mpfr::mpreal;
// Required precision of computations in decimal digits
const int digits = 50;  // Play with it to check different precisions

int num_sep(const set<int>& set_tau, const set<int>& set_upsilon, int k){
    set<int> mergedSet(set_tau.begin(), set_tau.end());
    mergedSet.insert(set_upsilon.begin(), set_upsilon.end());
    int num_cluster = 0;
    int cur_pos_tau = -k;
    for (const auto& occ: mergedSet){
        if (occ >= cur_pos_tau +k){
            num_cluster++;
            cur_pos_tau = occ;
        }
    }
    return num_cluster;
}

mpreal Var(int k, mpreal rr, const set<string>& ksp_set, const unordered_map<string, set<int>>& occurrence, const unordered_map<string, mpreal>& kmer_lower, mpreal E_upper){
    mpreal Var_upper;
    for (const auto& tau: ksp_set){
        for (const auto& upsilon: ksp_set){
            if (tau != upsilon){
                Var_upper += pow((1-rr),num_sep(occurrence.at(tau),occurrence.at(upsilon), k));
            }
        }
    }
    if ((ksp_set.size() - E_upper)>= 0.5){
        Var_upper += ((ksp_set.size() - E_upper) - (ksp_set.size() - E_upper)*(ksp_set.size() - E_upper));
    }
    else{
        Var_upper += 1.0/4.0;
    }
    
    return Var_upper;
}