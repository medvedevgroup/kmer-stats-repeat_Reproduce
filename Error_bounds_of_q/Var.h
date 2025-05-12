#ifndef VAR_H
#define VAR_H
#include "mpreal.h"
#include <string>
#include <set>
#include <unordered_map>
using mpfr::mpreal;
using namespace std;

mpreal Var(int k, mpreal rr, const set<string>& ksp_set, const unordered_map<string, set<int>>& occurrence, const unordered_map<string, mpreal>& kmer_lower, mpreal E_upper);
#endif