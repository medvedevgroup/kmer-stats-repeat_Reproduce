#ifndef NEWTON_METHOD_H
#define NEWTON_METHOD_H
#include "mpreal.h"
#include <string>
#include <set>
using mpfr::mpreal;
using namespace std;
mpreal Polynomial (const unordered_map<int,int>& coeffs_Map, mpreal q);
mpreal Newton_Method (const unordered_map<int, int>& coeffs_Map, const set<string>& ksp_set, int k, double epsilon, double I_obs);
//pair<mpreal,mpreal> Taylor_lower_bound_2nd(const unordered_map<int, int>& coeffs_Map, mpreal Var_upper);
// mpreal Taylor_lower_bound_2nd(const unordered_map<int, int>& coeffs_Map, mpreal q_high, mpreal q_low, mpreal Var_upper);
// //mpreal Taylor_lower_bound_3rd(int a1, int d);
// pair<mpreal,mpreal> Taylor_lower_bound_3rd(const unordered_map<int, int>& coeffs_Map, mpreal Var_upper);
// pair<pair<mpreal,mpreal>,pair<mpreal,mpreal>> min_gamma(const unordered_map<int, int>& coeffs_Map, int k);
pair<mpreal,mpreal> min_g1(const unordered_map<int, int>& coeffs_Map, int k, mpreal boundary);
// pair<mpreal,mpreal> max_g2(const unordered_map<int, int>& coeffs_Map, int k, mpreal left, mpreal right);
// pair<mpreal,mpreal> max_g3(const unordered_map<int, int>& coeffs_Map, int k, mpreal boundary);
#endif