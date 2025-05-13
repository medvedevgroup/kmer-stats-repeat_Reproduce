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


mpreal Polynomial (const unordered_map<int,int>& coeffs_Map, mpreal q){
    mpreal::set_default_prec(mpfr::digits2bits(digits));
    mpreal poly = 0;
    for (const auto& [exp, mul] : coeffs_Map) {
        poly += mul * pow(q , exp);
    }
    return poly;
}

mpreal Polynomial_prime (const unordered_map<int,int>& coeffs_Map, mpreal q){
    mpreal::set_default_prec(mpfr::digits2bits(digits));
    mpreal poly_prime = 0;
    for (const auto& [exp, mul] : coeffs_Map) {
        poly_prime += exp * mul * pow(q , exp -1);
    }
    return poly_prime;
}

mpreal Polynomial_nth_derivative (const unordered_map<int,int>& coeffs_Map, mpreal q, int n){
    mpreal::set_default_prec(mpfr::digits2bits(digits));
    mpreal poly_nth_prime = 0;
    for (const auto& [exp, mul] : coeffs_Map) {
        if (exp >= n){
            mpreal nth_coef = 1;
            for (int i = 0; i < n; i++){
                nth_coef *= (exp-i);
            }
            poly_nth_prime += nth_coef * mul * pow(q , exp - n);
        }
    }
    return poly_nth_prime;
}



mpreal Newton_Method (const unordered_map<int, int>& coeffs_Map, const set<string>& ksp_set, int k, double epsilon, double I_obs){
    mpreal::set_default_prec(mpfr::digits2bits(digits));
    //set<string> ksp_set = kspectrum(seq, k);
    
    // for (const auto& pair : coeffs_Map) {
    //     std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
    // }
    //cout << "a: " << (ksp_set.size() - I_obs) << endl;
    mpreal q_hat = 1 - pow(0.09, k); 
    mpreal f = Polynomial(coeffs_Map, q_hat) - (ksp_set.size() - I_obs);      // P(q) - a
    while (abs(f) > epsilon) {
        f = Polynomial(coeffs_Map, q_hat) - (ksp_set.size() - I_obs);      // P(q) - a
        mpreal fp = Polynomial_prime (coeffs_Map, q_hat);   // P'(q)
        q_hat -= f / fp;
        //cout << "f: " << f << endl;
        //cout << "fp: " << fp << endl;
        //cout << "q_hat: " << q_hat << endl;
    }


    return q_hat;
}


pair<mpreal,mpreal> min_g1(const unordered_map<int, int>& coeffs_Map, int k, mpreal boundary){
    mpreal::set_default_prec(mpfr::digits2bits(digits));
    mpreal g1_min = 0;
    mpreal min_q = 0;

    for (int i = 0; i <= 1000; i++){
        mpreal q = mpreal(i)/mpreal(1000) * boundary;
        mpreal F1 = -1 * Polynomial_nth_derivative(coeffs_Map, q,1);
        mpreal one_over_k =  mpreal(1.0)/mpreal(k);
        mpreal H = one_over_k * pow(1.0-q,one_over_k-1)  / F1;
        if (H < g1_min){
            g1_min = H;
            min_q = q;
        }
    }
    return {g1_min, min_q};
}
