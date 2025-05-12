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



mpreal Taylor_lower_bound_2nd(const unordered_map<int, int>& coeffs_Map, mpreal q_high, mpreal q_low, mpreal Var_upper){
    mpreal::set_default_prec(mpfr::digits2bits(digits));
    mpreal lower_bound = 0.5 * (Polynomial_nth_derivative(coeffs_Map, q_high,2) / pow(Polynomial_nth_derivative(coeffs_Map, q_low,1),3)) * Var_upper;
    return lower_bound;
}

mpreal Taylor_lower_bound_3rd(const unordered_map<int, int>& coeffs_Map){
    mpreal::set_default_prec(mpfr::digits2bits(digits));
    mpreal max_lower_bound=0;
    mpreal lower_bound=0;
    mpreal max_q = 0;
    for (int i = 0; i <= 1000; i++){
        mpreal q = mpreal(i) / 1000;
        mpreal F1 = Polynomial_nth_derivative(coeffs_Map, q,1);
        lower_bound = abs((Polynomial_nth_derivative(coeffs_Map, q,3) * F1) - (mpreal(3) * pow(Polynomial_nth_derivative(coeffs_Map, q,2),2)))
        / (mpreal(6) * pow(F1,5));
        // cout << q << ": " << lower_bound << endl;
        if (lower_bound > max_lower_bound){
            max_lower_bound = lower_bound;
            max_q = q;
        }
    }
    return lower_bound;
}

// pair<pair<mpreal,mpreal>,pair<mpreal,mpreal>> min_gamma(const unordered_map<int, int>& coeffs_Map, int k){
//     mpreal::set_default_prec(mpfr::digits2bits(digits));
//     mpreal gamma_min = 0;
//     mpreal min_r = 0;
//     mpreal gamma_max = 0;
//     mpreal max_r = 0;
//     for (int i = 0; i <= 990; i++){
//         mpreal r = mpreal(i) / 1000;
//         mpreal F1 = Polynomial_nth_derivative(coeffs_Map, r,1);
//         mpreal F2 = Polynomial_nth_derivative(coeffs_Map, r,2);
//         mpreal one_over_k =  mpreal(1.0)/mpreal(k);
//         mpreal H = (one_over_k * pow(1.0-r,one_over_k-2) * mpreal(0.5)) * ((1.0-one_over_k)*F1 - (1-r)*F2) / (pow(F1,3));
//         if (H < gamma_min){
//             gamma_min = H;
//             min_r = r;
//         }
//         if (H > gamma_max){
//             gamma_max = H;
//             max_r = r;
//         }
//     }
//     return {{gamma_min, min_r},{gamma_max, max_r}};
// }

// pair<mpreal,mpreal> min_g1(const unordered_map<int, int>& coeffs_Map, int k, mpreal boundary){
//     mpreal::set_default_prec(mpfr::digits2bits(digits));
//     mpreal g1_min = 0;
//     mpreal min_q = 0;

//     for (int i = 0; i <= 1000; i++){
//         mpreal q = mpreal(i)/mpreal(1000) * boundary;
//         mpreal F1 = -1 * Polynomial_nth_derivative(coeffs_Map, q,1);
//         mpreal one_over_k =  mpreal(1.0)/mpreal(k);
//         mpreal H = one_over_k * pow(1.0-q,one_over_k-1)  / F1;
//         if (H < g1_min){
//             g1_min = H;
//             min_q = q;
//         }
//     }
//     return {g1_min, min_q};
// }

// pair<mpreal,mpreal> max_g2(const unordered_map<int, int>& coeffs_Map, int k, mpreal left, mpreal right){
//     mpreal::set_default_prec(mpfr::digits2bits(digits));
//     mpreal g2_max = 0;
//     mpreal max_q = 0;

//     for (int i = 0; i <= 1000; i++){
//         mpreal q = mpreal(i)/mpreal(1000) * (right -left) + left;
//         mpreal F1 = Polynomial_nth_derivative(coeffs_Map, q,1);
//         mpreal F2 = Polynomial_nth_derivative(coeffs_Map, q,2);
//         mpreal one_over_k =  mpreal(1.0)/mpreal(k);
//         mpreal H = abs((one_over_k * pow(1.0-q,one_over_k-2) * mpreal(0.5)) * ((1.0-one_over_k)*F1 - (1-q)*F2) / (pow(F1,3)));
//         if (H > g2_max){
//             g2_max = H;
//             max_q = q;
//         }
//     }
//     return {g2_max,max_q };
// }

// pair<mpreal,mpreal> max_g3(const unordered_map<int, int>& coeffs_Map, int k, mpreal boundary){
//     mpreal::set_default_prec(mpfr::digits2bits(digits));
//     mpreal g3_max = 0;
//     mpreal max_q = 0;

//     for (int i = 0; i <= 1000; i++){
//         mpreal q = mpreal(i)/mpreal(1000) * boundary;
//         mpreal F1 = Polynomial_nth_derivative(coeffs_Map, q,1);
//         mpreal F2 = Polynomial_nth_derivative(coeffs_Map, q,2);
//         mpreal F3 = Polynomial_nth_derivative(coeffs_Map, q,2);
//         mpreal one_over_k =  mpreal(1.0)/mpreal(k);
//         mpreal part1 = (1-one_over_k)*(2-one_over_k)*pow(1.0/(-F1),3);
//         mpreal part2 = 3 * (1-one_over_k) * (1-q) * (F2) * pow(F1,-4);
//         mpreal part3 = pow((1-q),2) * (F3*F1 - 3* F2*F2) * pow(-F1,-5);
//         mpreal H = abs((1.0 / 6.0)*pow(1.0-q,one_over_k-3)*one_over_k * (part1+part2+part3));
//         if (H > g3_max){
//             g3_max = H;
//             max_q = q;
//         }
//     }


//     return {g3_max,max_q };
// }