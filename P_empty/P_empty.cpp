#include <iostream>
#include <random>
#include <string>
#include <set>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <map>
#include <unordered_map>
#include "Newton_method.h"
#include "DP_Prob_Empty.h"
#include <cmath>
using namespace std;


#include "mpreal.h"
using mpfr::mpreal;
// Required precision of computations in decimal digits
const int digits = 50;  // Play with it to check different precisions

int main (int argc, char* argv[]){
    mpreal::set_default_prec(mpfr::digits2bits(digits));
    srand(time(nullptr));
    string filename="";
    int L = 0;
    int k = 0;
    int start_pos = 0;
    double epsilon = 0.000001;
    int num_replicates = 0;
    double r = 0;
    int type = 0;

    vector<string> args(argv + 1, argv + argc);
    for (auto i = args.begin(); i != args.end(); ++i) {
        if (*i == "-h" || *i == "--help") {
            cout << "Syntax: -i [input-file] -l <length>[100] -s <start position>[5]  -k <kmer-size>[33] -e <presion>[0.001] -c <replicates>[30] -r <mutation rate> 0.1 -t <type> 0" << endl;
            return 0;
        } else if (*i == "-l") {
            L = stoull(*++i);
        } else if (*i == "-k") {
            k = stoi(*++i);
        } else if (*i == "-i") {
            filename = (*++i);
        } else if (*i == "-s") {
            start_pos = stoi(*++i);
        } else if (*i == "-e") {
            epsilon = stod(*++i);
        } else if (*i == "-c") {
            num_replicates = stoi(*++i);
        } else if (*i == "-r") {
            r = stod(*++i);
        } else if (*i == "-t") {
            type = stoi(*++i);
        }
    }
    // SimParams sim(k, L, r);
    // std::string input_String;
    // //input_String = readSubstring(filename, start_pos ,L + k - 1);

    // if (filename==""){
    //     if(type == 0){
    //         input_String = generateRandomString(L + k - 1);
    //     }else{
    //         input_String = generateStressTestString(L + k - 1, k, type);
    //     }
    // }else{
    //     input_String = readSubstring(filename, start_pos , L + k - 1);
    // }
    //cout << "input" << input_String << endl;
    // ResultAggregator res;
    // mpreal r_mp = sim.r;
    // mpreal one_minus_r = 1 - r_mp;        
    // mpreal rr = pow(one_minus_r, k);

    // set<string> kmer_set = kspectrum_update(input_String, k);
    // pair<double,pair<double,double>> est_error = cal_fun_f(sim, kmer_set, input_String, occurrence);
    // double estimator = est_error.first;
    // double theoretical_err_lower = est_error.second.first;
    // double theoretical_err_upper = est_error.second.second;
    // E_upper = estimator + theoretical_err_upper;
    // E_lower = estimator + theoretical_err_lower;
    // std::unordered_map<int, int> coeffs_Map; // a * x^b, b->a
    // int max_copy = 0;
    // for (const auto& pair : occurrence){
    //     int copy = pair.second.size();
    //     if (copy > max_copy){
    //         max_copy = copy;
    //     }
    //     coeffs_Map[copy]++;
    // }



    // mpreal q_high = Newton_Method (coeffs_Map, kmer_set, k, epsilon, E_lower.toDouble());
    // mpreal q_low = Newton_Method (coeffs_Map, kmer_set, k, epsilon, E_upper.toDouble());
    // mpreal boundary = Newton_Method (coeffs_Map, kmer_set, k, epsilon, 1.0);
    // pair<mpreal,mpreal> g1 = min_g1(coeffs_Map, k,boundary);
    // mpreal g1_min = g1.first;
    // mpreal min_q1 = g1.second;
    // cout << "min q1: " << min_q1 << endl;
    mpreal Pr_empty_upper = Dynamic_Programming_Pr_emtpy(L, k, r);

    // mpreal lower_bound_r, upper_bound_r;
    // lower_bound_r = (1.0 - pow(1.0 - q_low, 1.0/k))* (1- Pr_empty_upper) + g1_min * E_upper * Pr_empty_upper;
    // upper_bound_r = (1.0 - pow(1.0 - q_high, 1.0/k))- g1_min * (E_upper) + Pr_empty_upper;
    // cout << "g1 min: " << g1_min << " U_E: " << E_upper << endl;
    // std::cout.precision(10);
    // cout << lower_bound_r << "," << upper_bound_r << endl;
    cout << Pr_empty_upper << endl;
    // for(int i=0; i<num_replicates; i++){
    //     std::string mutatedString = generateMutatedString(input_String, r);
    //     int I_obs = intersect_size(kmer_set, mutatedString, k);
    //     res.isizes.push_back(I_obs);
    //     mpreal q_hat = Newton_Method (coeffs_Map, kmer_set, k, epsilon, I_obs);
    //     mpreal r_hat = 1.0 - pow(1.0 - q_hat, 1.0/k);
    //     cout << r_hat << endl; 
    // }
    
    
}