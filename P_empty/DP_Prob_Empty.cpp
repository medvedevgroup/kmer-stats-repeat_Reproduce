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

mpreal Dynamic_Programming_Pr_emtpy(int n, int k, mpreal r){ //number of kmer
    vector<mpreal> DP_array(n,mpreal(0));
    mpreal q = 1 - pow(1-r,k);
    DP_array[0]= q;
    for (int i = 1; i<k ; i++){
        DP_array[i] = ((1-pow(1-r,i))*q) + (pow(1-r,i)*(1-pow(1-r,k-i)));
    }
    for (int i = k; i < n;i++){
        mpreal Pr = 0;
        for (int j = 0; j<k; j++){
            Pr += DP_array[i-j-1] * r * pow(1-r,j);
        }
        DP_array[i] = Pr;
    }
    return DP_array[n-1];
}