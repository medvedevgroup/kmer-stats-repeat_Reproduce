#include <iostream>
#include <random>
#include <string>
#include <set>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <map>
#include <unordered_map>
#include "MurmurHash3.h"
using namespace std;



set<string> kspectrum_sketch(string s, int k, double theta, uint32_t seed){
    //uint32_t seed = 2016;
    // uint64_t largest_value = 0xFFFFFFFFFFFFFFFF;
	// uint64_t threshold = std::round((long double)(largest_value)/(long double)(scaled));
    set<string> kmers;
    for(int i=0; i<s.length()-k+1; i++){
        string kmer = s.substr(i, k);
        
        uint64_t hash[2];
        MurmurHash3_x64_128(kmer.c_str(), kmer.length(), seed, hash);
        if (static_cast<double>(hash[0])/(UINT64_MAX + 1.0) <= theta){
            kmers.insert(kmer);
        }
    }
    return kmers;
}

int intersect_size_sketch(const set<string>& set1, string t, int k, double theta, uint32_t seed){
    set<string> set2 = kspectrum_sketch(t, k,theta, seed);

    std::vector<string> intersect;

    std::set_intersection(set1.begin(), set1.end(),
                          set2.begin(), set2.end(),
                          std::back_inserter(intersect));

    //cout << "intsec: " << intersect.size() << endl;
    return intersect.size();
    
}
