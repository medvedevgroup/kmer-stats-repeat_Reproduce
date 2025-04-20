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
#include "sketch.h"
using namespace std;


#include "mpreal.h"
using mpfr::mpreal;
// Required precision of computations in decimal digits
const int digits = 50;  // Play with it to check different precisions
unordered_map<string, set<int>> occurrence;
set<string> kspectrum_update(string s, int k){
    set<string> kmers;
    for(int i=0; i<s.length()-k+1; i++){
        string kmer = s.substr(i, k);
        kmers.insert(kmer);
        occurrence[kmer].insert(i);
    }
    return kmers;
}

std::string generateRandomString(int length) {
    const std::string charset = "ACGT"; // Characters to choose from
    std::random_device rd;
    std::mt19937 gen(rd());

    //std::mt19937 gen(seed); to set seed
    std::uniform_int_distribution<> dis(0, charset.length() - 1);

    std::string randomString;
    randomString.reserve(length);

    for (int i = 0; i < length; ++i) {
        randomString += charset[dis(gen)];
    }

    return randomString;
}

string generateStressTestString(int length, int k, int type){
    std::string pattern = "ACGT";
    if(type == 1){
        pattern = generateRandomString(k); //say k=5, and randomstring is GTACTA: then GTACTAGTACTAGTACTA....
    }else if(type == 2){
        pattern = "ACGT";     //ACGTACGTACGT....
    }else if(type == 3){
        pattern = "AC";     //ACACACACACAC....
    }else if(type == 4){
        pattern = "A";   //AAAAAAAAAAAAA....
    }else if(type == 5){
        pattern = "A";   //AAAAAAAAAAAAA....AG
    }

    //cout << "Pattern: " << pattern << endl;
    std::string s = "";
    std::string result;
    for (int i = 0; i < length; ++i) {
        s += pattern[i % pattern.length()];
    }
    if (type==5){
        s += "G";
    }

    return s;
}

std::string generateMutatedString(string s, double r) {
    string notT = "ACG";
    string notG = "ACT";
    string notC = "AGT";
    string notA = "CGT";

    std::random_device rd;
    std::mt19937 gen(rd());
    
    std::string mutatedString;
    mutatedString.reserve(s.length());

    for (int i = 0; i < s.length(); ++i) {
        std::uniform_real_distribution<> dis1(0.0, 1.0);
        // std::uniform_real_distribution<> dis2(0.0, 1.0);
        // std::uniform_real_distribution<> dis3(0.0, 1.0);
        double random_value = dis1(gen);
        
        
        //cout<<i<<" "<<random_value<<" "<<(1-r)<<endl;
        
        if (random_value < (1-r)) {
            mutatedString += s[i]; // Same with probability (1-r) //stay same
        } else {
            std::uniform_int_distribution<> dis(0, 2); // Range: {0, 1, 2}
            if(s[i] == 'A'){
                mutatedString += notA[dis(gen)];
            }else if(s[i] == 'C'){
                mutatedString += notC[dis(gen)];
            }else if(s[i] == 'G'){
                mutatedString += notG[dis(gen)];
            }else if(s[i] == 'T'){
                mutatedString += notT[dis(gen)];
            }
            
        }
    }
    return mutatedString;
}

int intersect_size(const set<string>& set1, string t, int k){
    set<string> set2 = kspectrum(t, k);

    std::vector<string> intersect;

    std::set_intersection(set1.begin(), set1.end(),
                          set2.begin(), set2.end(),
                          std::back_inserter(intersect));
    return intersect.size();
}

string readSubstring(string filename, int start ,int n){
    ifstream file(filename); 
    string line;
    string sequence;

    while (getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') { 
            continue;
        } else {
            sequence += line; 
        }
    }

    file.close();
    string seq = sequence.substr(start,n);
    transform(seq.begin(), seq.end(), seq.begin(),::toupper);
    return seq;
}



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
    int skecth_replicates = 1;

    vector<string> args(argv + 1, argv + argc);
    for (auto i = args.begin(); i != args.end(); ++i) {
        if (*i == "-h" || *i == "--help") {
            cout << "Syntax: -i [input-file] -l <length>[100] -s <start position>[5]  -k <kmer-size>[33] -e <presion>[0.001] -c <replicates>[30] -r <mutation rate> 0.1 -t <type> 0 -z <skecth_replicates> 1" << endl;
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
        } else if (*i == "-z") {
            skecth_replicates = stoi(*++i);
        } 
    }

    std::string input_String;
    //input_String = readSubstring(filename, start_pos ,L + k - 1);

    if (filename==""){
        if(type == 0){
            input_String = generateRandomString(L + k - 1);
        }else{
            input_String = generateStressTestString(L + k - 1, k, type);
        }
    }else{
        input_String = readSubstring(filename, start_pos , L + k - 1);
    }
    //cout << "input" << input_String << endl;
    set<string> kmer_set = kspectrum_update(input_String, k);
    //set<string> kmer_set_sketch;
    // if (theta < 1.0){
    //     kmer_set_sketch = kspectrum_sketch(input_String, k, theta);
    // }

    std::unordered_map<int, int> coeffs_Map; // a * x^b, b->a
    for (const auto& pair : occurrence){
        int copy = pair.second.size();
        coeffs_Map[copy]++;
    }
    std::random_device rd;
    /*generate once mutation*/
    std::string mutatedString = generateMutatedString(input_String, r);
    int I_obs = intersect_size(kmer_set, mutatedString, k);
        
    //cout << "I_obs" << I_obs <<endl;
    mpreal q_hat = Newton_Method (coeffs_Map, kmer_set, k, epsilon, I_obs);
    //cout<< q_hat << endl;
    mpreal r_hat = 1.0 - pow(1.0 - q_hat, 1.0/k);

    for(int i=0; i<num_replicates; i++){
        int I_obs_sketch=0;

        vector<double> thresholds = {0.1, 0.01};
        for(int j=0; j< skecth_replicates ; j++){
            cout << r_hat;
            for (double theta: thresholds){
                uint32_t seed = rd();
                set<string> kmer_set_sketch = kspectrum_sketch(input_String, k, theta, seed);
                I_obs_sketch = round(static_cast<double>(intersect_size_sketch(kmer_set_sketch, mutatedString, k, theta, seed)) / (theta));
                mpreal q_hat_sketch = 0;
                mpreal r_hat_sketch = 0;
                if (I_obs_sketch < kmer_set.size() && I_obs_sketch > 0 ){
                    q_hat_sketch = Newton_Method (coeffs_Map, kmer_set, k, epsilon, I_obs_sketch);
                    r_hat_sketch = 1.0 - pow(1.0 - q_hat_sketch, 1.0/k);
                }
                else if (I_obs_sketch == 0){
                    q_hat_sketch = 1.0;
                    r_hat_sketch = 1.0;
                }
                cout << " , " << r_hat_sketch;
            }
            cout << endl;
        }
          
    }

    //cout<< filename << ":" << start_pos << " " << L <<" "<< k << " " << num_replicates << " " << q_hat.toDouble() << " " << r_hat.toDouble() << " " << r_prime << endl;
}