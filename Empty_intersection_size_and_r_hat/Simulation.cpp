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

unordered_map<string, set<int>> occurrence;
mpreal E_upper, E_lower;
string type_map[] = {"random", "stress_random_lenk", "stress_ACGT", "stress_AC", "stress_A", "sickcase"};

class SimParams{
    public:
        int k;
        unsigned long L;
        double r;

    SimParams(){
    }

    SimParams(int k, unsigned long L, double r){
        this->k = k;
        this->L = L;
        this->r = r;
    }
    double get_one_minus_q(){
        return pow(1.0-r, k);
        //(1 − r)^k
    }
    double get_q(){
        return 1 - pow(1.0-r, k);
        //1 − (1 − r)^k
    }
    double get_p(){
        return r/3.0/(1-r);
    }

    unsigned long get_n(){
        return L+k-1;
    }  
};
set<string> kspectrum_update(string s, int k){
    set<string> kmers;
    for(int i=0; i<s.length()-k+1; i++){
        string kmer = s.substr(i, k);
        kmers.insert(kmer);
        occurrence[kmer].insert(i);
    }
    return kmers;
}

set<string> kspectrum(string s, int k){
    set<string> kmers;
    for(int i=0; i<s.length()-k+1; i++){
        string kmer = s.substr(i, k);
        kmers.insert(kmer);
    }
    return kmers;
}

int hammingDist(string str1, string str2) 
{ 
    int i = 0, count = 0; 
    while (str1[i] != '\0') { 
        if (str1[i] != str2[i]) 
            count++; 
        i++; 
    } 
    return count; 
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

pair<mpreal,pair<mpreal,mpreal>> cal_fun_g(const SimParams& sim, const string& seq, const string& tau, const set<int>& tau_occ){
    // // // PHASE1: INIT
    int k = sim.k;
    mpreal r = sim.r;
    mpreal one_minus_r = 1 - r;        
    unsigned long L = sim.L;
    mpreal p = r / (3 * (1 - r));
    mpreal rr = pow(one_minus_r, k);
    int copy = tau_occ.size();
    int HD = 0;
    
    int num_cluster = 0;
    int cur_pos_tau = -k;
    mpreal union_bound = 0;
    for (const auto& occ: tau_occ){
        if (occ >= cur_pos_tau +k){
            num_cluster++;
            cur_pos_tau = occ;
        }
    }

    for (int i = 0; i < L; ++i) {
        if (tau_occ.find(i) == tau_occ.end()){ // if not tau
            string s_i = (seq.substr(i, k));
            HD = hammingDist(s_i, tau);
            union_bound += rr * pow(p, HD);
        }
    }
    mpreal fun_f_tau = (1 - pow((1-rr),copy));
    mpreal theoretical_err_tau_upper = std::min(union_bound,(pow((1-rr),num_cluster)));
    mpreal theoretical_err_tau_lower = (1 - pow((1-rr),num_cluster) - (1 - pow((1-rr),copy)));
    pair<mpreal,mpreal> tau_error = {theoretical_err_tau_lower, theoretical_err_tau_upper};
    pair<mpreal,pair<mpreal,mpreal>> est_tau_error = {fun_f_tau,tau_error};
    return est_tau_error;
}



pair<double,pair<double,double>> cal_fun_f (const SimParams& sim, const set<string>& ksp_set, const string& seq, const unordered_map<string, set<int>>& occurrence){
    
    mpreal fun_f = 0;
    mpreal theoretical_err_lower = 0;
    mpreal theoretical_err_upper = 0;
    for (auto tau : ksp_set){
        pair<mpreal,pair<mpreal,mpreal>> est_tau_error = cal_fun_g(sim, seq, tau, occurrence.at(tau));
        fun_f += est_tau_error.first;
        theoretical_err_lower += est_tau_error.second.first;
        theoretical_err_upper += est_tau_error.second.second;
    } 

    pair<double,double> error = {theoretical_err_lower.toDouble(), theoretical_err_upper.toDouble()};
    pair<double,pair<double,double>> est_error = {fun_f.toDouble(), error};
    return est_error;
}

class ResultAggregator{
    public:
        vector<int> isizes;

        double isize_mean;
        double isize_sd;

        vector<double> test_vals_xtau;
        double test_vals_xtau_mean;
        double test_vals_xtau_sd;

        //vector<string> labels = {"r", "L","k", "num_reps", "estimate", "mean", "sd", "variance", "abs_error", "rel_error", "fixed_tau"};
        // = {"num_reps", "estimate", "mean", "sd", "variance", "abs_error", "rel_error"};

        string labels[9] = {"r", "L","k", "num_reps", "estimate", "mean", "sd", "variance", "fixed_tau"};
        vector<string> values;
        
        template<typename T>
        void put_values(T v){
            values.push_back(to_string(v));
        }

        double relative_error(double estimate, double mean){
            return abs((estimate-mean)/mean);
        }

        void diffReport(int num_replicates, double estimate, double mean, double sd, double r, int L, int k, string fixed_tau){
            put_values(r);
            put_values(L);
            put_values(k);
            put_values(num_replicates);
            put_values(estimate);
            put_values(mean);
            put_values(sd);
            put_values(sd*sd);
            // put_values(abs(estimate-mean));
            // put_values(relative_error(estimate, mean));
            values.push_back(fixed_tau);

            for(int i = 0; i< 9; i++){
                std::cout << std::fixed << std::setprecision(2);
                cout<<labels[i]<<" "<<values[i]<<" ";
            }
            cout<<endl;

            /*
            for(int i = 0; i< labels.size(); i++){
                cout<<values[i]<<" ";
            }
            cout<<endl;
            */
            
            //cout<<"num_reps "<< num_replicates<< " estimate "<<estimate<<" mean "<<mean<<" sd "<<sd<<" variance "<<sd*sd<<" abs_error "<<abs(estimate-mean)<<" rel_error "<<relative_error(estimate, mean)<<endl;
        }

        double calculateMean(const std::vector<int>& numbers) {
            double sum = 0.0;
            for (const auto& num : numbers) {
                sum += num;
            }
            return sum / numbers.size();
        }

        double calculateMean(const std::vector<double>& numbers) {
            double sum = 0.0;
            for (const auto& num : numbers) {
                sum += num;
            }
            return sum / numbers.size();
        }

        double calculateStandardDeviation(const std::vector<int>& numbers) {
            double mean = calculateMean(numbers);
            double squaredDiffSum = 0.0;
            for (const auto& num : numbers) {
                double diff = num - mean;
                squaredDiffSum += diff * diff;
            }
            double variance = squaredDiffSum / numbers.size();
            return std::sqrt(variance);
        }

        double calculateStandardDeviation(const std::vector<double>& numbers) {
            double mean = calculateMean(numbers);
            double squaredDiffSum = 0.0;
            for (const auto& num : numbers) {
                double diff = num - mean;
                squaredDiffSum += diff * diff;
            }
            double variance = squaredDiffSum / numbers.size();
            return std::sqrt(variance);
        }

        void calculateMeanSD(const std::vector<int>& numbers, double& mean, double& sd) {
            mean = calculateMean(numbers);
            double squaredDiffSum = 0.0;
            for (const auto& num : numbers) {
                double diff = num - mean;
                squaredDiffSum += diff * diff;
            }
            double variance = squaredDiffSum / numbers.size();
            sd = std::sqrt(variance);
        }

        void calculateMeanSD(const std::vector<double>& numbers, double& mean, double& sd) {
            mean = calculateMean(numbers);
            double squaredDiffSum = 0.0;
            for (const auto& num : numbers) {
                double diff = num - mean;
                squaredDiffSum += diff * diff;
            }
            double variance = squaredDiffSum / numbers.size();
            sd = std::sqrt(variance);
        }
};

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
    SimParams sim(k, L, r);
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
    ResultAggregator res;
    mpreal r_mp = sim.r;
    mpreal one_minus_r = 1 - r_mp;        
    mpreal rr = pow(one_minus_r, k);

    set<string> kmer_set = kspectrum_update(input_String, k);
    // pair<double,pair<double,double>> est_error = cal_fun_f(sim, kmer_set, input_String, occurrence);
    // double estimator = est_error.first;
    // double theoretical_err_lower = est_error.second.first;
    // double theoretical_err_upper = est_error.second.second;
    // E_upper = estimator + theoretical_err_upper;
    // E_lower = estimator + theoretical_err_lower;
    std::unordered_map<int, int> coeffs_Map; // a * x^b, b->a
    int max_copy = 0;
    for (const auto& pair : occurrence){
        int copy = pair.second.size();
        if (copy > max_copy){
            max_copy = copy;
        }
        coeffs_Map[copy]++;
    }



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
    for(int i=0; i<num_replicates; i++){
        std::string mutatedString = generateMutatedString(input_String, r);
        int I_obs = intersect_size(kmer_set, mutatedString, k);
        res.isizes.push_back(I_obs);
        mpreal q_hat = Newton_Method (coeffs_Map, kmer_set, k, epsilon, I_obs);
        mpreal r_hat = 1.0 - pow(1.0 - q_hat, 1.0/k);
        cout << r_hat << endl; 
    }
    
    
}