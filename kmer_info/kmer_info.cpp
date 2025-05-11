#include <iostream>
#include <filesystem>
#include <random>
#include <string>
#include <set>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <map>
#include <unordered_map>
using namespace std;

namespace fs = std::filesystem;

bool ensure_directory_exists(const std::string& dir_path) {
    try {
        if (!fs::exists(dir_path)) {
            bool created = fs::create_directories(dir_path);
            
            if (created) {
                std::cout << "Directory created: " << dir_path << std::endl;
                return true;
            } else {
                std::cerr << "Failed to create directory: " << dir_path << std::endl;
                return false;
            }
        } else if (!fs::is_directory(dir_path)) {
            std::cerr << "Path exists but is not a directory: " << dir_path << std::endl;
            return false;
        } else {
            std::cout << "Directory already exists: " << dir_path << std::endl;
            return true;
        }
    }
    catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
        return false;
    }
}

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

int hammingDist(const string& str1, const string& str2) 
{ 
    int i = 0, count = 0; 
    while (str1[i] != '\0') { 
        if (str1[i] != str2[i]) 
            count++; 
        i++; 
    } 
    return count; 
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
    srand(time(nullptr));
    string filename="";
    string outfolder="";
    int L = 0;
    int k = 0;
    int start_pos = 0;


    vector<string> args(argv + 1, argv + argc);
    for (auto i = args.begin(); i != args.end(); ++i) {
        if (*i == "-h" || *i == "--help") {
            cout << "Syntax: -i [input-file] -l <length>[100] -s <start position>[5]  -k <kmer-size>[33] -o <output folder>" << endl;
            return 0;
        } else if (*i == "-l") {
            L = stoull(*++i);
        } else if (*i == "-k") {
            k = stoi(*++i);
        } else if (*i == "-i") {
            filename = (*++i);
        } else if (*i == "-s") {
            start_pos = stoi(*++i);
        } else if (*i == "-o") {
            outfolder = (*++i);
        }
    }

    std::string input_String;
    input_String = readSubstring(filename, start_pos ,L + k - 1);
    //cout << "input" << input_String << endl;
    set<string> kmer_set = kspectrum_update(input_String, k);

    std::map<int, int> coeffs_Map; // a * x^b, b->a
    for (const auto& pair : occurrence){
        int copy = pair.second.size();
        coeffs_Map[copy]++;
    }
    

    std::map<int, int> diff_Map; // occ - sep : counts
    for (const auto& tau : occurrence){
        int num_cluster = 0;
        int cur_pos_tau = -k;
        int occ_count = 0;
        for (const auto& occ: tau.second){
            occ_count ++;
            if (occ >= cur_pos_tau +k){
                num_cluster++;
                cur_pos_tau = occ;
            }
        }
        diff_Map[occ_count - num_cluster] ++;
    }

    std::map<int, int> HD_Map;
    std::vector<std::string> kmer_vec(kmer_set.begin(), kmer_set.end());

    for (size_t i = 0; i < kmer_vec.size(); ++i) {
        const auto& tau = kmer_vec[i];
        for (size_t j = i + 1; j < kmer_vec.size(); ++j) {
            const auto& upsilon = kmer_vec[j];
            int HD = hammingDist(tau, upsilon);
            HD_Map[HD]++;
        }
    }


    if (ensure_directory_exists("./"+ outfolder)){
        std::ofstream occfile("./"+ outfolder + "/occ_counts.txt");
        for (const auto& pair : coeffs_Map){
            occfile << pair.first << " : " << pair.second << std::endl;
        }
        occfile.close();
    
        std::ofstream difffile("./"+ outfolder + "/diff_counts.txt");
        for (const auto& pair : diff_Map){
            difffile << pair.first << " : " << pair.second << std::endl;
        }
        difffile.close();

        std::ofstream HDfile("./"+ outfolder + "/HD_counts.txt");
        for (const auto& pair : HD_Map){
            HDfile << pair.first << " : " << pair.second << std::endl;
        }
        HDfile.close();
    }
    
    


    //cout<< filename << ":" << start_pos << " " << L <<" "<< k << " " << num_replicates << " " << q_hat.toDouble() << " " << r_hat.toDouble() << " " << r_prime << endl;
}