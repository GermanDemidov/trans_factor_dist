//
//  main.cpp
//  trans_factors_distrib
//
//  Created by Герман Демидов on 10.02.15.
//  Copyright (c) 2015 gdemidov. All rights reserved.
//

#include <iostream>
#include "parser.h"

class Solver {
    std::map<char, double> background_probabilities;
public:
    Solver(std::string fasta_filename, std::string pcm_filename) {
        // initialization of background frequencies
        background_probabilities.insert(std::make_pair('A', 0.28768819278776));
        background_probabilities.insert(std::make_pair('C', 0.21231180721224));
        background_probabilities.insert(std::make_pair('G', 0.21231180721224));
        background_probabilities.insert(std::make_pair('T', 0.28768819278776));

        
        // parse FASTA-file, creates hash map "sequences" inside fasta_sequences
        // key - fasta id, value - fasta sequence
        Parser_fasta fasta_sequences(fasta_filename);
        // fasta_sequences.debug_print();
        
        // parse PCM-file, creates PWM - hash map
        // key - TF name, value - vector of doubles
        Parser_pcm pcm_for_TFs(pcm_filename);
        
        // test of pcm_for_TFs
        std::vector<std::string> protein_names = pcm_for_TFs.return_protein_names();
        std::copy(protein_names.begin(), protein_names.end(), std::ostream_iterator<std::string>(std::cout, " "));
        std::vector<std::map<char, double>> profile_for_hkb = pcm_for_TFs.return_profile_for_protein(protein_names[protein_names.size() - 1]);
        std::cout << profile_for_hkb[1].at('C') << "\n"; // 9.326495726495676
        std::cout << profile_for_hkb[9].at('T') << "\n"; // 34.74119658119628
        
        
    }
};

int main(int argc, const char * argv[]) {
    // insert code here...
    std::string pcm_file_input = argv[1];
    std::string fasta_file_input = argv[2];
    std::cout << "Starting with PCM " << pcm_file_input << "\n" << "and input fasta " << fasta_file_input << "...\n";
    
    Solver solve(fasta_file_input, pcm_file_input);
    
    
    
    return 0;
}
