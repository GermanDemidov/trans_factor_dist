//
//  parser.cpp
//  trans_factors_distrib
//
//  Created by Герман Демидов on 10.03.15.
//  Copyright (c) 2015 gdemidov. All rights reserved.
//

#include "parser.h"


// FASTA parser section
Parser_fasta::Parser_fasta(std::string filename) {
    std::ifstream fasta_file (filename);
    if (fasta_file.is_open())
    {
        std::string line;
        std::string current_seq_id = "";
        std::string current_seq = "";
        while (getline (fasta_file, line) )
        {
            if (line[0] == '>') {
                if (current_seq_id != "") {
                    sequences.insert(make_pair(current_seq_id, current_seq));
                    current_seq = "";
                }
                current_seq_id = line;
            } else {
                current_seq += line;
            }
        }
        if (current_seq_id != "") {
            sequences.insert(make_pair(current_seq_id, current_seq));
        }
        fasta_file.close();
    } else {
        std::cout << "Please, check your fasta file.\n";
    }
};

void Parser_fasta::debug_print() {
    for (std::map<std::string, std::string>::iterator it = sequences.begin(); it != sequences.end(); ++it) {
        std::cout << it->first << "\n" << it->second << "\n";
    }
}



// PCM parser section
Parser_pcm::Parser_pcm(std::string filename) {
    nucleotides.push_back('A');
    nucleotides.push_back('C');
    nucleotides.push_back('G');
    nucleotides.push_back('T');

    std::ifstream pcm_file (filename);
    if (pcm_file.is_open()) {
        std::string line;
        std::string id_of_tf = "";
        std::vector<std::map<char, double>> current_profile;
        while (getline(pcm_file, line)) {
            std::stringstream ss(line);
            if (line[0] == '>') {
                if (id_of_tf != "") {
                    pcm.insert(make_pair(id_of_tf, current_profile));
                    current_profile.clear();
                }
                id_of_tf = "";
                int length_of_profile = 0;
                ss >> id_of_tf;
                ss >> length_of_profile;
                protein_names.push_back(id_of_tf);
            } else {
                std::map<char, double> current_profile_for_position;
                double current_probability = 0.0;
                for (std::vector<char>::iterator it = nucleotides.begin(); it != nucleotides.end(); ++it) {
                    ss >> current_probability;
                    current_profile_for_position.insert(std::make_pair(*it, current_probability));
                }
                current_profile.push_back(current_profile_for_position);
            }
        }
        if (id_of_tf != "") {
            pcm.insert(make_pair(id_of_tf, current_profile));
        }
        pcm_file.close();
    } else {
        std::cout << "Please, check your PCM file.\n";
    }
}

std::vector<std::string> Parser_pcm::return_protein_names() {
    return protein_names;
}

std::vector<std::map<char, double>> Parser_pcm::return_profile_for_protein(std::string tf_name) {
    return pcm.at(tf_name);
}


// PWM calculation section
// formulas used: 
void Parser_pcm::calculate_pwm(std::map<char, double> background_probabilites) {
    
}





