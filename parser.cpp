//
//  parser.cpp
//  trans_factors_distrib
//
//  Created by Герман Демидов on 10.03.15.
//  Copyright (c) 2015 gdemidov. All rights reserved.
//

#include "parser.h"
#include <math.h>


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
                std::transform(line.begin(), line.end(),line.begin(), ::toupper);
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

std::string& Parser_fasta::get_sequence_by_name(std::string name) {
    std::string founded_name;
    std::string tmp_best_candidate = "";
    for(std::map<std::string, std::string>::iterator iter = sequences.begin(); iter != sequences.end(); ++iter)
    {
        std::string current_name =  iter->first;
        std::string tmp_candidate = "";
        for (int i = 0; i < name.size(); i++) {
            if (name[i] == current_name[i]) {
                tmp_candidate += name[i];
                if (tmp_candidate.size() > tmp_best_candidate.size()) {
                    tmp_best_candidate = tmp_candidate;
                    founded_name = current_name;
                }
            }
        }
    }
    return sequences[founded_name];
}

std::map<std::string, std::string> Parser_fasta::return_sequences() {
    return sequences;
}

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
            } else if (line[0] != '<') {
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

std::vector<std::map<char, double>> Parser_pcm::return_profile_for_protein_pcm(std::string tf_name) {
    return pcm.at(tf_name);
}

std::vector<std::map<char, double>> Parser_pcm::return_profile_for_protein_pwm(std::string tf_name) {
    return pwm.at(tf_name);
}


// PWM calculation section
// formulas used: 
void Parser_pcm::calculate_pwm(std::map<char, double>& background_probabilites) {
    for (std::map<std::string, std::vector<std::map<char, double>>>::iterator it = pcm.begin(); it != pcm.end(); ++it) {
        std::vector<std::map<char, double>> tmp_vector_for_pwm;
        
        double W = 0.0;
        double a = 0.0;

        for (std::vector<char>::iterator nucleot = nucleotides.begin(); nucleot != nucleotides.end(); ++nucleot) {
                W += it->second[0].at(*nucleot);
            }
        a = log(W);
        for (std::vector<std::map<char, double>>::iterator it_v = it->second.begin(); it_v != it->second.end(); ++it_v) {
            std::map<char, double> tmp_map_for_pwm;
            for (std::vector<char>::iterator nucleot = nucleotides.begin(); nucleot != nucleotides.end(); ++nucleot) {
                double value = log(
                                   (it_v->at(*nucleot) + a * background_probabilites.at(*nucleot)) /
                                   ((W + a) * background_probabilites.at(*nucleot))
                                   );
                // debug print
                // uncomment if needed
                // std::cout << it->first << " " << *nucleot << " " << it_v->at(*nucleot) << " " << a << " " <<background_probabilites.at(*nucleot) << " " << W << "\n" << value << "\n\n";
                tmp_map_for_pwm.insert(std::make_pair(*nucleot, value));
            }
            tmp_vector_for_pwm.push_back(tmp_map_for_pwm);
        }
        pwm.insert(std::make_pair(it->first, tmp_vector_for_pwm));
    }
}




Parser_dnase_acc::Parser_dnase_acc(std::string input_file) {
    std::ifstream pcm_file (input_file);
    int start = -1;
    int end = -1;
    short current_state = -1;
    short previous_state = -1;
    int counter = -1;
    std::cout << input_file << "\n";
    if (pcm_file.is_open()) {
        while (pcm_file >> current_state) {
            counter++;
            if (previous_state == 0 && current_state == 1) {
                start = counter;
                no_acc_intervals.push_back(std::make_pair(end, start));
            }
            if (previous_state == 1 && current_state == 0) {
                end = counter;
                open_acc_intervals.push_back(std::make_pair(start, end));
            }
            previous_state = current_state;
        }
    }

    if (current_state == 0) no_acc_intervals.push_back(std::make_pair(end, counter));
    else no_acc_intervals.push_back(std::make_pair(start, counter));
};

bool Parser_dnase_acc::is_in_interval(std::pair<int, int> coords_of_tf_binding) {
    for (auto& coords : open_acc_intervals) {
        if (coords.first <= coords_of_tf_binding.first) {
            if (coords.second >= coords_of_tf_binding.second) {
                return true;
            }
        }
    }
    return false;
}

bool Parser_dnase_acc::is_in_interval(int x) {
    for (auto& coords : open_acc_intervals) {
        if (coords.first <= x && coords.second >= x) {
                return true;
        }
    }
    return false;
}

std::vector<std::pair<int, int>> Parser_dnase_acc::return_no_acc_intervals() {
    return no_acc_intervals;
}




