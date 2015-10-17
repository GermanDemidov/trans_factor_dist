//
//  Transcription_factors_in_cell.cpp
//  trans_factors_distrib
//
//  Created by Герман Демидов on 24.04.15.
//  Copyright (c) 2015 gdemidov. All rights reserved.
//

#include "transcription_factors_in_cell.h"

Transcription_factors_in_cell::Transcription_factors_in_cell(std::vector<std::string>& prot_names, std::vector<Transcription_factor>& tfs_to_start) {
    tfs = tfs_to_start;
    protein_names = prot_names;
    for (int i = 0; i < protein_names.size(); i++) {
        number_of_unbinded_TFs[protein_names[i]] = 0;
    }
    for (int i = 0; i < tfs.size(); ++i) {
        number_of_unbinded_TFs[tfs[i].type_name]++;
        indexes_of_unbinded_TFs.insert(i);
    }
    
    Choose_element_from_array<int> rand_elem_int;
    choose_for_int_vect = rand_elem_int;
}

int Transcription_factors_in_cell::choose_next_unbinded_DNA_to_interact() {
    int answer = -1;
    if (indexes_of_unbinded_TFs.size() > 0) {
        answer = pick_a_number(0, indexes_of_unbinded_TFs.size() - 1);
        std::set<int>::iterator it = indexes_of_unbinded_TFs.begin();
        for (; answer != 0; answer--) it++;
        return *it;
    }
    return answer;
}

int Transcription_factors_in_cell::choose_next_binded_DNA_to_interact() {
    int answer = -1;
    if (!indexes_of_binded_TFs.empty()) {
        answer = pick_a_number(0, indexes_of_binded_TFs.size() - 1);
        std::set<int>::iterator it = indexes_of_binded_TFs.begin();
        for (; answer != 0; answer--) it++;
        return *it;
    } else {
        std::cerr << "Binded to DNA list is empty and you try to interact with them" << "\n";
        return -1;
    }
    return *indexes_of_binded_TFs.begin();
}

double Transcription_factors_in_cell::generate_next_characteristic_time_for_binded_TFs() {
    return generate_next_time(200.0 * indexes_of_binded_TFs.size());
}

double Transcription_factors_in_cell::generate_next_characteristic_time_for_unbinded_TFs() {
    int sum_of_unbinded_TFs = 0;
    for (std::map<std::string, int>::iterator it = number_of_unbinded_TFs.begin(); it != number_of_unbinded_TFs.end(); ++it) {
        sum_of_unbinded_TFs += it->second;
    }
    return generate_next_time(sum_of_unbinded_TFs);
}

void Transcription_factors_in_cell::bind_tf_to_dna(int i, bool bind_to_forward, int coord_of_binding) {
    tfs[i].bind_to_dna(bind_to_forward, coord_of_binding);
    indexes_of_unbinded_TFs.erase(i);
    indexes_of_binded_TFs.insert(i);
    number_of_unbinded_TFs[tfs[i].type_name] -= 1;
}

void Transcription_factors_in_cell::unbind_tf_from_dna(int i) {
    tfs[i].unbind_from_dna();
    indexes_of_binded_TFs.erase(i);
    number_of_unbinded_TFs[tfs[i].type_name] += 1;
    indexes_of_unbinded_TFs.insert(i);
}

int Transcription_factors_in_cell::number_of_binded_dna() {
    return (int) indexes_of_binded_TFs.size();
}

std::string Transcription_factors_in_cell::get_type_of_TF(int i) {
    return tfs[i].type_name;
}

int Transcription_factors_in_cell::get_size_of_TF(int i) {
    return tfs[i].get_size();
}

Transcription_factor& Transcription_factors_in_cell::get_tf_by_index(int i) {
    return tfs[i];
}

int Transcription_factors_in_cell::num_of_binded_TFs() {
    return indexes_of_binded_TFs.size();
}

int Transcription_factors_in_cell::num_of_unbinded_TFs() {
    int answer = 0;
    for (int i = 0; i < protein_names.size(); ++i) {
        answer += protein_names[i].size();
    }
    return answer;
}


std::map<std::string, std::vector<double> > Transcription_factors_in_cell::return_average_sliding_lengths_vectors() {
    std::map<std::string, std::vector<double> > map_with_averages;
    for (int i = 0; i < tfs.size(); i++) {
        if (!tfs[i].sliding_length.empty()) {
            if (map_with_averages.find(tfs[i].type_name) == map_with_averages.end() ) {
                std::vector<double> tmp_vect;
                map_with_averages.insert(std::make_pair(tfs[i].type_name, tmp_vect));
            }
            double average = 0.0;
            for (int j = 0; j < tfs[i].sliding_length.size(); j++) {
                    average += tfs[i].sliding_length[j];
            }
            average /= tfs[i].sliding_length.size();
            if (average >= 1) {
                map_with_averages[tfs[i].type_name].push_back(average);
            }
        }
    }
    return map_with_averages;
}

void Transcription_factors_in_cell::null_everything(std::vector<std::string> prot_names) {
    protein_names = prot_names;
    for (int i = 0; i < protein_names.size(); i++) {
        number_of_unbinded_TFs[protein_names[i]] = 0;
    }
    indexes_of_unbinded_TFs.clear();
    for (int i = 0; i < tfs.size(); ++i) {
        number_of_unbinded_TFs[tfs[i].type_name]++;
        indexes_of_unbinded_TFs.insert(i);
    }
    indexes_of_binded_TFs.clear();
    for (int i = 0; i < tfs.size(); i++) {
        tfs[i].null_everything();
    }
    Choose_element_from_array<int> rand_elem_int;
    choose_for_int_vect = rand_elem_int;
}

