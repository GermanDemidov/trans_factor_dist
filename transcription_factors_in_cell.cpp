//
//  Transcription_factors_in_cell.cpp
//  trans_factors_distrib
//
//  Created by Герман Демидов on 24.04.15.
//  Copyright (c) 2015 gdemidov. All rights reserved.
//

#include "Transcription_factors_in_cell.h"

Transcription_factors_in_cell::Transcription_factors_in_cell(std::vector<std::string> prot_names, std::vector<Transcription_factor>& tfs_to_start) {
    tfs = tfs_to_start;
    protein_names = prot_names;
    for (int i = 0; i < protein_names.size(); i++) {
        number_of_unbinded_TFs[protein_names[i]] = 0;
    }
    for (int i = 0; i < tfs.size(); ++i) {
        number_of_unbinded_TFs[tfs[i].type_name]++;
    }
    
    Choose_element_from_array<int> rand_elem_int;
    choose_for_int_vect = rand_elem_int;
}

int Transcription_factors_in_cell::choose_next_unbinded_DNA_to_interact() {
    int answer = -1;
    bool succesfully_found_unbinded = false;
    while (!succesfully_found_unbinded) {
        answer = pick_a_number(0, tfs.size() - 1);
        if (!tfs[answer].is_binded()) succesfully_found_unbinded = true;
    }
    return answer;
}

int Transcription_factors_in_cell::choose_next_binded_DNA_to_interact() {
    int answer = -1;
    if (!indexes_of_binded_TFs.empty()) {
        answer = pick_a_number(0, indexes_of_binded_TFs.size() - 1);
        std::set<int>::iterator it = indexes_of_binded_TFs.begin();
        for (; answer != 0; answer--) it++;
        std::cout << "IT" << *it << "\n";
        return *it;
    } else {
        std::cerr << "Binded to DNA list is empty and you try to interact with them" << "\n";
        return -1;
    }
    return *indexes_of_binded_TFs.begin();
}

double Transcription_factors_in_cell::generate_next_characteristic_time_for_binded_TFs() {
    return generate_next_time(0.1 * indexes_of_binded_TFs.size());
}

double Transcription_factors_in_cell::generate_next_characteristic_time_for_unbinded_TFs() {
    int sum_of_unbinded_TFs = 0;
    for (std::map<std::string, int>::iterator it = number_of_unbinded_TFs.begin(); it != number_of_unbinded_TFs.end(); ++it) {
        sum_of_unbinded_TFs += it->second;
    }
    return generate_next_time(sum_of_unbinded_TFs);
}

void Transcription_factors_in_cell::bind_tf_to_dna(int i, bool bind_to_forward) {
    tfs[i].bind_to_dna(bind_to_forward);
    indexes_of_binded_TFs.insert(i);
    std::cout << indexes_of_binded_TFs.size() << " BINDED SUCCESFULLY\n";
    number_of_unbinded_TFs[tfs[i].type_name] -= 1;
}

void Transcription_factors_in_cell::unbind_tf_from_dna(int i) {
    tfs[i].unbind_from_dna();
    indexes_of_binded_TFs.erase(i);
    number_of_unbinded_TFs[tfs[i].type_name] += 1;
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


