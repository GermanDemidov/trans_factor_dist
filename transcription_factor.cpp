//
//  transcription_factor.cpp
//  trans_factors_distrib
//
//  Created by Герман Демидов on 11.03.15.
//  Copyright (c) 2015 gdemidov. All rights reserved.
//

#include "transcription_factor.h"


Transcription_factor::Transcription_factor(int id_of_tf, std::vector<std::map<char, double>>motif_of_tf, std::string type_of_tf): motif(motif_of_tf), id_of_current_tf(id_of_tf), type_name(type_of_tf) {
    binded_to_dna_specifically = false;
    binded_to_dna = false;
    coordinate_in_sequence = 0;
    binded_to_forward_dna = false;
    
}

double Transcription_factor::calculate_weight_of_binding(std::string& DNA) {
    double weight_in_current_pos = 0.0;
    
    //std::cout << motif[0].at('A') << " MOTIF \n";
    for (int i = 0; i < motif.size(); ++i) {
        weight_in_current_pos += motif[i].at(DNA.at(coordinate_in_sequence + i));
    }
    
    return weight_in_current_pos;
}

void Transcription_factor::change_coordinate_in_sequence(int i) {
    coordinate_in_sequence = i;
}

int Transcription_factor::get_size() {
    return (int)motif.size();
}

bool Transcription_factor::is_binded() {
    return binded_to_dna;
}

void Transcription_factor::bind_to_dna(bool bind_to_forward) {
    binded_to_dna = true;
    if (bind_to_forward)
        binded_to_forward_dna = true;
    else binded_to_forward_dna = false;
}

void Transcription_factor::unbind_from_dna() {
    binded_to_dna = false;
}

bool Transcription_factor::is_binded_to_forward() {
    return binded_to_forward_dna;
}

int Transcription_factor::get_coordinate_in_sequence() {
    return coordinate_in_sequence;
}

void Transcription_factor::set_binded_to_dna_specifically() {
    binded_to_dna_specifically = true;
}
