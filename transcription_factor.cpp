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
    binded_to_dna_non_specifically = false;
    coordinate_in_sequence = 0;
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
