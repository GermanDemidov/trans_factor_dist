//
//  transcription_factor.cpp
//  trans_factors_distrib
//
//  Created by Герман Демидов on 11.03.15.
//  Copyright (c) 2015 gdemidov. All rights reserved.
//

#include "transcription_factor.h"

Transcription_factor::Transcription_factor(int id_of_tf, const std::vector<std::map<char, double>>& motif_of_tf): motif(motif_of_tf), id_of_current_tf(id_of_tf) {
    binded_to_dna_specifically = false;
    binded_to_dna_non_specifically = false;
    coordinate_in_sequence = 0;
    
}