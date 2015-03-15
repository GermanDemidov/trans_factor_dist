//
//  transcription_factor.h
//  trans_factors_distrib
//
//  Created by Герман Демидов on 11.03.15.
//  Copyright (c) 2015 gdemidov. All rights reserved.
//

#ifndef __trans_factors_distrib__transcription_factor__
#define __trans_factors_distrib__transcription_factor__

#include <stdio.h>
#include <vector>
#include <map>

class Transcription_factor {
public:
    Transcription_factor(int, const std::vector<std::map<char, double>>&);
    double calculate_weight_of_binding(std::string);
    void calculate_next_event();
    
private:
    int coordinate_in_sequence;
    bool binded_to_dna_non_specifically;
    bool binded_to_dna_specifically;
    const int id_of_current_tf;
    const std::vector<std::map<char, double>>& motif;
};

#endif /* defined(__trans_factors_distrib__transcription_factor__) */
