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
#include <string>
#include <iostream>

class Transcription_factor {
public:
    std::string type_name;
    Transcription_factor(int, std::vector<std::map<char, double>>, std::string);
    double calculate_weight_of_binding(std::string&);
    void calculate_next_event();
    void change_coordinate_in_sequence(int);
    int get_size();
    bool is_binded();
    int get_coordinate_in_sequence();
    
    void set_binded_to_dna_specifically();
    
    // change bool field of binding
    void bind_to_dna(bool, int);
    void unbind_from_dna();
    
    bool is_binded_to_forward();
    bool is_binded_specifically();
    std::vector<int> sliding_length;
    double time_of_binding = 0.0;
    double time_of_unbinding = 0.0;
    
    void null_everything();


private:
    int start_coord = 0;
    int finish_coord = 0;
    bool binded_to_forward_dna = false;
    int coordinate_in_sequence;
    bool binded_to_dna;
    bool binded_to_dna_specifically;
    int id_of_current_tf;
    std::vector<std::map<char, double>> motif;

};

#endif /* defined(__trans_factors_distrib__transcription_factor__) */
